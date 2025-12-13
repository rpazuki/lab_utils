##############################################################
#  labUtils: Metabolic network mapping utilities             #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
#                                                            #
##############################################################

import difflib
import logging
import warnings
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd

from labUtils.utils import (
    find_molecular_weight,
    find_molecular_weight_by_id,
    get_compound_by_id,
    get_compound_by_name,
    ions_from_smiles,
    mg_to_mmol,
    od600_to_gCDW,
)

# pylint: disable=import-outside-toplevel


# and enumerator for medium source
class MediumSource(Enum):
    UNSTATED = "unstated"
    SUPPLEMENT = "supplement"
    MEDIUM = "medium"
    FIXED = "fixed"


def build_mappings(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: dict[str, str],
    supplement_column: str,
    custom_mapping_file: str | Path | None = None,
    custom_mapping: dict[str, str | dict[str, Any]] | None = None,
    fixed_exchanges: list[str] | None = None,
    separator: str = ";",
    fuzzy_threshold: float = 0.6,
    verbose: bool = False,
) -> pd.DataFrame:
    ###########################################################################
    ## Level 1: Base mapping from supplement_to_exchange_map
    # Create a dataframe from organism's supplement to exchange mapping
    supp_exchange = (
        (
            supp.strip().lower(),
            "",
            exchange,
            0.0,  # mass_per_litre_mg
            0.0,  # mmol_concentration
            0.0,  # flux_upper_bound
            MediumSource.UNSTATED.value,
        )
        for supp, exchange in supplement_to_exchange_map.items()
        if supp is not None
    )

    df_mapping = pd.DataFrame(
        list(supp_exchange),
        columns=[
            "name",
            "iupac_name",
            "exchange_reaction",
            "mass_per_litre",
            "mmol_concentration",
            "flux_upper_bound",
            "source",
        ],
    )
    #############################################################################
    # Collect unique supplements from the dataframe
    unique_supplements = set()
    if supplement_column in growth_rates_df.columns:
        for val in growth_rates_df[supplement_column].dropna().astype(str):
            parts = [s.strip() for s in val.split(separator) if s.strip()]
            for p in parts:
                unique_supplements.add(p.lower())

    def find_in_df_mapping(df_mapping, name: str) -> tuple[pd.DataFrame, bool]:
        df = df_mapping.loc[df_mapping["name"] == name, :]
        if not df.empty:
            return df, True
        df = df_mapping.loc[df_mapping["iupac_name"] == name, :]
        if not df.empty:
            return df, False
        return pd.DataFrame(), False

    def update_mapping_df(df_mapping, row):
        # Update or append to df_mapping
        existing_idx, found_by_name = find_in_df_mapping(df_mapping, row["name"])
        if not existing_idx.empty:
            if not found_by_name:
                original_name = existing_idx["name"].values[0]
                iupac_name = row["name"]
                row["name"] = original_name
                row["iupac_name"] = iupac_name
            df_mapping.loc[existing_idx.index[0]] = row
            #     print("=" * 50)
            #     print(f" Updated mapping by name for '{row}'.")
            #     print("=" * 50)
            # elif row["iupac_name"] != "":
            existing_idx, found_by_name = find_in_df_mapping(df_mapping, row["iupac_name"])
            if not existing_idx.empty:
                if found_by_name:
                    original_name = existing_idx["name"].values[0]
                    row["name"] = original_name
                df_mapping.loc[existing_idx.index[0]] = row
                # print("=" * 50)
                # print(f" Updated mapping by iupac name for '{row}'.")
                # print("=" * 50)
            else:
                df_mapping = pd.concat([df_mapping, pd.DataFrame([row])], ignore_index=True)
                # print("=" * 50)
                # print(f" Add mapping for '{row}'.")
                # print("=" * 50)
        else:
            df_mapping = pd.concat([df_mapping, pd.DataFrame([row])], ignore_index=True)
            # print("=" * 50)
            # print(f" Add mapping for '{row}'.")
            # print("=" * 50)
        return df_mapping

    def create_new_rows_for_salts(df_mapping, salt_name, properties):
        rows = []
        salt_name = salt_name.strip().lower()
        if properties.get("pubchem_name", None):
            compound = get_compound_by_name(properties["pubchem_name"])  # type: ignore
        elif properties.get("pubchem_id", None):
            compound = get_compound_by_id(properties["pubchem_id"])  # type: ignore
        elif salt_name:
            compound = get_compound_by_name(salt_name)
        else:
            raise ValueError(f"The exchange name:'{salt_name}' is marked as salt but no valid name or ID was provided.")
        compound_smiles: str = compound.smiles  # type: ignore
        #
        compounds = ions_from_smiles(compound_smiles)
        if verbose:
            print(f"\nProcessing salt '{salt_name}' with SMILES:'{compound_smiles}' into ions:")
            for ion in compounds:
                print(f"  Ion: {ion}")

        # assert len(compounds) == len(properties.get("exchange_name", [])), "Mismatch between ions and exchange names"
        total_molar_mass = sum(ion["molar_mass"] for ion in compounds)
        for ion_info in compounds:
            # ion_name = ion_info["smiles_component"]
            # there are two iupac names obtained from PubChem: iupac_name & iupac_name_secondary
            # We first lookup the exsiting names in the df_mapping. If any one of them exist, use it.
            # Otherwise, we use the iupac_name from ion_info
            iupac_name = ion_info["iupac_name"]
            df, _ = find_in_df_mapping(df_mapping, iupac_name)
            if not df.empty:
                ion_name = iupac_name
                ion_exchange = df["exchange_reaction"].values[0]
            else:
                iupac_name_secondary = ion_info.get("iupac_name_secondary", "Unknown")
                df, _ = find_in_df_mapping(df_mapping, iupac_name_secondary)
                if not df.empty:
                    ion_name = iupac_name_secondary
                    ion_exchange = df["exchange_reaction"].values[0]
                else:
                    # ion_name = ion_info["smiles_component"]
                    raise ValueError(
                        f"Could not find ion name:'{iupac_name}' or "
                        f"its iupac name '{iupac_name_secondary}' for salt '{salt_name}' in mapping dataframe."
                    )

            ratio = ion_info["molar_mass"] / total_molar_mass if total_molar_mass > 0 else 0.0
            pubchem_cid = ion_info.get("pubchem_cid", None)
            # copy the properties to avoid modifying the original
            row_properties = properties.copy()
            # Adjust mass_per_litre based on ratio
            mass_per_litre = properties.get("mass_per_litre", 0.0)
            row_properties["mass_per_litre"] = float(mass_per_litre) * ratio  # type: ignore
            row_properties["name"] = ion_name.strip().lower()
            # row_properties["exchange_name"] = ion_exchange
            row_properties["exchange_name"] = ion_exchange
            if pubchem_cid is not None:
                row_properties["pubchem_id"] = pubchem_cid
            row = create_new_row(ion_name, row_properties)
            rows.append(row)

        if verbose:
            print("=" * 100)
            print(f" Created mapping for salt '{salt_name}'")
            for r in rows:
                print(f"   ion '{r['name']}': {r}")
            print("=" * 100)
        return rows

    def create_new_row(supp, properties):
        row = {"name": supp.strip().lower()}
        # Extract exchange_name from nested dict
        exchange_name: str = properties.get("exchange_name") if isinstance(properties, dict) else None  # type: ignore
        if exchange_name:
            row["exchange_reaction"] = exchange_name
        else:
            raise ValueError(f"Missing exchange_name for '{supp}' in file mapping '{str(custom_mapping_file)}'.")

        # Extract iupac_name
        row["iupac_name"] = properties.get("iupac_name", "")
        # Extract and convert mass_per_litre to mmol if requested
        mass_per_litre = properties.get("mass_per_litre") if isinstance(properties, dict) else 0.0
        row["mass_per_litre"] = float(mass_per_litre)  # type: ignore
        if "pubchem_name" in properties:
            molecular_weight = find_molecular_weight(properties["pubchem_name"])  # type: ignore
        elif "pubchem_id" in properties:
            molecular_weight = find_molecular_weight_by_id(properties["pubchem_id"])  # type: ignore
        else:
            molecular_weight = find_molecular_weight(row["name"])
        mmol = mg_to_mmol(float(mass_per_litre), molecular_weight) if molecular_weight else 0.0  # type: ignore
        row["mmol_concentration"] = mmol  # type: ignore

        # Extract flux_upper_bound if present
        if "flux_upper_bound" in properties and isinstance(properties, dict):
            row["flux_upper_bound"] = float(properties["flux_upper_bound"])  # type: ignore
        else:
            row["flux_upper_bound"] = 0.0  # type: ignore

        # update source. unique_supplements from metadata has precedence
        # over file/custom mapping. The default is UNSTATED
        if row["name"] in unique_supplements:
            row["source"] = MediumSource.SUPPLEMENT.value
        elif "source" in properties:
            row["source"] = str(properties["source"])  # type: ignore
        else:
            row["source"] = MediumSource.UNSTATED.value

        return row

    ############################################################################
    # Level 2: File mapping (overrides base)
    #
    # Load mappings from YAML file if provided, or use default (done once!)
    file_mapping: dict[str, dict[str, Any]] = {}
    if custom_mapping_file is None:
        # Default to custom_exchange_mapping.yaml in yamls directory
        default_yaml_path = Path(__file__).parent / "yamls" / "custom_exchange_mapping.yaml"
        if default_yaml_path.exists():
            custom_mapping_file = default_yaml_path

    if custom_mapping_file is not None and str(custom_mapping_file).strip():
        try:
            import yaml

            yaml_path = Path(custom_mapping_file)
            if yaml_path.exists():
                with open(yaml_path, encoding="utf-8") as f:
                    yaml_data = yaml.safe_load(f)
                    if yaml_data and "custom_mapping" in yaml_data:
                        file_mapping = yaml_data["custom_mapping"] or {}
        except Exception as e:
            logging.warning(f"Failed to load custom mapping file {custom_mapping_file}: {e}")
    #########
    for name, properties in file_mapping.items():
        if name is None:
            continue
        # salts are handled by selecting their solute ions
        if properties.get("is_salt", False):  # type: ignore
            rows = create_new_rows_for_salts(df_mapping, name, properties)
            # Update or append to df_mapping
            for row in rows:
                df_mapping = update_mapping_df(df_mapping, row)
        else:
            row = create_new_row(name, properties)
            # Update or append to df_mapping
            df_mapping = update_mapping_df(df_mapping, row)
    ############################################################################
    # Level 3: Custom mapping parameter (highest precedence, overrides all)
    if custom_mapping:
        for name, properties in custom_mapping.items():
            if name is None:
                continue

            # salts are handled by selecting their solute ions
            if properties.get("is_salt", False):  # type: ignore
                rows = create_new_rows_for_salts(df_mapping, name, properties)
                # Update or append to df_mapping
                for row in rows:
                    df_mapping = update_mapping_df(df_mapping, row)
            else:
                row = create_new_row(name, properties)
                # Update or append to df_mapping
                df_mapping = update_mapping_df(df_mapping, row)
    #####################################################################################
    # Level 4: unique_supplements are provided from metadate. If any is missing, add them
    #          with default properties
    map_keys = list(df_mapping["name"].values)
    fuzzy_matches = []
    unmapped = []
    for name in unique_supplements:
        df, _ = find_in_df_mapping(df_mapping, name)
        if df.empty:
            # Fuzzy match against SBML/manual keys
            matches = difflib.get_close_matches(name, map_keys, n=1, cutoff=fuzzy_threshold)
            if matches:
                best = matches[0]
                df_row = df_mapping.loc[df_mapping["name"] == best].iloc[0]
                row = {
                    "name": name,
                    "iupac_name": df_row["iupac_name"],
                    "exchange_reaction": df_row["exchange_reaction"],
                    "mass_per_litre": df_row["mass_per_litre"],
                    "mmol_concentration": df_row["mmol_concentration"],
                    "flux_upper_bound": df_row["flux_upper_bound"],
                    "source": MediumSource.SUPPLEMENT.value,
                }
                # Update or append to df_mapping
                df_mapping = update_mapping_df(df_mapping, row)
                fuzzy_matches.append((name, best, row))
            else:
                message = (
                    f"The exchange mapping element '{name}' found in metadata '{growth_rates_df.Name}', \n"
                    f" but missing in mappings file '{custom_mapping_file}', \n"
                    f" and organism's mapping provided as 'supplement_to_exchange_map' argument."
                )
                logging.warning(message)
                warnings.warn(message, UserWarning, stacklevel=2)
                unmapped.append(name)
    if verbose:
        print("\n" + "=" * 70)
        print("Medium source to Exchange Reaction Mapping")
        print("=" * 70)

        if fuzzy_matches:
            print(f"\nFuzzy matches ({len(fuzzy_matches)}):")
            for name, matched_key, r in fuzzy_matches:
                print(f"  {name:30s} -> {r['exchange_reaction']:50s} (matched via '{matched_key}')")

        if unmapped:
            print(f"\nUnmapped supplements ({len(unmapped)}):")
            for name in unmapped:
                print(f"  {name}")

        print("\n" + "=" * 70 + "\n")

    if fixed_exchanges:
        # Mark baseline exchanges
        df_mapping.loc[df_mapping["exchange_reaction"].isin(fixed_exchanges), "source"] = MediumSource.FIXED.value

    #######################################################################################
    # Sanity check: each row must be eather is_supplement or is_baseline, both cannot be True
    # df_subsets = df_mapping[(df_mapping["is_supplement"] & df_mapping["is_baseline"])]
    # if not df_subsets.empty:
    #     for _, row in df_subsets.iterrows():
    #         message = (
    #             f"Medium source '{row['name']}' cannot be both is_supplement:'True' and is_baseline:'True'. \n"
    #             f"Check mappings file '{custom_mapping_file}' and  metadata."
    #         )
    #         logging.warning(message)
    #         warnings.warn(message, UserWarning, stacklevel=2)
    #####################################################################################
    return df_mapping


def build_supplement_flux_dataframe(
    growth_rates_df: pd.DataFrame,
    mappings_df: pd.DataFrame,
    supplement_column: str = "supplements",
    growth_rate_column: str = "mu_max",
    success_column: str = "success",
    max_od600_column: str = "max_value",
    max_time_column: str = "max_time",
    total_volume_column: str = "total_volume_uL",
    od600_conversion_rate: float = 0.4,
    separator: str = ";",
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    # Create a dictionary of all sources that are either supplement, medium, or fixed
    all_exchanges = mappings_df.loc[
        (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
        | (mappings_df["source"] == MediumSource.MEDIUM.value)
        | (mappings_df["source"] == MediumSource.FIXED.value),
        "exchange_reaction",
    ]
    mediums_exchanges = mappings_df.loc[(mappings_df["source"] == MediumSource.MEDIUM.value), "exchange_reaction"]

    def calculate_flux(mapping_row):
        rxn = mapping_row["exchange_reaction"].values[0]
        mmol_value = mapping_row["mmol_concentration"].values[0]
        flux_upper_bound = mapping_row["flux_upper_bound"].values[0]
        # ensure column exists
        if rxn not in r_mmol_per_gCWD_per_time:
            r_mmol_per_gCWD_per_time[rxn] = 0.0
        # Add mmol value for this supplement
        max_od = row.get(max_od600_column, 0.0)
        max_time = row.get(max_time_column, 0.0)
        avg_od600 = max_od / 2
        # od600 to gCDW per litre
        gCWD_per_litre = od600_to_gCDW(avg_od600, od600_conversion_rate)  # noqa: N806
        # gCWD per liter to gCWD
        gCWD = (  # noqa: N806
            gCWD_per_litre * (row.get(total_volume_column, 1.0) / 1e6)  # convert uL to L  # noqa: N806
        )
        return rxn, flux_upper_bound, mmol_value / (gCWD * max_time) if gCWD > 0 and max_time > 0 else 0.0

    # Build rows with mmol values
    rows_mmol: list[dict[str, float]] = []
    for _, row in growth_rates_df.iterrows():
        # Skip rows where success column exists and is False
        if success_column in growth_rates_df.columns and not row[success_column]:
            continue
        # create a row's dictionary for all exchanges
        r_mmol_per_gCWD_per_time: dict[str, float] = dict.fromkeys(all_exchanges, 0.0)  # noqa: N806
        # parse supplements for this row
        supps_val = row.get(supplement_column, None)
        if pd.notna(supps_val):
            for part in [s.strip().lower() for s in str(supps_val).split(separator) if s.strip()]:
                mapping_row = mappings_df.loc[(mappings_df["name"] == part) | (mappings_df["iupac_name"] == part)]
                if mapping_row.empty:
                    continue

                rxn, flux_upper_bound, flux = calculate_flux(mapping_row)
                # flux in mmol / gCWD / time
                r_mmol_per_gCWD_per_time[rxn] = flux
                if r_mmol_per_gCWD_per_time[rxn] == 0.0 and flux_upper_bound > 0.0:
                    r_mmol_per_gCWD_per_time[rxn] = flux_upper_bound
        for ex in mediums_exchanges:
            mapping_row = mappings_df.loc[mappings_df["exchange_reaction"] == ex]
            if mapping_row.empty:
                continue

            rxn, flux_upper_bound, flux = calculate_flux(mapping_row)
            # flux in mmol / gCWD / time
            r_mmol_per_gCWD_per_time[rxn] = flux
            if r_mmol_per_gCWD_per_time[rxn] == 0.0 and flux_upper_bound > 0.0:
                r_mmol_per_gCWD_per_time[rxn] = flux_upper_bound
        rows_mmol.append(r_mmol_per_gCWD_per_time)

    result_df_mmol = pd.DataFrame(rows_mmol)
    # ensure fixed columns are present (even if empty)
    fixed_exchanges = mappings_df.loc[mappings_df["source"] == MediumSource.FIXED.value]
    for ex, flux_upper_bound in fixed_exchanges[["exchange_reaction", "flux_upper_bound"]].itertuples(index=False):
        if ex not in result_df_mmol.columns:
            if flux_upper_bound > 0.0:
                result_df_mmol[ex] = flux_upper_bound  # noqa: N806
            else:
                result_df_mmol[ex] = 0.0
        else:
            # Make sure they are not supplements, since some of the supplement values
            # are none-zero
            if result_df_mmol[ex].sum() == 0.0:
                result_df_mmol[ex] = flux_upper_bound

    # Reorder columns: exchanges sorted, then growth rate last
    exch_cols = sorted([c for c in result_df_mmol.columns if c != growth_rate_column])
    # Apply suffix to exchange columns if requested
    if exchange_suffix:
        rename_map = {col: f"{col}{exchange_suffix}" for col in exch_cols}
        result_df_mmol = result_df_mmol.rename(columns=rename_map)
        exch_cols = [f"{col}{exchange_suffix}" for col in exch_cols]
    # Add growth rate column
    result_df_mmol[growth_rate_column] = growth_rates_df.loc[
        growth_rates_df[success_column] if success_column in growth_rates_df.columns else growth_rates_df.index,
        growth_rate_column,
    ].values

    # Reorder columns: exchanges sorted, then growth rate last
    result_df_mmol = result_df_mmol[exch_cols + [growth_rate_column]]

    return result_df_mmol


def build_AMN_inputs_dataframe(  # noqa: N802
    growth_rates_df: pd.DataFrame,
    mappings_df: pd.DataFrame,
    supplement_column: str = "supplements",
    growth_rate_column: str = "mu_max",
    success_column: str = "success",
    separator: str = ";",
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    # Create a dictionary of all sources that are either supplement, medium, or fixed
    all_exchanges = mappings_df.loc[
        (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
        | (mappings_df["source"] == MediumSource.MEDIUM.value)
        | (mappings_df["source"] == MediumSource.FIXED.value),
        "exchange_reaction",
    ]
    # Build rows
    rows: list[dict[str, int | float | None]] = []
    for _, row in growth_rates_df.iterrows():
        # Skip rows where success column exists and is False
        if success_column in growth_rates_df.columns and not row[success_column]:
            continue
        # create a row's dictionary for all exchanges
        new_row: dict[str, int | float | None] = dict.fromkeys(all_exchanges, 0)
        # parse supplements for this row
        supps_val = row.get(supplement_column)
        if pd.notna(supps_val):
            for part in [s.strip().lower() for s in str(supps_val).split(separator) if s.strip()]:
                names = mappings_df.loc[(mappings_df["name"] == part) | (mappings_df["iupac_name"] == part)]
                if names.empty:
                    continue
                rxn = names["exchange_reaction"].values[0]
                # ensure column exists
                # if rxn not in row:
                #     row[rxn] = 0
                new_row[rxn] = 1

        # add growth rate as final column (preserve column name passed)
        new_row[growth_rate_column] = row.get(growth_rate_column)

        rows.append(new_row)

    result_df = pd.DataFrame(rows)
    # ensure medium columns are present (even if empty)
    medium_exchanges = mappings_df.loc[mappings_df["source"] == MediumSource.MEDIUM.value]
    for ex in medium_exchanges["exchange_reaction"]:
        result_df[ex] = 1
    # ensure fixed columns are present (even if empty)
    fixed_exchanges = mappings_df.loc[mappings_df["source"] == MediumSource.FIXED.value]
    for ex in fixed_exchanges["exchange_reaction"]:
        result_df[ex] = 1

    # Reorder columns: exchanges sorted, then growth rate last
    exch_cols = sorted([c for c in result_df.columns if c != growth_rate_column])

    # Apply suffix to exchange columns if requested
    if exchange_suffix:
        rename_map = {col: f"{col}{exchange_suffix}" for col in exch_cols}
        result_df = result_df.rename(columns=rename_map)
        exch_cols = [f"{col}{exchange_suffix}" for col in exch_cols]

    result_df = result_df[exch_cols + [growth_rate_column]]

    return result_df


def build_AMN_levels_dataframe(  # noqa: N802
    exchange_matrix: pd.DataFrame,
    mappings_df: pd.DataFrame,
    flux_df: pd.DataFrame | None = None,
    growth_rate_column: str = "mu_max",
    default_level: int = 1,
    default_max_value: int = 1000,
    sbml_bounds: dict[str, tuple[int, int]] | None = None,
    custom_bounds: dict[str, tuple[int, int]] | None = None,
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    # Create a dictionary of all flux upper bound where exchange are either supplement, medium, or fixed
    flux_upper_bounds = (
        mappings_df.loc[
            (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
            | (mappings_df["source"] == MediumSource.MEDIUM.value)
            | (mappings_df["source"] == MediumSource.FIXED.value),
            ["exchange_reaction", "flux_upper_bound"],
        ]
        .set_index("exchange_reaction")["flux_upper_bound"]
        .to_dict()
    )

    if flux_df is not None:
        # Override with flux_df values if provided and greater than zero
        for col in flux_df.columns:
            trimmed_col = col.replace(exchange_suffix, "") if exchange_suffix else col
            if trimmed_col in flux_upper_bounds and flux_df[col].max() > 0.0:
                flux_upper_bounds[trimmed_col] = flux_df[col].max()

    # Get all exchange columns (exclude growth rate column)
    exchange_cols = [col for col in exchange_matrix.columns if col != growth_rate_column]
    # Create the template dataframe
    template_data: dict[str, list[str | float]] = {"name": ["level", "max_value", "ratio_drawing"]}

    # Add columns for each exchange reaction with values
    for col in exchange_cols:
        trimmed_col = col.replace(exchange_suffix, "") if exchange_suffix else col
        # Precedence: custom_bounds > mappings_df.flux_upper_bound > sbml_bounds > defaults
        if custom_bounds and col in custom_bounds:
            # Highest priority: custom bounds
            level_val, max_val = custom_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        elif trimmed_col in flux_upper_bounds and flux_upper_bounds.get(trimmed_col, 0.0) > 0.0:
            # Second priority: Use upper bound from mappings_df.flux_upper_bound in fluxes if available
            max_val = flux_upper_bounds.get(trimmed_col, 0)
            template_data[col] = [default_level, max_val, 0]
        elif sbml_bounds and col in sbml_bounds:
            # Third priority: SBML bounds
            level_val, max_val = sbml_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        else:
            # Lowest priority: default values
            template_data[col] = [default_level, default_max_value, 0]

    return pd.DataFrame(template_data)


def parse_sbml_exchanges(sbml_path: str | Path) -> dict[str, str]:
    """
    Parse SBML file to extract exchange reactions and metabolite names.

    Parameters
    ----------
    sbml_path : Union[str, Path]
        Path to SBML XML file

    Returns
    -------
    Dict[str, str]
        Mapping from common metabolite names to exchange reaction IDs

    Notes
    -----
    This function attempts to use COBRApy if available, otherwise falls back
    to XML parsing. Exchange reactions are identified by the "EX_" prefix.

    Examples
    --------
    >>> mapping = parse_sbml_exchanges("iML1515.xml")
    >>> print(mapping["glucose"])  # EX_glc__D_e
    """
    try:
        import cobra

        # Suppress COBRApy INFO logging
        logging.getLogger("cobra").setLevel(logging.WARNING)

        model = cobra.io.read_sbml_model(str(sbml_path))

        mapping = {}
        for reaction in model.reactions:
            if reaction.id.startswith("EX_"):
                # Get metabolite names from the reaction
                for metabolite in reaction.metabolites:
                    # Use the metabolite's common name or ID
                    name = metabolite.name if metabolite.name else metabolite.id
                    # Clean up the name
                    name = name.replace("_e", "").replace("_", " ").strip()
                    mapping[name.lower()] = reaction.id

        return mapping

    except ImportError:
        print("Warning: COBRApy not available, attempting XML parsing")
        return _parse_sbml_exchanges_fallback(sbml_path)


def _parse_sbml_exchanges_fallback(sbml_path: str | Path) -> dict[str, str]:
    """
    Fallback XML parser for SBML files when COBRApy is not available.

    Parameters
    ----------
    sbml_path : Union[str, Path]
        Path to SBML XML file

    Returns
    -------
    Dict[str, str]
        Mapping from common metabolite names to exchange reaction IDs
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(str(sbml_path))
    root = tree.getroot()

    # Find namespace
    namespace = {"sbml": root.tag.split("}")[0].strip("{")} if "}" in root.tag else {}

    mapping = {}

    # Parse reactions
    for reaction in root.findall(".//sbml:reaction", namespace):
        rxn_id = reaction.get("id", "")
        if rxn_id.startswith("EX_"):
            rxn_name = reaction.get("name", "")
            if rxn_name:
                # Clean up the name
                name = rxn_name.replace("_e", "").replace("_", " ").strip()
                mapping[name.lower()] = rxn_id

    return mapping


def parse_sbml_exchange_bounds(
    sbml_path: str | Path,
    default_level: int = 1,
) -> dict[str, tuple[int, int]]:
    """
    Parse SBML file to extract exchange reaction bounds (upper bounds).

    This function extracts the upper bound values for exchange reactions from an SBML file
    and returns them in a format suitable for use with create_exchange_bounds_template.

    Parameters
    ----------
    sbml_path : Union[str, Path]
        Path to SBML XML file
    default_level : int
        Default value to use for the "level" parameter (default: 1)
        The level value is not stored in SBML, so this default is used for all reactions

    Returns
    -------
    Dict[str, Tuple[int, int]]
        Mapping from exchange reaction IDs to (level, max_value) tuples.
        Example: {"EX_glc__D_e": (1, 1000), "EX_o2_e_i": (1, 500)}

    Notes
    -----
    This function attempts to use COBRApy if available, otherwise falls back
    to XML parsing. Exchange reactions are identified by the "EX_" prefix.

    The upper_bound from SBML is used as max_value. If the upper bound is
    infinite or very large (>10000), it's capped at 1000 for practical use.

    Examples
    --------
    >>> # Parse bounds from SBML
    >>> bounds = parse_sbml_exchange_bounds("iML1515.xml")
    >>> print(bounds["EX_glc__D_e"])  # (1, 1000)
    >>>
    >>> # Use with create_exchange_bounds_template
    >>> matrix = create_supplement_exchange_matrix(growth_df, mapping)
    >>> template = create_exchange_bounds_template(
    ...     matrix,
    ...     sbml_bounds=bounds
    ... )
    """
    try:
        import cobra

        # Suppress COBRApy INFO logging
        logging.getLogger("cobra").setLevel(logging.WARNING)

        model = cobra.io.read_sbml_model(str(sbml_path))

        bounds_map = {}
        for reaction in model.reactions:
            if reaction.id.startswith("EX_"):
                # Get upper bound (max_value)
                upper_bound = reaction.upper_bound

                # Cap very large or infinite bounds at 1000
                if upper_bound > 10000 or upper_bound == float("inf"):
                    max_value = 1000
                else:
                    max_value = int(upper_bound)

                # Store as (level, max_value) tuple
                bounds_map[reaction.id] = (default_level, max_value)

        return bounds_map

    except ImportError:
        print("Warning: COBRApy not available, attempting XML parsing")
        return _parse_sbml_exchange_bounds_fallback(sbml_path, default_level)


def _parse_sbml_exchange_bounds_fallback(
    sbml_path: str | Path,
    default_level: int = 1,
) -> dict[str, tuple[int, int]]:
    """
    Fallback XML parser for SBML exchange bounds when COBRApy is not available.

    Parameters
    ----------
    sbml_path : Union[str, Path]
        Path to SBML XML file
    default_level : int
        Default value to use for the "level" parameter

    Returns
    -------
    Dict[str, Tuple[int, int]]
        Mapping from exchange reaction IDs to (level, max_value) tuples
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(str(sbml_path))
    root = tree.getroot()

    # Find namespace
    namespace = {"sbml": root.tag.split("}")[0].strip("{")} if "}" in root.tag else {}

    bounds_map = {}

    # Parse reactions
    for reaction in root.findall(".//sbml:reaction", namespace):
        rxn_id = reaction.get("id", "")
        if rxn_id.startswith("EX_"):
            # Try to find upper bound in kinetic law or default to 1000
            upper_bound = 1000

            # Look for upperFluxBound attribute or parameter
            kinetic_law = reaction.find(".//sbml:kineticLaw", namespace)
            if kinetic_law is not None:
                for param in kinetic_law.findall(".//sbml:parameter", namespace):
                    param_id = param.get("id", "")
                    if "upper" in param_id.lower() or "ub" in param_id.lower():
                        try:
                            value = float(param.get("value", "1000"))
                            if value > 10000 or value == float("inf"):
                                upper_bound = 1000
                            else:
                                upper_bound = int(value)
                        except (ValueError, TypeError):
                            upper_bound = 1000

            bounds_map[rxn_id] = (default_level, upper_bound)

    return bounds_map


def load_default_iml1515_mapping() -> dict[str, str]:
    """
    Load a default mapping for iML1515 E. coli model.

    Returns
    -------
    Dict[str, str]
        Default mapping from common supplement names to exchange reactions

    Notes
    -----
    This is a curated mapping for common supplements used in E. coli experiments.
    """
    return {
        # Sugars / Carbon sources
        "glucose": "EX_glc__D_e",
        "fructose": "EX_fru_e_i",
        "galactose": "EX_gal_e_i",
        "ribose": "EX_rib__D_e_i",
        "maltose": "EX_malt_e_i",
        "melibiose": "EX_melib_e_i",
        "trehalose": "EX_tre_e_i",
        "glycerol": "EX_glyc_e_i",
        # Organic acids
        "acetate": "EX_ac_e_i",
        "pyruvate": "EX_pyr_e_i",
        "succinate": "EX_succ_e_i",
        "lactate": "EX_lac__D_e_i",
        "d-lactate": "EX_lac__D_e_i",
        # Amino acids
        "alanine": "EX_ala__L_e_i",
        "l-alanine": "EX_ala__L_e_i",
        "proline": "EX_pro__L_e_i",
        "l-proline": "EX_pro__L_e_i",
        "threonine": "EX_thr__L_e_i",
        "l-threonine": "EX_thr__L_e_i",
        "glycine": "EX_gly_e_i",
        # Nucleobases
        "adenine": "EX_ade_e",
        "uracil": "EX_ura_e",
        "guanine": "EX_gua_e",
        "cytosine": "EX_csn_e",
        "thymine": "EX_thymd_e",
    }


def load_minimal_media_exchanges() -> list[str]:
    """
    Load standard minimal media exchange reactions for E. coli.

    Returns
    -------
    List[str]
        List of exchange reaction IDs that are typically always available
        in minimal media (ions, minerals, gases, etc.)

    Notes
    -----
    These are components that are usually present in all media conditions
    and represent the minimal requirements for growth.
    """
    return [
        "EX_pi_e_i",  # Phosphate
        "EX_co2_e_i",  # Carbon dioxide
        "EX_fe3_e_i",  # Iron (III)
        "EX_h_e_i",  # Proton
        "EX_mn2_e_i",  # Manganese
        "EX_fe2_e_i",  # Iron (II)
        "EX_zn2_e_i",  # Zinc
        "EX_mg2_e_i",  # Magnesium
        "EX_ca2_e_i",  # Calcium
        "EX_ni2_e_i",  # Nickel
        "EX_cu2_e_i",  # Copper
        "EX_sel_e_i",  # Selenium
        "EX_cobalt2_e_i",  # Cobalt
        "EX_h2o_e_i",  # Water
        "EX_mobd_e_i",  # Molybdate
        "EX_so4_e_i",  # Sulfate
        "EX_nh4_e_i",  # Ammonium
        "EX_k_e_i",  # Potassium
        "EX_na1_e_i",  # Sodium
        "EX_cl_e_i",  # Chloride
        "EX_o2_e_i",  # Oxygen
        "EX_tungs_e_i",  # Tungsten
        "EX_slnt_e_i",  # Selenite
    ]
