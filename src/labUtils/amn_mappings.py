##############################################################
#  labUtils: Metabolic network mapping utilities             #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
#                                                            #
##############################################################

from __future__ import annotations

import difflib
import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd

from labUtils.utils import (
    find_molecular_weight,
    find_molecular_weight_by_id,
    g_to_mol,
    get_compound_by_name,
    od600_to_gCDW,  # noqa: F401  (re-export expected by some callers)
)

# pylint: disable=import-outside-toplevel


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class MediumSource(Enum):
    UNSTATED = "unstated"
    SUPPLEMENT = "supplement"
    MEDIUM = "medium"
    FIXED = "fixed"


class RecordOrigin(Enum):
    SMBL = "smbl"
    UPDATED = "updated"


class SearchResult(Enum):
    NAME = "name"
    IUPAC = "iupac"
    FUZZY_NAME = "fuzzy_name"
    FUZZY_IUPAC = "fuzzy_iupac"
    NOT_FOUND = "not_found"


# ---------------------------------------------------------------------------
# Output containers
# ---------------------------------------------------------------------------


@dataclass
class MappingTable:
    """Typed wrapper around the mapping dataframe produced by `build_mapping_table`.

    Attributes
    ----------
    df
        The full mapping dataframe (same shape as `build_mappings` returns).
    organism
        Optional organism tag (informational only; not used for routing).
    yaml_source
        Path of the YAML file the mapping was overlaid from, if any.
    """

    df: pd.DataFrame
    organism: str | None = None
    yaml_source: Path | None = None

    def __iter__(self):  # backwards-compat helper: `df, = MappingTable(...)`
        return iter([self.df])


# ---------------------------------------------------------------------------
# iML1515 / BiGG defaults
# ---------------------------------------------------------------------------


IML1515_DEFAULT_SUPPLEMENT_MAP: dict[str, str] = {
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


IML1515_MINIMAL_MEDIA: list[str] = [
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


# ---------------------------------------------------------------------------
# Public: build_mappings (and dataclass-returning sibling)
# ---------------------------------------------------------------------------


def build_mappings(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: dict[str, str],
    supplement_column: str = "supplements",
    custom_mapping_file: str | Path | None = None,
    custom_mapping: dict[str, str | dict[str, Any]] | None = None,
    separator: str = ";",
    fuzzy_threshold: float = 0.6,
    halt_on_error: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """Build the supplement-to-exchange mapping dataframe.

    See module docstring / characterization tests for the precise contract; this
    function preserves the existing behavior (column order, precedence, sanity
    checks) and is consumed by `build_supplement_flux_dataframe`,
    `build_AMN_inputs_dataframe`, and `build_AMN_levels_dataframe`.
    """
    # Level 1 — seed from supplement_to_exchange_map.
    logging.info("Function build_mappings is started")
    # ===================================================================================
    ## Level 1: Base mapping from supplement_to_exchange_map
    # Create a dataframe from organism's supplement to exchange mapping
    supp_exchange = (
        (
            supp.strip().lower(),
            "",
            set(),
            exchange,
            0.0,  # mass_per_litre
            0.0,  # mmol_concentration
            0.0,  # flux_upper_bound
            MediumSource.UNSTATED.value,
            RecordOrigin.SMBL.value,
        )
        for supp, exchange in supplement_to_exchange_map.items()
        if supp is not None
    )
    df_mapping = pd.DataFrame(
        list(supp_exchange),
        columns=[
            "name",
            "iupac_name",
            "other_names",
            "exchange_reaction",
            "mass_per_litre",
            "mmol_concentration",
            "flux_upper_bound",
            "source",
            "record_origin",
        ],
    )

    # Level 2 — overlay YAML file (default or explicit).
    file_mapping = _load_yaml_mappings(custom_mapping_file)
    df_mapping = _apply_input_mapping(file_mapping, df_mapping, fuzzy_threshold, halt_on_error, verbose)

    # Level 3 — overlay caller-provided custom_mapping (highest precedence).
    # ===================================================================================
    def alert(msg: str, halt_on_error: bool = True):
        if halt_on_error:
            logging.error(msg)
            raise ValueError(msg)
        if verbose:
            # warnings.warn((msg), stacklevel=2)
            logging.warning(msg)
        else:
            logging.info(msg)

    # ===================================================================================
    def search_by_name_in_df_mapping(df_mapping, row) -> tuple[pd.DataFrame, SearchResult]:
        # Search the row by its name, iupac_name, or fuzzy matching
        # Level 1: Try name
        name = row["name"]
        df = df_mapping.loc[df_mapping["name"] == name, :]
        if not df.empty:
            return df, SearchResult.NAME
        # Level 2: Try iupac_name
        iupac_name = row["iupac_name"].strip()
        if iupac_name != "":
            df = df_mapping.loc[df_mapping["iupac_name"] == iupac_name, :]
            if not df.empty:
                return df, SearchResult.IUPAC
        # Level 3: Try fuzzy matching
        map_keys = list(df_mapping["name"].values)
        # names vs names
        matches = difflib.get_close_matches(name, map_keys, n=1, cutoff=fuzzy_threshold)
        if matches:
            best = matches[0]
            return df_mapping.loc[df_mapping["name"] == best, :], SearchResult.FUZZY_NAME
        # names vs iupac_names
        if iupac_name != "":
            matches = difflib.get_close_matches(iupac_name, map_keys, n=1, cutoff=fuzzy_threshold)
            if matches:
                best = matches[0]
                return df_mapping.loc[df_mapping["name"] == best, :], SearchResult.FUZZY_NAME

        # iupac_names vs iupac_names
        map_keys = list(df_mapping["iupac_name"].values)
        if iupac_name != "":
            matches = difflib.get_close_matches(iupac_name, map_keys, n=1, cutoff=fuzzy_threshold)
            if matches:
                best = matches[0]
                return df_mapping.loc[df_mapping["iupac_name"] == best, :], SearchResult.FUZZY_IUPAC
        # iupac_names vs names
        matches = difflib.get_close_matches(name, map_keys, n=1, cutoff=fuzzy_threshold)
        if matches:
            best = matches[0]
            return df_mapping.loc[df_mapping["name"] == best, :], SearchResult.FUZZY_IUPAC

        return pd.DataFrame(), SearchResult.NOT_FOUND

    def update_mapping_df_by_row(df_mapping, row, halt_on_not_found: bool = True):
        def set_entire_row(df_mapping, new_row):
            # update the entire row
            df_mapping.loc[existing_rows.index[0]] = new_row
            return df_mapping

        def add_new_row(df_mapping, new_row):
            df_mapping = pd.concat([df_mapping, pd.DataFrame([new_row])], ignore_index=True)
            return df_mapping

        def update_row_by_existing(row, existing_row):
            row["other_names"] |= {row["name"], row["iupac_name"]}
            row["other_names"] = {item for item in row["other_names"] if item != ""}
            row["name"] = existing_row["name"]
            row["iupac_name"] = existing_row["iupac_name"]
            row["exchange_reaction"] = existing_row["exchange_reaction"]
            return row

        # The df_mapping contains all exchange_reactions of the organism. So,
        # there should be at least one in them that matches the new row.
        # After making sure there is one, based on the name search, we decide how to update it.
        #
        #
        # Check the exchange_reactions existing in the df_mapping
        if (
            row["exchange_reaction"].strip() != ""
            and df_mapping.loc[df_mapping["exchange_reaction"] == row["exchange_reaction"]].empty
        ):
            alert(
                f"There is no exchange reaction:'{row['exchange_reaction']}' in the organism. "
                f"This is occured for mapping element:'{row['name']}' (iupac name:'{row['iupac_name']}'). ",
                halt_on_not_found,
            )
            # If the exchange reaction is not found, we cannot proceed.
            return df_mapping
        # When exchange reaction exists, update if the name matches,
        # otherwise, add it as a new record.
        if row["exchange_reaction"].strip() != "":
            existing_rows, search_result = search_by_name_in_df_mapping(df_mapping, row)
            if existing_rows.empty or search_result in {
                SearchResult.NOT_FOUND,
                SearchResult.FUZZY_NAME,
                SearchResult.FUZZY_IUPAC,
            }:
                df_mapping = add_new_row(df_mapping, row)
                return df_mapping
            else:
                if existing_rows.shape[0] > 1:
                    alert(
                        f"Multiple existing rows found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}'). Using the first one.",
                        halt_on_error=False,
                    )
                existing_row = existing_rows.iloc[0]
                row = update_row_by_existing(row, existing_row)
                df_mapping = set_entire_row(df_mapping, row)
                return df_mapping

        ###############################################################
        #       Search by name, iupac_name, or fuzzy matching
        existing_rows, search_result = search_by_name_in_df_mapping(df_mapping, row)

        if existing_rows.empty and row["exchange_reaction"].strip() == "":
            alert(
                f"Could not found a match for name:'{row['name']}' (iupac_name:'{row['iupac_name']}').",
                halt_on_error=halt_on_not_found,
            )
            return df_mapping
        elif existing_rows.empty:  # Empty means not found, so, add new row
            return add_new_row(df_mapping, row)
        # there is at least one existing row
        if existing_rows.shape[0] > 1:
            alert(
                f"Multiple existing rows found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}'). Using the first one.",
                halt_on_error=False,
            )
        existing_row = existing_rows.iloc[0]
        match search_result:
            # If there is a match by name or iupac_name, replace it by the new row.
            # However, there exchange_reaction must be the same, otherwise, it is a mistake in the mapping.
            case SearchResult.NAME | SearchResult.IUPAC:
                if verbose:
                    alert(
                        f"exact match found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}') to "
                        f"existing mapping name:'{existing_row['name']}', iupac_name:'{existing_row['iupac_name']}', exchange reaction:'{existing_row['exchange_reaction']}'.",
                        halt_on_error=False,
                    )
                row = update_row_by_existing(row, existing_row)
                df_mapping = set_entire_row(df_mapping, row)
            case SearchResult.FUZZY_NAME | SearchResult.FUZZY_IUPAC:
                alert(
                    f"fuzzy match found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}') to "
                    f"existing mapping name:'{existing_row['name']}', iupac_name:'{existing_row['iupac_name']}', exchange reaction:'{existing_row['exchange_reaction']}'.",
                    halt_on_error=False,
                )
                row = update_row_by_existing(row, existing_row)
                df_mapping = set_entire_row(df_mapping, row)
            case SearchResult.NOT_FOUND:
                alert(
                    f"Could not found a match for name:'{row['name']}' (iupac_name:'{row['iupac_name']}').",
                    halt_on_error=halt_on_not_found,
                )
        logging.info("Function update_mapping_df_by_row finished successfully")
        return df_mapping

    def create_new_row(name, properties):

        row = {"name": name.strip().lower()}
        # Extract exchange_name from nested dict
        exchange_name: str = properties.get("exchange_name", "")
        row["exchange_reaction"] = exchange_name
        # Extract iupac_name
        row["iupac_name"] = properties.get("iupac_name", "")
        # other names
        row["other_names"] = set(properties.get("other_names", []))
        # Extract and convert mass_per_litre to mmol if requested
        mass_per_litre = properties.get("mass_per_litre") if isinstance(properties, dict) else 0.0
        row["mass_per_litre"] = float(mass_per_litre)  # type: ignore
        # mmol_per_liter / mol_per_liter override the PubChem-derived calculation entirely.
        # Precedence: mmol_per_liter > mol_per_liter > mass_per_litre → MW lookup.
        if "mmol_per_liter" in properties:
            row["mmol_concentration"] = float(properties["mmol_per_liter"])  # type: ignore
        elif "mol_per_liter" in properties:
            row["mmol_concentration"] = float(properties["mol_per_liter"]) * 1000.0  # type: ignore
        else:
            if "pubchem_name" in properties:
                molecular_weight = find_molecular_weight(properties["pubchem_name"])  # type: ignore
            elif "pubchem_id" in properties:
                molecular_weight = find_molecular_weight_by_id(properties["pubchem_id"])  # type: ignore
            elif "iupac_name" in properties and properties["iupac_name"] != "":
                molecular_weight = find_molecular_weight(properties["iupac_name"])  # type: ignore
            else:
                molecular_weight = find_molecular_weight(row["name"])
            mol = g_to_mol(float(mass_per_litre), molecular_weight) if molecular_weight else 0.0  # type: ignore
            mmol = mol * 1000.0  # convert mol to mmol
            row["mmol_concentration"] = mmol  # type: ignore

        # Extract flux_upper_bound if present
        if "flux_upper_bound" in properties and isinstance(properties, dict):
            row["flux_upper_bound"] = float(properties["flux_upper_bound"])  # type: ignore
        else:
            row["flux_upper_bound"] = 0.0  # type: ignore

        if "source" in properties:
            row["source"] = str(properties["source"])  # type: ignore
        else:
            row["source"] = MediumSource.UNSTATED.value

        # set record_origin
        row["record_origin"] = RecordOrigin.UPDATED.value

        return row

    def create_new_rows_for_compounds(df_mapping, compound_name, properties):
        rows = []
        compound_name = compound_name.strip().lower()
        if len(properties.get("compounds", [])) == 0:
            raise ValueError(
                f"The exchange name:'{compound_name}' is marked as compound but no compounds were provided."
            )
        compounds = []
        for subcompound in properties.get("compounds", []):
            df = df_mapping.loc[(df_mapping["name"] == subcompound) | (df_mapping["iupac_name"] == subcompound), :]
            if df.empty:
                alert(
                    f"Could not find compound name:'{subcompound}' for compound '{compound_name}' in mappings.",
                    halt_on_error=halt_on_error,
                )
                continue
            if df.shape[0] > 1:
                alert(
                    f"Multiple entries found for compound name:'{subcompound}' for compound '{compound_name}'. Using the first one.",
                    halt_on_error=False,
                )
            existing_row = df.iloc[0]
            if existing_row["iupac_name"]:
                compound = get_compound_by_name(existing_row["iupac_name"])
            elif existing_row["name"]:
                compound = get_compound_by_name(existing_row["name"])
            else:
                raise ValueError(
                    f"The exchange name:'{existing_row['name']}' is marked as compound but no valid name or ID was provided for searching pubchem db."
                )

            compounds.append({"compound": compound, "mapping_row": existing_row})

        # total_molar_mass = sum(info["molar_mass"] for compound_entry in compounds for info in compound_entry["info"])
        #  molecular_weight in g/mol
        total_molar_mass = sum(compound_entry["compound"].molecular_weight for compound_entry in compounds)
        for compound_entry in compounds:
            compound = compound_entry["compound"]
            existing_row = compound_entry["mapping_row"]
            #  molecular_weight in g/mol
            ratio = compound.molecular_weight / total_molar_mass if total_molar_mass > 0 else 0.0
            # Adjust mass_per_litre based on ratio.
            # IMPORTANT: The mass_per_litre is additive for compounds.
            mass_per_litre = properties.get("mass_per_litre", 0.0)
            corrected_mass_per_litre = float(mass_per_litre) * ratio + existing_row.get("mass_per_litre", 0.0)
            # New row properties
            row_properties = {
                "name": existing_row["name"],
                "iupac_name": existing_row["iupac_name"],
                "other_names": existing_row["other_names"],
                "exchange_name": existing_row["exchange_reaction"],
                "mass_per_litre": corrected_mass_per_litre,
                "mmol_concentration": 0.0,
                "flux_upper_bound": 0.0,
                "source": MediumSource.MEDIUM.value,
                "record_origin": RecordOrigin.UPDATED.value,
            }
            row = create_new_row(existing_row["name"], row_properties)
            rows.append(row)

        return rows

    def update_df_mappings(input_mapping: dict[str, str | dict[str, Any]], df_mapping):
        for name, properties in input_mapping.items():
            if name is None:
                continue
            # compounds are handled differently
            if properties.get("is_compound", False):  # type: ignore
                rows = create_new_rows_for_compounds(df_mapping, name, properties)
                # Update or append to df_mapping
                for row in rows:
                    df_mapping = update_mapping_df_by_row(df_mapping, row, halt_on_not_found=halt_on_error)
            else:
                row = create_new_row(name, properties)
                # Update or append to df_mapping
                df_mapping = update_mapping_df_by_row(df_mapping, row, halt_on_not_found=halt_on_error)
        return df_mapping

    # ===================================================================================
    # Level 2: File mapping (overrides base)
    #
    # Load mappings from YAML file if provided, or use default (done once!)
    file_mapping: dict[str, str | dict[str, Any]] = {}
    if custom_mapping_file is None:
        # Default to custom_exchange_mapping.yaml in yamls directory
        default_yaml_path = Path(__file__).parent / "yamls" / "custom_exchange_mapping.yaml"
        if default_yaml_path.exists():
            custom_mapping_file = default_yaml_path

    if custom_mapping_file is not None and str(custom_mapping_file).strip():
        try:
            import yaml

            yaml_path = Path(__file__).parent / "yamls" / custom_mapping_file
            if yaml_path.exists():
                with open(yaml_path, encoding="utf-8") as f:
                    yaml_data = yaml.safe_load(f)
                    if yaml_data and "custom_mapping" in yaml_data:
                        file_mapping = yaml_data["custom_mapping"] or {}
        except Exception as e:
            logging.warning(f"Failed to load custom mapping file {custom_mapping_file}: {e}")
    #
    # Use the file mapping to update the dataframe
    df_mapping = update_df_mappings(file_mapping, df_mapping)
    # ===================================================================================
    # Level 3: Custom mapping parameter (highest precedence, overrides all)
    if custom_mapping:
        df_mapping = _apply_input_mapping(custom_mapping, df_mapping, fuzzy_threshold, halt_on_error, verbose)

    # Level 4 — mark observed supplements from growth_rates_df.
    unique_supplements: set[str] = set()
    if supplement_column in growth_rates_df.columns:
        for val in growth_rates_df[supplement_column].dropna().astype(str):
            parts = [s.strip() for s in val.split(separator) if s.strip()]
            for p in parts:
                unique_supplements.add(p.lower())

    unique_supplements_mapping: dict[str, dict[str, Any]] = {}
    for name in unique_supplements:
        row = {
            "name": name,
            "iupac_name": "",
            "other_names": set(),
            "exchange_name": "",
            "mass_per_litre": 0.0,
            "mmol_concentration": 0.0,
            "flux_upper_bound": 0.0,
            "source": MediumSource.SUPPLEMENT.value,
            "record_origin": RecordOrigin.UPDATED.value,
        }
        existing_rows, _ = _search_mapping_by_name(df_mapping, row, fuzzy_threshold)
        if not existing_rows.empty:
            df_mapping.loc[existing_rows.index[0], "source"] = MediumSource.SUPPLEMENT.value
            df_mapping.loc[existing_rows.index[0], "record_origin"] = RecordOrigin.UPDATED.value
            df_mapping.at[existing_rows.index[0], "other_names"].add(name)  # type: ignore[union-attr]
        else:
            unique_supplements_mapping[name] = _build_mapping_row(name, row)
    df_mapping = _apply_input_mapping(unique_supplements_mapping, df_mapping, fuzzy_threshold, halt_on_error, verbose)

    # Sanity check 1 — at most one FIXED record per exchange reaction.
    fixed_exchanges = df_mapping.loc[df_mapping["source"] == MediumSource.FIXED.value, "exchange_reaction"].unique()
    for exchange in fixed_exchanges:
        df_fixed = df_mapping.loc[
            (df_mapping["exchange_reaction"] == exchange) & (df_mapping["source"] == MediumSource.FIXED.value), :
        ]
        if df_fixed.shape[0] > 1:
            _log_alert(
                f"Sanity check warning: Multiple FIXED records found for exchange reaction:'{exchange}': \n{df_fixed.to_string()}. "
                f"Please ensure only one FIXED record per exchange reaction.",
                halt_on_error,
                verbose,
            )

    # Sanity check 2 — within UPDATED records, one MEDIUM/FIXED per exchange.
    df_sanity = df_mapping[
        (
            (df_mapping["record_origin"] == RecordOrigin.UPDATED.value)
            & ((df_mapping["source"] == MediumSource.MEDIUM.value) | (df_mapping["source"] == MediumSource.FIXED.value))
        )
    ]
    for exchange_reaction, group in df_sanity.groupby("exchange_reaction"):
        if group.shape[0] > 1 and len(group["source"].unique()) != group.shape[0]:
            _log_alert(
                f"Sanity check warning: Multiple sources found for exchange reaction:'{exchange_reaction}': \n{group.to_string()}. "
                f"Please ensure only one of them is defined as SUPPLEMENT, MEDIUM, or FIXED.",
                halt_on_error=False,
                verbose=verbose,
            )

    # Normalize other_names: sorted, semicolon-joined string.
    df_mapping["other_names"] = df_mapping["other_names"].apply(lambda x: sorted(x))
    df_mapping["other_names"] = df_mapping["other_names"].apply(lambda x: ";".join(x))
    df_mapping = df_mapping.sort_values(by=["exchange_reaction"]).reset_index(drop=True)

    if verbose:
        df_report = df_mapping[
            (
                (df_mapping["record_origin"] == RecordOrigin.UPDATED.value)
                & (df_mapping["source"] != MediumSource.UNSTATED.value)
            )
        ]
        df_report = df_report.sort_values(by=["source", "exchange_reaction"]).reset_index(drop=True)
        logging.info(f"Final updated mapping dataframe:\n{df_report.to_string()}")
    # ===================================================================================
    logging.info("Function build_mappings finished successfully")
    return df_mapping


def build_mapping_table(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: dict[str, str],
    *,
    organism: str | None = None,
    supplement_column: str = "supplements",
    custom_mapping_file: str | Path | None = None,
    custom_mapping: dict[str, str | dict[str, Any]] | None = None,
    separator: str = ";",
    fuzzy_threshold: float = 0.6,
    halt_on_error: bool = False,
    verbose: bool = False,
) -> MappingTable:
    """Dataclass-returning wrapper around `build_mappings`.

    Same behavior, but returns a `MappingTable` so downstream callers can carry
    along the organism tag and yaml source metadata.
    """
    df = build_mappings(
        growth_rates_df,
        supplement_to_exchange_map,
        supplement_column=supplement_column,
        custom_mapping_file=custom_mapping_file,
        custom_mapping=custom_mapping,
        separator=separator,
        fuzzy_threshold=fuzzy_threshold,
        halt_on_error=halt_on_error,
        verbose=verbose,
    )
    yaml_source: Path | None = None
    if custom_mapping_file is not None and str(custom_mapping_file).strip():
        yaml_source = Path(__file__).parent / "yamls" / custom_mapping_file
        if not yaml_source.exists():
            yaml_source = None
    return MappingTable(df=df, organism=organism, yaml_source=yaml_source)


# ---------------------------------------------------------------------------
# Public: build_supplement_flux_dataframe
# ---------------------------------------------------------------------------


def build_supplement_flux_dataframe(
    growth_rates_df: pd.DataFrame,
    mappings_df: pd.DataFrame,
    growth_rate_column: str = "mu_max",
    success_column: str = "success",
    max_od600_column: str = "max_value",
    max_time_column: str = "max_time",
    supplement_column: str = "supplements",
    total_volume_column: str = "total_volume_uL",
    max_growth_rate_column: str = "mv_mu_max_value",
    od600_conversion_rate: float = 0.4,  # noqa: ARG001 (kept for backwards compatibility)
    separator: str = ";",
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    logging.info("Function build_supplement_flux_dataframe is started")
    # ==================================================================================
    # the precedence order is SUPPLEMENT > MEDIUM > FIXED, and alll UNSTATED are ignored
    #

    # ==================================================================================
    # Create dictionary and lists for supplements, mediums, fixed exchanges

    # Create a dictionary of all sources that are either supplement, medium, or fixed
    valid_mappings = mappings_df.loc[
        (
            (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
            | (mappings_df["source"] == MediumSource.MEDIUM.value)
            | (mappings_df["source"] == MediumSource.FIXED.value)
        )
        & (mappings_df["source"] != MediumSource.UNSTATED.value),
        :,
    ]
    all_exchanges = valid_mappings["exchange_reaction"]
    supplement_exchanges = supplement_mappings["exchange_reaction"]
    mediums_exchanges = medium_mappings["exchange_reaction"]

    columns: dict[str, float] = dict.fromkeys(all_exchanges, 0.0)

    rows: list[dict[str, float]] = []
    for _, existing_row in growth_rates_df.iterrows():
        if success_column in growth_rates_df.columns and not existing_row[success_column]:
            continue
        new_row: dict[str, float] = columns.copy()
        for column in new_row:
            if column in fixed_exchanges["exchange_reaction"].values:
                flux_upper_bound = fixed_exchanges[fixed_exchanges["exchange_reaction"] == column][
                    "flux_upper_bound"
                ].values[0]
                new_row[column] = flux_upper_bound
            if column in mediums_exchanges.values:
                mapping_row = medium_mappings.loc[medium_mappings["exchange_reaction"] == column]
                if mapping_row.empty:
                    continue
                exchange_reaction, _, flux = _calculate_flux_row(
                    mapping_row, existing_row, max_time_column, max_growth_rate_column,
                )
                new_row[exchange_reaction] = flux
            if column in supplement_exchanges.values:
                mapping_row = supplement_mappings.loc[supplement_mappings["exchange_reaction"] == column]
                if mapping_row.empty:
                    continue
                supps_vals = existing_row.get(supplement_column, [])
                parts = [s.strip().lower() for s in str(supps_vals).split(separator) if s.strip()]
                names = [
                    supplement_mappings.loc[
                        (supplement_mappings["name"] == part) | (supplement_mappings["iupac_name"] == part),
                        "exchange_reaction",
                    ].values
                    for part in parts
                ]
                names = {item for sublist in names for item in sublist}
                if column in names:
                    mapping_row = pd.concat([mapping_row.loc[mapping_row["name"] == part, :] for part in parts])
                    if mapping_row.shape[0] > 1:
                        logging.warning(
                            f"Multiple mapping rows found for supplement exchange:'{column}' for parts:{parts}. Using the first one.\n"
                            f"{mapping_row.to_string()}"
                        )
                    exchange_reaction, _, flux = _calculate_flux_row(
                        mapping_row, existing_row, max_time_column, max_growth_rate_column,
                    )
                    new_row[exchange_reaction] = flux
                else:
                    new_row[column] = 0
        rows.append(new_row)

    result_df = pd.DataFrame(rows)
    exch_cols = sorted([c for c in result_df.columns if c != growth_rate_column])
    if exchange_suffix:
        result_df = result_df.rename(columns={c: f"{c}{exchange_suffix}" for c in exch_cols})
        exch_cols = [f"{c}{exchange_suffix}" for c in exch_cols]
    result_df[growth_rate_column] = growth_rates_df.loc[
        growth_rates_df[success_column] if success_column in growth_rates_df.columns else growth_rates_df.index,
        growth_rate_column,
    ].values
    result_df = result_df[exch_cols + [growth_rate_column]]

    logging.info("Function build_supplement_flux_dataframe finished successfully")
    return result_df


# ---------------------------------------------------------------------------
# Public: build_AMN_inputs_dataframe
# ---------------------------------------------------------------------------


def build_AMN_inputs_dataframe(  # noqa: N802
    growth_rates_df: pd.DataFrame,
    mappings_df: pd.DataFrame,
    supplement_column: str = "supplements",
    growth_rate_column: str = "mu_max",
    success_column: str = "success",
    separator: str = ";",
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    logging.info("Function build_AMN_inputs_dataframe is started")
    # ==================================================================================
    # the precedence order is SUPPLEMENT > MEDIUM > FIXED, and alll UNSTATED are ignored
    #

    # ==================================================================================
    # Create dictionary and lists for supplements, mediums, fixed exchanges

    # Create a dictionary of all sources that are either supplement, medium, or fixed
    valid_mappinsgs = mappings_df.loc[
        (
            (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
            | (mappings_df["source"] == MediumSource.MEDIUM.value)
            | (mappings_df["source"] == MediumSource.FIXED.value)
        )
        & (mappings_df["source"] != MediumSource.UNSTATED.value),
        :,
    ]
    all_exchanges = valid_mappinsgs["exchange_reaction"]
    supplement_mappings = valid_mappinsgs.loc[(valid_mappinsgs["source"] == MediumSource.SUPPLEMENT.value), :]
    # createlists of exchanges
    supplement_exchanges = valid_mappinsgs.loc[
        (valid_mappinsgs["source"] == MediumSource.SUPPLEMENT.value), "exchange_reaction"
    ]
    mediums_exchanges = valid_mappinsgs.loc[
        (valid_mappinsgs["source"] == MediumSource.MEDIUM.value), "exchange_reaction"
    ]
    fixed_exchanges = valid_mappinsgs.loc[
        valid_mappinsgs["source"] == MediumSource.FIXED.value, ["exchange_reaction", "flux_upper_bound"]
    ]
    # ==================================================================================
    # Build dataframe columns for all exchanges and initialize to zero
    columns: dict[str, float] = dict.fromkeys(all_exchanges, 0.0)
    rows: list[dict[str, float]] = []
    for _, existing_row in growth_rates_df.iterrows():
        if success_column in growth_rates_df.columns and not existing_row[success_column]:
            continue
        new_row: dict[str, float] = columns.copy()
        for column in new_row:
            new_row[column] = 0
            if column in fixed_exchanges["exchange_reaction"].values:
                new_row[column] = 1
            if column in mediums_exchanges.values:
                new_row[column] = 1
            if column in supplement_exchanges.values:
                supps_vals = existing_row.get(supplement_column, [])
                parts = (s.strip().lower() for s in str(supps_vals).split(separator) if s.strip())
                names = [
                    supplement_mappings.loc[
                        (supplement_mappings["name"] == part) | (supplement_mappings["iupac_name"] == part),
                        "exchange_reaction",
                    ].values
                    for part in parts
                ]
                names = {item for sublist in names for item in sublist}
                if column in names:
                    new_row[column] = 1
                else:
                    new_row[column] = 0
        rows.append(new_row)

    result_df = pd.DataFrame(rows)
    exch_cols = sorted([c for c in result_df.columns if c != growth_rate_column])
    if exchange_suffix:
        result_df = result_df.rename(columns={c: f"{c}{exchange_suffix}" for c in exch_cols})
        exch_cols = [f"{c}{exchange_suffix}" for c in exch_cols]
    result_df[growth_rate_column] = growth_rates_df.loc[
        growth_rates_df[success_column] if success_column in growth_rates_df.columns else growth_rates_df.index,
        growth_rate_column,
    ].values
    result_df = result_df[exch_cols + [growth_rate_column]]

    logging.info("Function build_AMN_inputs_dataframe finished successfully")
    return result_df


# ---------------------------------------------------------------------------
# Public: build_AMN_levels_dataframe + multi-level generator
# ---------------------------------------------------------------------------


def build_AMN_levels_dataframe(  # noqa: N802
    exchange_matrix: pd.DataFrame,
    mappings_df: pd.DataFrame,
    flux_df: pd.DataFrame,
    growth_rate_column: str = "mu_max",
    default_level: int = 1,
    default_variable_level: int = 100,
    default_max_value: int = 1000,
    sbml_bounds: dict[str, tuple[int, int]] | None = None,
    custom_bounds: dict[str, tuple[int, int]] | None = None,
    exchange_suffix: str | None = None,
) -> pd.DataFrame:
    logging.info("Function build_AMN_levels_dataframe is started")
    def trim(item: str) -> str:
        if not exchange_suffix:
            return item
        ln = len(exchange_suffix)
        return item[:-ln] if item[-ln:] == exchange_suffix else item

    # Pick one row per exchange via precedence SUPPLEMENT > MEDIUM > FIXED.
    single_record_mappings_df = mappings_df.loc[
        (
            (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
            | (mappings_df["source"] == MediumSource.MEDIUM.value)
            | (mappings_df["source"] == MediumSource.FIXED.value)
        )
        & (mappings_df["source"] != MediumSource.UNSTATED.value),
        :,
    ]
    single_record_mappings_df = single_record_mappings_df.sort_values(
        by=["exchange_reaction", "source"],
        key=lambda x: x.map(
            {
                MediumSource.SUPPLEMENT.value: 0,
                MediumSource.MEDIUM.value: 1,
                MediumSource.FIXED.value: 2,
            }
        ),
    )
    single_record_mappings_df = single_record_mappings_df.drop_duplicates(subset=["exchange_reaction"], keep="first")
    flux_upper_bounds = (
        single_record_mappings_df[["exchange_reaction", "flux_upper_bound"]]
        .set_index("exchange_reaction")["flux_upper_bound"]
        .to_dict()
    )
    flux_sources = (
        single_record_mappings_df[["exchange_reaction", "source"]].set_index("exchange_reaction")["source"].to_dict()
    )
    for col in flux_df.columns:
        trimmed_col = trim(col)
        if trimmed_col in flux_upper_bounds and flux_df[col].max() > 0.0:
            flux_upper_bounds[trimmed_col] = flux_df[col].max()

    exchange_cols = [col for col in exchange_matrix.columns if col != growth_rate_column]
    template_data: dict[str, list[str | float]] = {"name": ["level", "max_value", "ratio_drawing"]}
    for col in exchange_cols:
        trimmed_col = trim(col)
        if custom_bounds and col in custom_bounds:
            level_val, max_val = custom_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        elif trimmed_col in flux_upper_bounds:
            max_val = flux_upper_bounds.get(trimmed_col, 0)
            if flux_sources.get(trimmed_col, "") in {MediumSource.FIXED.value, MediumSource.MEDIUM.value}:
                template_data[col] = [1, max_val, 0]
            else:
                template_data[col] = [default_variable_level, max_val, 0]
        elif sbml_bounds and col in sbml_bounds:
            level_val, max_val = sbml_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        else:
            template_data[col] = [default_level, default_max_value, 0]

    logging.info("Function build_AMN_levels_dataframe finished successfully")
    return pd.DataFrame(template_data)


def build_AMN_levels_for_different_levels(  # noqa: N802
    df_AMN_levels: pd.DataFrame,  # noqa: N803
    fixed_levels: list[int] | None = None,
    fixed_max_values: list[int] | None = None,
    default_level: int = 1,
    default_max_value: int = 20,
):
    """Generate multiple AMN levels dataframes for different fixed levels and max values."""
    logging.info("Function build_AMN_levels_for_different_levels is started")
    new_dfs = []
    if fixed_levels is None:
        logging.info("Function build_AMN_levels_for_different_levels finished successfully")
        return new_dfs

    if fixed_max_values is None:
        fixed_max_values = [default_max_value]

    if fixed_levels is None:  # pragma: no cover (kept for behavioral parity)
        fixed_levels = [default_level]

    for level in fixed_levels:
        for max_value in fixed_max_values:
            df_level = df_AMN_levels.copy().reset_index(drop=True)
            for col in df_level.columns:
                if col == "name":
                    continue
                if df_level.at[0, col] == 1:
                    df_level.at[1, col] = max_value
                else:
                    df_level.at[0, col] = level
            new_dfs.append(df_level)

    logging.info("Function build_AMN_levels_for_different_levels finished successfully")
    return new_dfs


# ---------------------------------------------------------------------------
# Public: SBML parsing
# ---------------------------------------------------------------------------


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
    logging.info("Function parse_sbml_exchanges is started")
    try:
        import cobra

        logging.getLogger("cobra").setLevel(logging.WARNING)
        model = cobra.io.read_sbml_model(str(sbml_path))
        mapping: dict[str, str] = {}
        for reaction in model.reactions:
            if reaction.id.startswith("EX_"):
                for metabolite in reaction.metabolites:
                    name = metabolite.name if metabolite.name else metabolite.id
                    name = name.replace("_e", "").replace("_", " ").strip()
                    mapping[name.lower()] = reaction.id

                # logging.info("Function labUtils.amn_mappings.parse_sbml_exchanges finished successfully")
        return mapping
    except ImportError:
        logging.warning("COBRApy not available, attempting XML parsing for SBML exchanges")
        logging.info("Function parse_sbml_exchanges finished with error: COBRApy not available, using fallback XML parsing")
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
    logging.info("Function parse_sbml_exchange_bounds is started")
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
                max_value = 1000 if upper_bound > 10000 or upper_bound == float("inf") else int(upper_bound)
                # Store as (level, max_value) tuple
                bounds_map[reaction.id] = (default_level, max_value)

        logging.debug("Function parse_sbml_exchange_bounds finished successfully")
        return bounds_map

    except ImportError:
        logging.warning("COBRApy not available, attempting XML parsing for SBML exchanges")
        logging.debug("Function parse_sbml_exchange_bounds finished successfully")
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
    namespace = {"sbml": root.tag.split("}")[0].strip("{")} if "}" in root.tag else {}

    mapping: dict[str, str] = {}
    for reaction in root.findall(".//sbml:reaction", namespace):
        rxn_id = reaction.get("id", "")
        if rxn_id.startswith("EX_"):
            rxn_name = reaction.get("name", "")
            if rxn_name:
                name = rxn_name.replace("_e", "").replace("_", " ").strip()
                mapping[name.lower()] = rxn_id
    return mapping


def _parse_sbml_exchange_bounds_fallback(
    sbml_path: str | Path, default_level: int = 1
) -> dict[str, tuple[int, int]]:
    """XML-only fallback for `parse_sbml_exchange_bounds`."""
    import xml.etree.ElementTree as ET

    tree = ET.parse(str(sbml_path))
    root = tree.getroot()
    namespace = {"sbml": root.tag.split("}")[0].strip("{")} if "}" in root.tag else {}

    bounds_map: dict[str, tuple[int, int]] = {}
    for reaction in root.findall(".//sbml:reaction", namespace):
        rxn_id = reaction.get("id", "")
        if rxn_id.startswith("EX_"):
            upper_bound = 1000
            kinetic_law = reaction.find(".//sbml:kineticLaw", namespace)
            if kinetic_law is not None:
                for param in kinetic_law.findall(".//sbml:parameter", namespace):
                    param_id = param.get("id", "")
                    if "upper" in param_id.lower() or "ub" in param_id.lower():
                        try:
                            value = float(param.get("value", "1000"))
                            upper_bound = 1000 if value > 10000 or value == float("inf") else int(value)
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
    logging.debug("Function load_default_iml1515_mapping finished successfully")
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
    logging.debug("Function load_minimal_media_exchanges finished successfully")
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
