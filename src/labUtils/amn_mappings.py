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
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, cast

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
    organism_sbml_path: str | Path | None = None,
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
    logging.info("Function build_mappings is started")
    # Level 1 — seed from supplement_to_exchange_map.
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

    sbml_lookup: dict[str, str] = {}
    sbml_exchange_ids: set[str] = set()
    resolution_stats = MappingResolutionStats()

    # Track direct resolutions from seed map.
    resolution_stats.direct = len(supplement_to_exchange_map)

    if organism_sbml_path is not None and str(organism_sbml_path).strip():
        try:
            sbml_mapping = parse_sbml_exchanges(organism_sbml_path)
            sbml_lookup = {
                _normalize_mapping_key(name): exchange
                for name, exchange in sbml_mapping.items()
                if str(name).strip()
            }
            sbml_exchange_ids = set(sbml_mapping.values())
            logging.info(
                "Loaded %d SBML exchange names and %d SBML exchange ids for mapping fallback.",
                len(sbml_lookup),
                len(sbml_exchange_ids),
            )
        except Exception as exc:
            logging.warning("Failed to load SBML exchange mapping for fallback resolution: %s", exc)

    df_mapping = _harmonize_seed_mappings_with_sbml(
        df_mapping,
        sbml_lookup,
        sbml_exchange_ids,
        fuzzy_threshold,
        resolution_stats=resolution_stats,
    )

    # Level 2 — overlay YAML file (default or explicit).
    file_mapping = _load_yaml_mappings(custom_mapping_file)
    df_mapping = _apply_input_mapping(
        file_mapping,
        df_mapping,
        fuzzy_threshold,
        halt_on_error,
        verbose,
        sbml_lookup=sbml_lookup,
        sbml_exchange_ids=sbml_exchange_ids,
        resolution_stats=resolution_stats,
    )

    # Level 3 — overlay caller-provided custom_mapping (highest precedence).
    if custom_mapping:
        df_mapping = _apply_input_mapping(
            custom_mapping,
            df_mapping,
            fuzzy_threshold,
            halt_on_error,
            verbose,
            sbml_lookup=sbml_lookup,
            sbml_exchange_ids=sbml_exchange_ids,
            resolution_stats=resolution_stats,
        )

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
    df_mapping = _apply_input_mapping(
        unique_supplements_mapping,
        df_mapping,
        fuzzy_threshold,
        halt_on_error,
        verbose,
        sbml_lookup=sbml_lookup,
        sbml_exchange_ids=sbml_exchange_ids,
    )

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

    logging.info("Function build_mappings finished successfully")
    resolution_stats.log_summary(verbose=verbose)
    return df_mapping


def build_mapping_table(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: dict[str, str],
    *,
    organism: str | None = None,
    supplement_column: str = "supplements",
    custom_mapping_file: str | Path | None = None,
    custom_mapping: dict[str, str | dict[str, Any]] | None = None,
    organism_sbml_path: str | Path | None = None,
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
        organism_sbml_path=organism_sbml_path,
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
    """Convert mappings + per-experiment growth kinetics into an exchange-flux dataframe.

    Output: one row per successful growth experiment, columns are sorted exchange IDs
    plus a trailing growth-rate column. See characterization tests for exact contract.
    """
    logging.info("Function build_supplement_flux_dataframe is started")
    valid_mappings, supplement_mappings, medium_mappings, fixed_exchanges = _split_by_source(mappings_df)
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
    """Binary input matrix indicating which exchanges are active per experiment."""
    logging.info("Function build_AMN_inputs_dataframe is started")
    valid_mappings, supplement_mappings, medium_mappings, fixed_exchanges = _split_by_source(mappings_df)
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
    new_dfs: list[pd.DataFrame] = []
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
    """Parse SBML to extract `metabolite_common_name -> exchange_reaction_id` mapping.

    Uses COBRApy if available; otherwise falls back to lightweight XML parsing.
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
        return mapping
    except ImportError:
        logging.warning("COBRApy not available, attempting XML parsing for SBML exchanges")
        return _parse_sbml_exchanges_fallback(sbml_path)


def parse_sbml_exchange_bounds(
    sbml_path: str | Path,
    default_level: int = 1,
) -> dict[str, tuple[int, int]]:
    """Parse SBML to extract `exchange_id -> (level, max_value)` for AMN templates."""
    logging.info("Function parse_sbml_exchange_bounds is started")
    try:
        import cobra

        logging.getLogger("cobra").setLevel(logging.WARNING)
        model = cobra.io.read_sbml_model(str(sbml_path))
        bounds_map: dict[str, tuple[int, int]] = {}
        for reaction in model.reactions:
            if reaction.id.startswith("EX_"):
                upper_bound = reaction.upper_bound
                max_value = 1000 if upper_bound > 10000 or upper_bound == float("inf") else int(upper_bound)
                bounds_map[reaction.id] = (default_level, max_value)
        logging.debug("Function parse_sbml_exchange_bounds finished successfully")
        return bounds_map
    except ImportError:
        logging.warning("COBRApy not available, attempting XML parsing for SBML exchanges")
        return _parse_sbml_exchange_bounds_fallback(sbml_path, default_level)


# ---------------------------------------------------------------------------
# Public: iML1515 defaults (curated)
# ---------------------------------------------------------------------------


def load_default_iml1515_mapping() -> dict[str, str]:
    """Curated supplement -> exchange mapping for the iML1515 E. coli model."""
    logging.debug("Function load_default_iml1515_mapping finished successfully")
    return dict(IML1515_DEFAULT_SUPPLEMENT_MAP)


def load_minimal_media_exchanges() -> list[str]:
    """Standard minimal-media exchange reactions for E. coli."""
    logging.debug("Function load_minimal_media_exchanges finished successfully")
    return list(IML1515_MINIMAL_MEDIA)


# ---------------------------------------------------------------------------
# Private helpers — logging / alerts
# ---------------------------------------------------------------------------


def _log_alert(msg: str, halt_on_error: bool = True, verbose: bool = False) -> None:
    """Single source of truth for the verbose/halt/warning logging dance."""
    if halt_on_error:
        logging.error(msg)
        raise ValueError(msg)
    if verbose:
        logging.warning(msg)
    else:
        logging.info(msg)


# ---------------------------------------------------------------------------
# Private helpers — mapping construction
# ---------------------------------------------------------------------------


def _search_mapping_by_name(
    df_mapping: pd.DataFrame, row: dict[str, Any], fuzzy_threshold: float
) -> tuple[pd.DataFrame, SearchResult]:
    """Search the mapping dataframe by exact name, exact iupac, then fuzzy."""
    name = row["name"]
    df = df_mapping.loc[df_mapping["name"] == name, :]
    if not df.empty:
        return cast(pd.DataFrame, df), SearchResult.NAME

    iupac_name = row["iupac_name"].strip()
    if iupac_name != "":
        df = df_mapping.loc[df_mapping["iupac_name"] == iupac_name, :]
        if not df.empty:
            return cast(pd.DataFrame, df), SearchResult.IUPAC

    map_keys = list(df_mapping["name"].values)
    matches = difflib.get_close_matches(name, map_keys, n=1, cutoff=fuzzy_threshold)
    if matches:
        return cast(pd.DataFrame, df_mapping.loc[df_mapping["name"] == matches[0], :]), SearchResult.FUZZY_NAME
    if iupac_name != "":
        matches = difflib.get_close_matches(iupac_name, map_keys, n=1, cutoff=fuzzy_threshold)
        if matches:
            return cast(pd.DataFrame, df_mapping.loc[df_mapping["name"] == matches[0], :]), SearchResult.FUZZY_NAME

    map_keys = list(df_mapping["iupac_name"].values)
    if iupac_name != "":
        matches = difflib.get_close_matches(iupac_name, map_keys, n=1, cutoff=fuzzy_threshold)
        if matches:
            return cast(pd.DataFrame, df_mapping.loc[df_mapping["iupac_name"] == matches[0], :]), SearchResult.FUZZY_IUPAC
    matches = difflib.get_close_matches(name, map_keys, n=1, cutoff=fuzzy_threshold)
    if matches:
        return cast(pd.DataFrame, df_mapping.loc[df_mapping["name"] == matches[0], :]), SearchResult.FUZZY_IUPAC

    return pd.DataFrame(), SearchResult.NOT_FOUND


def _merge_row_into_existing(row: dict[str, Any], existing_row: pd.Series) -> dict[str, Any]:
    """Update `row` in place with synonyms from existing, then mirror existing's identity."""
    row["other_names"] |= {row["name"], row["iupac_name"]}
    row["other_names"] = {item for item in row["other_names"] if item != ""}
    row["name"] = existing_row["name"]
    row["iupac_name"] = existing_row["iupac_name"]
    row["exchange_reaction"] = existing_row["exchange_reaction"]
    return row


def _append_mapping_row(df_mapping: pd.DataFrame, row: dict[str, Any]) -> pd.DataFrame:
    """Append one row to `df_mapping` without using concat on empty/all-NA frames."""
    row_df = pd.DataFrame([row], columns=df_mapping.columns)
    if df_mapping.empty:
        return row_df
    df_mapping.loc[len(df_mapping)] = row_df.iloc[0]
    return df_mapping


@dataclass
class MappingResolutionStats:
    """Tracking for mapping resolution stages."""
    direct: int = 0
    sbml_exact: int = 0
    sbml_fuzzy: int = 0
    unresolved: int = 0

    def log_summary(self, verbose: bool = False) -> None:
        """Emit a summary report of resolution stages."""
        total = self.direct + self.sbml_exact + self.sbml_fuzzy + self.unresolved
        if total == 0:
            return
        pct_direct = (self.direct / total) * 100 if total > 0 else 0
        pct_sbml_exact = (self.sbml_exact / total) * 100 if total > 0 else 0
        pct_sbml_fuzzy = (self.sbml_fuzzy / total) * 100 if total > 0 else 0
        pct_unresolved = (self.unresolved / total) * 100 if total > 0 else 0

        msg = (
            f"Mapping resolution summary: "
            f"direct={self.direct}({pct_direct:.0f}%), "
            f"sbml_exact={self.sbml_exact}({pct_sbml_exact:.0f}%), "
            f"sbml_fuzzy={self.sbml_fuzzy}({pct_sbml_fuzzy:.0f}%), "
            f"unresolved={self.unresolved}({pct_unresolved:.0f}%) [total={total}]"
        )
        if verbose:
            logging.warning(msg)
        else:
            logging.info(msg)


def _normalize_mapping_key(value: str) -> str:
    """Normalize free-text compound keys for exact/fuzzy matching."""
    return " ".join(str(value).strip().lower().replace("_", " ").replace("-", " ").split())


def _resolve_exchange_from_sbml(
    row: dict[str, Any],
    sbml_lookup: dict[str, str],
    fuzzy_threshold: float,
) -> tuple[str | None, str | None]:
    """Resolve a mapping row to an SBML exchange by exact/fuzzy name lookup."""
    candidates: list[str] = []
    for field in ("name", "iupac_name"):
        val = str(row.get(field, "") or "").strip()
        if val:
            candidates.append(val)

    other_names = row.get("other_names")
    if isinstance(other_names, (set, list, tuple)):
        for val in other_names:
            txt = str(val).strip()
            if txt:
                candidates.append(txt)
    elif isinstance(other_names, str) and other_names.strip():
        candidates.extend([s.strip() for s in other_names.split(";") if s.strip()])

    for candidate in candidates:
        normalized = _normalize_mapping_key(candidate)
        if normalized in sbml_lookup:
            return sbml_lookup[normalized], "sbml_exact"

    if not candidates or not sbml_lookup:
        return None, None

    sbml_keys = list(sbml_lookup.keys())
    for candidate in candidates:
        normalized = _normalize_mapping_key(candidate)
        matches = difflib.get_close_matches(normalized, sbml_keys, n=1, cutoff=fuzzy_threshold)
        if matches:
            return sbml_lookup[matches[0]], "sbml_fuzzy"

    return None, None


def _harmonize_seed_mappings_with_sbml(
    df_mapping: pd.DataFrame,
    sbml_lookup: dict[str, str],
    sbml_exchange_ids: set[str],
    fuzzy_threshold: float,
    resolution_stats: MappingResolutionStats | None = None,
) -> pd.DataFrame:
    """Resolve seed exchanges that are not present in the organism SBML."""
    if df_mapping.empty or not sbml_lookup or not sbml_exchange_ids:
        return df_mapping

    for idx, existing_row in df_mapping.iterrows():
        exchange_reaction = str(existing_row.get("exchange_reaction", "") or "").strip()
        if not exchange_reaction or exchange_reaction in sbml_exchange_ids:
            continue

        row = {
            "name": str(existing_row.get("name", "") or ""),
            "iupac_name": str(existing_row.get("iupac_name", "") or ""),
            "other_names": existing_row.get("other_names", set()),
        }
        resolved_exchange, resolution_stage = _resolve_exchange_from_sbml(row, sbml_lookup, fuzzy_threshold)
        if resolved_exchange is None:
            continue

        df_mapping.at[idx, "exchange_reaction"] = resolved_exchange
        logging.info(
            "Resolved seed mapping exchange for '%s' from '%s' to '%s' via %s.",
            row["name"],
            exchange_reaction,
            resolved_exchange,
            resolution_stage,
        )
        if resolution_stats is not None:
            resolution_stats.direct = max(0, resolution_stats.direct - 1)
            if resolution_stage == "sbml_exact":
                resolution_stats.sbml_exact += 1
            elif resolution_stage == "sbml_fuzzy":
                resolution_stats.sbml_fuzzy += 1

    return df_mapping


def _upsert_mapping_row(
    df_mapping: pd.DataFrame,
    row: dict[str, Any],
    halt_on_not_found: bool,
    fuzzy_threshold: float,
    verbose: bool,
    sbml_lookup: dict[str, str] | None = None,
    sbml_exchange_ids: set[str] | None = None,
    resolution_stats: MappingResolutionStats | None = None,
) -> pd.DataFrame:
    """Insert or update a single mapping row.

    Reproduces the original update_mapping_df_by_row logic verbatim:
    - Validates that the row's exchange exists in df_mapping (unless empty).
    - If found by exact name/iupac, merges with existing.
    - If found by fuzzy or not at all, appends as new.
    """
    sbml_lookup = sbml_lookup or {}
    sbml_exchange_ids = sbml_exchange_ids or set()
    exchange_reaction = row["exchange_reaction"].strip()

    if exchange_reaction != "":
        present_in_mapping = not df_mapping.loc[df_mapping["exchange_reaction"] == exchange_reaction].empty
        present_in_sbml = exchange_reaction in sbml_exchange_ids
        if not present_in_mapping and not present_in_sbml:
            resolved_exchange, resolution_stage = _resolve_exchange_from_sbml(row, sbml_lookup, fuzzy_threshold)
            if resolved_exchange is not None:
                logging.info(
                    "Resolved mapping exchange for '%s' from '%s' to '%s' via %s.",
                    row["name"],
                    exchange_reaction,
                    resolved_exchange,
                    resolution_stage,
                )
                if resolution_stats is not None:
                    if resolution_stage == "sbml_exact":
                        resolution_stats.sbml_exact += 1
                    elif resolution_stage == "sbml_fuzzy":
                        resolution_stats.sbml_fuzzy += 1
                row["exchange_reaction"] = resolved_exchange
            else:
                if resolution_stats is not None:
                    resolution_stats.unresolved += 1
                _log_alert(
                    f"Could not validate mapping exchange '{exchange_reaction}' for mapping element "
                    f"'{row['name']}' (iupac_name='{row['iupac_name']}') against current mapping seed "
                    "or organism SBML.",
                    halt_on_not_found,
                    verbose,
                )
                return df_mapping

    if row["exchange_reaction"].strip() != "":
        existing_rows, search_result = _search_mapping_by_name(df_mapping, row, fuzzy_threshold)
        if existing_rows.empty or search_result in {
            SearchResult.NOT_FOUND,
            SearchResult.FUZZY_NAME,
            SearchResult.FUZZY_IUPAC,
        }:
            df_mapping = _append_mapping_row(df_mapping, row)
            return df_mapping
        if existing_rows.shape[0] > 1:
            _log_alert(
                f"Multiple existing rows found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}'). Using the first one.",
                halt_on_error=False,
                verbose=verbose,
            )
        existing_row = existing_rows.iloc[0]
        row = _merge_row_into_existing(row, existing_row)
        df_mapping.loc[existing_rows.index[0]] = row
        return df_mapping

    # Branch where row's exchange_reaction is empty: search by name/iupac, then update or warn.
    existing_rows, search_result = _search_mapping_by_name(df_mapping, row, fuzzy_threshold)

    if existing_rows.empty and row["exchange_reaction"].strip() == "":
        resolved_exchange, resolution_stage = _resolve_exchange_from_sbml(row, sbml_lookup, fuzzy_threshold)
        if resolved_exchange is not None:
            row["exchange_reaction"] = resolved_exchange
            logging.info(
                "Resolved missing exchange for '%s' to '%s' via %s.",
                row["name"],
                resolved_exchange,
                resolution_stage,
            )
            if resolution_stats is not None:
                if resolution_stage == "sbml_exact":
                    resolution_stats.sbml_exact += 1
                elif resolution_stage == "sbml_fuzzy":
                    resolution_stats.sbml_fuzzy += 1
            df_mapping = _append_mapping_row(df_mapping, row)
            return df_mapping

        if resolution_stats is not None:
            resolution_stats.unresolved += 1
        _log_alert(
            f"Could not resolve mapping for name '{row['name']}' (iupac_name='{row['iupac_name']}') "
            "after checking existing mappings and organism SBML.",
            halt_on_error=halt_on_not_found,
            verbose=verbose,
        )
        return df_mapping
    if existing_rows.empty:
        df_mapping = _append_mapping_row(df_mapping, row)
        return df_mapping

    if existing_rows.shape[0] > 1:
        _log_alert(
            f"Multiple existing rows found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}'). Using the first one.",
            halt_on_error=False,
            verbose=verbose,
        )
    existing_row = existing_rows.iloc[0]
    match search_result:
        case SearchResult.NAME | SearchResult.IUPAC:
            if verbose:
                _log_alert(
                    f"exact match found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}') to "
                    f"existing mapping name:'{existing_row['name']}', iupac_name:'{existing_row['iupac_name']}', exchange reaction:'{existing_row['exchange_reaction']}'.",
                    halt_on_error=False,
                    verbose=verbose,
                )
            row = _merge_row_into_existing(row, existing_row)
            df_mapping.loc[existing_rows.index[0]] = row
        case SearchResult.FUZZY_NAME | SearchResult.FUZZY_IUPAC:
            _log_alert(
                f"fuzzy match found for name:'{row['name']}' (iupac_name:'{row['iupac_name']}') to "
                f"existing mapping name:'{existing_row['name']}', iupac_name:'{existing_row['iupac_name']}', exchange reaction:'{existing_row['exchange_reaction']}'.",
                halt_on_error=False,
                verbose=verbose,
            )
            row = _merge_row_into_existing(row, existing_row)
            df_mapping.loc[existing_rows.index[0]] = row
        case SearchResult.NOT_FOUND:
            _log_alert(
                f"Could not resolve mapping for name '{row['name']}' (iupac_name='{row['iupac_name']}') "
                "after checking existing mappings and organism SBML.",
                halt_on_error=halt_on_not_found,
                verbose=verbose,
            )
    return df_mapping


def _build_mapping_row(name: str, properties: Any) -> dict[str, Any]:
    """Build a single mapping row from a `name`+`properties` pair (PubChem call hidden here)."""
    row: dict[str, Any] = {"name": name.strip().lower()}
    properties_dict = properties if isinstance(properties, dict) else {}
    if isinstance(properties, str):
        exchange_name = properties
    else:
        exchange_name = str(properties_dict.get("exchange_name", ""))
    row["exchange_reaction"] = exchange_name
    row["iupac_name"] = str(properties_dict.get("iupac_name", ""))
    other_names = properties_dict.get("other_names", [])
    if isinstance(other_names, str):
        row["other_names"] = {s.strip() for s in other_names.split(";") if s.strip()}
    else:
        row["other_names"] = set(other_names)

    mass_per_litre = properties_dict.get("mass_per_litre", 0.0)
    if mass_per_litre is None:
        mass_per_litre = 0.0
    row["mass_per_litre"] = float(mass_per_litre)

    # mmol_per_liter / mol_per_liter override the PubChem-derived calculation entirely.
    # Precedence: mmol_per_liter > mol_per_liter > mass_per_litre → MW lookup.
    mmol_per_liter = properties_dict.get("mmol_per_liter")
    mol_per_liter = properties_dict.get("mol_per_liter")
    if mmol_per_liter is not None:
        row["mmol_concentration"] = float(mmol_per_liter)
    elif mol_per_liter is not None:
        row["mmol_concentration"] = float(mol_per_liter) * 1000.0
    else:
        pubchem_name = properties_dict.get("pubchem_name")
        pubchem_id = properties_dict.get("pubchem_id")
        if pubchem_name is not None:
            molecular_weight = find_molecular_weight(str(pubchem_name))
        elif pubchem_id is not None:
            molecular_weight = find_molecular_weight_by_id(str(pubchem_id))
        elif row["iupac_name"] != "":
            molecular_weight = find_molecular_weight(row["iupac_name"])
        else:
            molecular_weight = find_molecular_weight(row["name"])

        mol = g_to_mol(float(mass_per_litre), molecular_weight) if molecular_weight else 0.0
        row["mmol_concentration"] = mol * 1000.0

    if "flux_upper_bound" in properties_dict:
        row["flux_upper_bound"] = float(properties_dict["flux_upper_bound"])
    else:
        row["flux_upper_bound"] = 0.0

    if "source" in properties_dict:
        row["source"] = str(properties_dict["source"])
    else:
        row["source"] = MediumSource.UNSTATED.value

    row["record_origin"] = RecordOrigin.UPDATED.value
    return row


def _explode_compound(
    df_mapping: pd.DataFrame,
    compound_name: str,
    properties: dict[str, Any],
    halt_on_error: bool,
    verbose: bool,
    sbml_lookup: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Split a compound into per-subcompound mapping rows weighted by molecular weight."""
    compound_name = compound_name.strip().lower()
    if len(properties.get("compounds", [])) == 0:
        raise ValueError(
            f"The exchange name:'{compound_name}' is marked as compound but no compounds were provided."
        )
    compounds = []
    sbml_lookup = sbml_lookup or {}
    for subcompound in properties.get("compounds", []):
        subcompound_name = str(subcompound).strip().lower()
        df = df_mapping.loc[
            (df_mapping["name"] == subcompound_name) | (df_mapping["iupac_name"] == subcompound_name), :
        ]
        if df.empty:
            candidate_names = [str(v) for v in df_mapping["name"].fillna("").values if str(v).strip()]
            candidate_iupac = [str(v) for v in df_mapping["iupac_name"].fillna("").values if str(v).strip()]
            fuzzy_pool = candidate_names + candidate_iupac
            matches = difflib.get_close_matches(subcompound_name, fuzzy_pool, n=1, cutoff=0.8)
            if matches:
                df = df_mapping.loc[
                    (df_mapping["name"] == matches[0]) | (df_mapping["iupac_name"] == matches[0]), :
                ]

        if df.empty:
            normalized = _normalize_mapping_key(subcompound_name)
            sbml_exchange = sbml_lookup.get(normalized)
            if sbml_exchange is not None:
                try:
                    compounds.append({
                        "compound": get_compound_by_name(subcompound_name),
                        "mapping_row": pd.Series({
                            "name": subcompound_name,
                            "iupac_name": "",
                            "other_names": set(),
                            "exchange_reaction": sbml_exchange,
                            "mass_per_litre": 0.0,
                        }),
                    })
                    logging.info(
                        "Resolved compound component '%s' via SBML exchange '%s'.",
                        subcompound_name,
                        sbml_exchange,
                    )
                    continue
                except Exception as exc:
                    logging.info(
                        "SBML component fallback for '%s' found exchange '%s' but compound lookup failed: %s",
                        subcompound_name,
                        sbml_exchange,
                        exc,
                    )
            _log_alert(
                f"Could not resolve compound component '{subcompound_name}' for compound '{compound_name}' "
                "after checking mapping rows and SBML fallbacks.",
                halt_on_error=halt_on_error,
                verbose=verbose,
            )
            continue
        if df.shape[0] > 1:
            _log_alert(
                f"Multiple entries found for compound name:'{subcompound}' for compound '{compound_name}'. Using the first one.",
                halt_on_error=False,
                verbose=verbose,
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

    total_molar_mass = sum(c["compound"].molecular_weight for c in compounds)
    rows = []
    for entry in compounds:
        compound = entry["compound"]
        existing_row = entry["mapping_row"]
        ratio = compound.molecular_weight / total_molar_mass if total_molar_mass > 0 else 0.0
        mass_per_litre = properties.get("mass_per_litre", 0.0)
        corrected_mass_per_litre = float(mass_per_litre) * ratio + existing_row.get("mass_per_litre", 0.0)
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
        rows.append(_build_mapping_row(existing_row["name"], row_properties))
    return rows


def _apply_input_mapping(
    input_mapping: dict[str, Any],
    df_mapping: pd.DataFrame,
    fuzzy_threshold: float,
    halt_on_error: bool,
    verbose: bool,
    sbml_lookup: dict[str, str] | None = None,
    sbml_exchange_ids: set[str] | None = None,
    resolution_stats: MappingResolutionStats | None = None,
) -> pd.DataFrame:
    """Iterate an input mapping (dict) and upsert each entry into df_mapping."""
    for name, properties in input_mapping.items():
        if name is None:
            continue
        if isinstance(properties, dict) and properties.get("is_compound", False):
            rows = _explode_compound(
                df_mapping,
                name,
                properties,
                halt_on_error,
                verbose,
                sbml_lookup=sbml_lookup,
            )
            for row in rows:
                df_mapping = _upsert_mapping_row(
                    df_mapping,
                    row,
                    halt_on_error,
                    fuzzy_threshold,
                    verbose,
                    sbml_lookup=sbml_lookup,
                    sbml_exchange_ids=sbml_exchange_ids,
                    resolution_stats=resolution_stats,
                )
        else:
            row = _build_mapping_row(name, properties)
            df_mapping = _upsert_mapping_row(
                df_mapping,
                row,
                halt_on_error,
                fuzzy_threshold,
                verbose,
                sbml_lookup=sbml_lookup,
                sbml_exchange_ids=sbml_exchange_ids,
                resolution_stats=resolution_stats,
            )
    return df_mapping


def _load_yaml_mappings(custom_mapping_file: str | Path | None) -> dict[str, Any]:
    """Resolve a yaml file (default if None) and load its `custom_mapping` block."""
    if custom_mapping_file is None:
        default_yaml_path = Path(__file__).parent / "yamls" / "custom_exchange_mapping.yaml"
        if default_yaml_path.exists():
            custom_mapping_file = default_yaml_path

    if custom_mapping_file is None or not str(custom_mapping_file).strip():
        return {}
    try:
        import yaml

        yaml_path = Path(__file__).parent / "yamls" / custom_mapping_file
        if not yaml_path.exists():
            return {}
        with open(yaml_path, encoding="utf-8") as f:
            yaml_data = yaml.safe_load(f)
        if yaml_data and "custom_mapping" in yaml_data:
            return yaml_data["custom_mapping"] or {}
    except Exception as e:
        logging.warning(f"Failed to load custom mapping file {custom_mapping_file}: {e}")
    return {}


# ---------------------------------------------------------------------------
# Private helpers — flux/inputs precedence
# ---------------------------------------------------------------------------


def _split_by_source(
    mappings_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Filter `mappings_df` into (valid, supplement, medium, fixed_with_bounds) subsets.

    Preserves the existing semantics: valid = SUPPLEMENT|MEDIUM|FIXED (excludes UNSTATED).
    `fixed_with_bounds` includes `flux_upper_bound` for downstream use.
    """
    valid_mappings = mappings_df.loc[
        (
            (mappings_df["source"] == MediumSource.SUPPLEMENT.value)
            | (mappings_df["source"] == MediumSource.MEDIUM.value)
            | (mappings_df["source"] == MediumSource.FIXED.value)
        )
        & (mappings_df["source"] != MediumSource.UNSTATED.value),
        :,
    ]
    supplement_mappings = valid_mappings.loc[valid_mappings["source"] == MediumSource.SUPPLEMENT.value, :]
    medium_mappings = valid_mappings.loc[valid_mappings["source"] == MediumSource.MEDIUM.value, :]
    fixed_exchanges = mappings_df.loc[
        mappings_df["source"] == MediumSource.FIXED.value, ["exchange_reaction", "flux_upper_bound"]
    ]
    return valid_mappings, supplement_mappings, medium_mappings, fixed_exchanges


def _calculate_flux_row(
    mapping_row: pd.DataFrame,
    growth_row: pd.Series,
    max_time_column: str,
    max_growth_rate_column: str,
    index: int = 0,
) -> tuple[str, float, float]:
    """Return `(exchange_reaction, flux_upper_bound, flux)` for one mapping row.

    Flux formula: `mmol / (max_time * od600_at_gr)` (matches the lab pipeline).
    """
    exchange_reaction = mapping_row["exchange_reaction"].values[index]
    mmol_value = mapping_row["mmol_concentration"].values[index]
    flux_upper_bound = mapping_row["flux_upper_bound"].values[index]
    max_time = growth_row.get(max_time_column, 0.0)
    od600_at_gr = growth_row.get(max_growth_rate_column, 0.0)
    flux = (
        mmol_value / (max_time * od600_at_gr)
        if max_time > 0 and od600_at_gr > 0
        else 0.0
    )
    return exchange_reaction, flux_upper_bound, flux


# ---------------------------------------------------------------------------
# Private helpers — SBML XML fallbacks
# ---------------------------------------------------------------------------


def _parse_sbml_exchanges_fallback(sbml_path: str | Path) -> dict[str, str]:
    """XML-only fallback when COBRApy is missing. Same EX_ filter as the cobra path."""
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
