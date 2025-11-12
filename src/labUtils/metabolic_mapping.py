##############################################################
#  labUtils: Metabolic network mapping utilities             #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2025                          #
#  License: MIT                                              #
#                                                            #
#  Example usage:                                            #
#  From growth rates dataframe                               #
#                                                            #
# from labUtils.metabolic_mapping import (                   #
#     create_supplement_exchange_matrix,                     #
#     parse_sbml_exchanges,                                  #
#     load_minimal_media_exchanges                           #
# )                                                          #
#                                                            #
# # Parse SBML model first                                   #
# mapping = parse_sbml_exchanges("path/to/iML1515.xml")      #
#                                                            #
# # Create exchange matrix                                   #
# df = ... # output of fit_modified_gompertz_per_series      #
# matrix = create_supplement_exchange_matrix(                #
#     df,                                                    #
#     supplement_to_exchange_map=mapping,                    #
#     supplement_column="supplements",                       #
#     growth_rate_column="mu_max",                           #
#     baseline_exchanges=load_minimal_media_exchanges(),     #
#     separator=";",                                         #
#     fuzzy_threshold=0.6                                    #
# )                                                          #
#                                                            #
# Or use curated mapping without SBML:                       #
#                                                            #
# from labUtils.metabolic_mapping import (                   #
#     load_default_iml1515_mapping                           #
# )                                                          #
#                                                            #
# mapping = load_default_iml1515_mapping()                   #
# matrix = create_supplement_exchange_matrix(                #
#     df,                                                    #
#     supplement_to_exchange_map=mapping,                    #
#     baseline_exchanges=load_minimal_media_exchanges()      #
# )                                                          #
##############################################################

import difflib
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

# pylint: disable=import-outside-toplevel


def create_supplement_exchange_matrix(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: Dict[str, str],
    supplement_column: str = "supplements",
    growth_rate_column: str = "mu_max",
    success_column: str = "success",
    baseline_exchanges: Optional[List[str]] = None,
    custom_mapping: Optional[Dict[str, str]] = None,
    separator: str = ";",
    fuzzy_threshold: float = 0.6,
) -> pd.DataFrame:
    """
    Create a binary matrix mapping supplements to exchange reactions.

    This function takes a growth rates dataframe (output from fit_modified_gompertz_per_series)
    and creates a matrix where:
    - Each row corresponds to a row in the input dataframe
    - Each column corresponds to an exchange reaction in the metabolic model
    - Values are 1 if the metabolite is present in supplements, 0 otherwise
    - The last column contains the growth rate

    Parameters
    ----------
    growth_rates_df : pd.DataFrame
        Output from fit_modified_gompertz_per_series containing growth rates and metadata
    supplement_to_exchange_map : Dict[str, str]
        Mapping from supplement names to exchange reaction IDs.
        Example: {"glucose": "EX_glc__D_e", "ribose": "EX_rib__D_e_i"}
        Use parse_sbml_exchanges() to extract this from an SBML file,
        or use load_default_iml1515_mapping() for a curated E. coli mapping.
    supplement_column : str
        Name of the column containing supplement information (default: "supplements")
    success_column : str
        Name of the column containing success status (default: "success")
    growth_rate_column : str
        Name of the column containing growth rate values (default: "mu_max")
    baseline_exchanges : Optional[List[str]]
        List of exchange reactions that are always present (minimal media components).
        These will be set to 1 for all rows.
    custom_mapping : Optional[Dict[str, str]]
        User-defined mapping that takes precedence over supplement_to_exchange_map.
        Use this to override specific mappings or add custom supplement names.
        Example: {"my_glucose": "EX_glc__D_e", "special_carbon": "EX_custom_e"}
    separator : str
        Character used to separate multiple supplements in the supplement column (default: ";")
    fuzzy_threshold : float
        Threshold for fuzzy matching of supplement names (default: 0.6).
        Higher values require closer matches.

    Returns
    -------
    pd.DataFrame
        Matrix with exchange reactions as columns and growth rate as the last column

    Examples
    --------
    >>> # Parse SBML file first
    >>> mapping = parse_sbml_exchanges("iML1515.xml")
    >>> # Or use default mapping
    >>> mapping = load_default_iml1515_mapping()
    >>>
    >>> # Create matrix
    >>> baseline = load_minimal_media_exchanges()
    >>> matrix = create_supplement_exchange_matrix(
    ...     growth_rates_df,
    ...     supplement_to_exchange_map=mapping,
    ...     baseline_exchanges=baseline
    ... )
    >>>
    >>> # With custom overrides
    >>> custom = {"weird_glucose": "EX_glc__D_e"}
    >>> matrix = create_supplement_exchange_matrix(
    ...     growth_rates_df,
    ...     supplement_to_exchange_map=mapping,
    ...     custom_mapping=custom,
    ...     baseline_exchanges=baseline
    ... )

    Notes
    -----
    - Supplement names are case-insensitive and whitespace is stripped
    - Unique supplements are extracted from the dataframe and mapped automatically
    - If exact match fails, fuzzy matching is attempted with configurable threshold
    - Unmapped supplements trigger a warning with details
    - Baseline exchanges are always set to 1 for all rows
    - Growth rate column name is preserved as specified in growth_rate_column parameter
    - Custom mapping takes precedence over supplement_to_exchange_map for matching
    """
    # Normalize mapping keys (lowercase, stripped)
    # Start with base mapping, then apply custom mapping (which takes precedence)
    normalized_map: Dict[str, str] = {}
    for k, v in supplement_to_exchange_map.items():
        if k is None:
            continue
        normalized_map[k.strip().lower()] = v

    # Apply custom mapping overrides
    if custom_mapping:
        for k, v in custom_mapping.items():
            if k is None:
                continue
            normalized_map[k.strip().lower()] = v

    # Collect unique supplements from the dataframe
    unique_supplements = set()
    if supplement_column in growth_rates_df.columns:
        for val in growth_rates_df[supplement_column].dropna().astype(str):
            parts = [s.strip() for s in val.split(separator) if s.strip()]
            for p in parts:
                unique_supplements.add(p.lower())

    # Map supplements -> exchange reaction ids using exact then fuzzy matching
    supplement_to_rxn: Dict[str, str] = {}
    unmapped: List[str] = []
    map_keys = list(normalized_map.keys())
    for supp in sorted(unique_supplements):
        if supp in normalized_map:
            supplement_to_rxn[supp] = normalized_map[supp]
            continue

        # Fuzzy match against SBML/manual keys
        if map_keys:
            matches = difflib.get_close_matches(supp, map_keys, n=1, cutoff=fuzzy_threshold)
            if matches:
                best = matches[0]
                supplement_to_rxn[supp] = normalized_map[best]
                warnings.warn(f"Supplement '{supp}' fuzzy-matched to '{best}' -> {normalized_map[best]}")
                continue

        unmapped.append(supp)

    if unmapped:
        warnings.warn(f"The following supplements could not be mapped to exchange reactions: {unmapped}")

    # Define all exchange columns: mapped reactions + baseline_exchanges
    mapped_rxns = sorted(set(supplement_to_rxn.values()))
    all_exchanges = sorted(set(mapped_rxns + (baseline_exchanges or [])))

    # Build rows
    rows: List[Dict[str, Union[int, float, None]]] = []
    for _, row in growth_rates_df.iterrows():
        # Skip rows where success column exists and is False
        if success_column in growth_rates_df.columns and not row[success_column]:
            continue

        r: Dict[str, Union[int, float, None]] = {ex: 0 for ex in all_exchanges}

        # baseline -> always 1 if provided
        if baseline_exchanges:
            for ex in baseline_exchanges:
                if ex in r:
                    r[ex] = 1
                else:
                    # if baseline exchange not yet in columns, add it
                    r[ex] = 1

        # parse supplements for this row
        supps_val = row.get(supplement_column, None)
        if pd.notna(supps_val):
            for part in [s.strip().lower() for s in str(supps_val).split(separator) if s.strip()]:
                rxn = supplement_to_rxn.get(part)
                if rxn:
                    # ensure column exists
                    if rxn not in r:
                        r[rxn] = 0
                    r[rxn] = 1

        # add growth rate as final column (preserve column name passed)
        r[growth_rate_column] = row.get(growth_rate_column, None)

        rows.append(r)

    result_df = pd.DataFrame(rows)

    # ensure baseline columns are present (even if empty)
    if baseline_exchanges:
        for ex in baseline_exchanges:
            if ex not in result_df.columns:
                result_df[ex] = 1

    # Reorder columns: exchanges sorted, then growth rate last
    exch_cols = sorted([c for c in result_df.columns if c != growth_rate_column])
    result_df = result_df[exch_cols + [growth_rate_column]]

    return result_df


def get_supplement_mapping(
    growth_rates_df: pd.DataFrame,
    supplement_to_exchange_map: Dict[str, str],
    supplement_column: str = "supplements",
    separator: str = ";",
    fuzzy_threshold: float = 0.6,
) -> Dict[str, Any]:
    """
    Get the supplement->reaction mapping for auditing purposes.

    This utility function returns the mapping used by create_supplement_exchange_matrix,
    including information about exact matches, fuzzy matches, and unmapped supplements.

    Parameters
    ----------
    growth_rates_df : pd.DataFrame
        Growth rates dataframe containing supplement information
    supplement_to_exchange_map : Dict[str, str]
        Mapping from supplement names to exchange reaction IDs
        Use parse_sbml_exchanges() to extract this from an SBML file
    supplement_column : str
        Name of the column containing supplement information (default: "supplements")
    separator : str
        Character used to separate multiple supplements (default: ";")
    fuzzy_threshold : float
        Threshold for fuzzy matching (default: 0.6)

    Returns
    -------
    Dict[str, Dict[str, str]]
        Dictionary with three keys:
        - "exact_matches": Dict mapping supplements to reactions (exact matches)
        - "fuzzy_matches": Dict mapping supplements to reactions with match info
        - "unmapped": List of supplements that couldn't be mapped

    Examples
    --------
    >>> mapping = parse_sbml_exchanges("iML1515.xml")
    >>> mapping_info = get_supplement_mapping(growth_df, mapping)
    >>> print("Exact matches:", mapping_info["exact_matches"])
    >>> print("Fuzzy matches:", mapping_info["fuzzy_matches"])
    >>> print("Unmapped:", mapping_info["unmapped"])
    """
    # Normalize mapping keys (lowercase, stripped)
    normalized_map: Dict[str, str] = {}
    for k, v in supplement_to_exchange_map.items():
        if k is None:
            continue
        normalized_map[k.strip().lower()] = v

    # Collect unique supplements from the dataframe
    unique_supplements = set()
    if supplement_column in growth_rates_df.columns:
        for val in growth_rates_df[supplement_column].dropna().astype(str):
            parts = [s.strip() for s in val.split(separator) if s.strip()]
            for p in parts:
                unique_supplements.add(p.lower())

    # Map supplements -> exchange reaction ids
    exact_matches: Dict[str, str] = {}
    fuzzy_matches: Dict[str, str] = {}
    unmapped: List[str] = []

    map_keys = list(normalized_map.keys())

    for supp in sorted(unique_supplements):
        if supp in normalized_map:
            exact_matches[supp] = normalized_map[supp]
            continue

        # Fuzzy match
        if map_keys:
            matches = difflib.get_close_matches(supp, map_keys, n=1, cutoff=fuzzy_threshold)
            if matches:
                best = matches[0]
                fuzzy_matches[supp] = f"{normalized_map[best]} (matched from '{best}')"
                continue

        unmapped.append(supp)

    return {
        "exact_matches": exact_matches,
        "fuzzy_matches": fuzzy_matches,
        "unmapped": unmapped
    }


def parse_sbml_exchanges(sbml_path: Union[str, Path]) -> Dict[str, str]:
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


def _parse_sbml_exchanges_fallback(sbml_path: Union[str, Path]) -> Dict[str, str]:
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


def create_exchange_bounds_template(
    exchange_matrix: pd.DataFrame,
    growth_rate_column: str = "mu_max",
    default_level: int = 1,
    default_max_value: int = 1000,
    sbml_bounds: Optional[Dict[str, Tuple[int, int]]] = None,
    custom_bounds: Optional[Dict[str, Tuple[int, int]]] = None,
) -> pd.DataFrame:
    """
    Create a template dataframe for exchange reaction bounds configuration.

    This function generates a configuration template with predefined values for
    level, max_value, and ratio_drawing parameters for all exchange reactions
    in the input matrix (excluding the growth rate column).

    Parameters
    ----------
    exchange_matrix : pd.DataFrame
        Output from create_supplement_exchange_matrix containing exchange reactions
        as columns and growth rate as the last column
    growth_rate_column : str
        Name of the growth rate column to exclude (default: "mu_max")
    default_level : int
        Default value for the "level" row (default: 1)
    default_max_value : int
        Default value for the "max_value" row (default: 1000)
    sbml_bounds : Optional[Dict[str, Tuple[int, int]]]
        Bounds extracted from SBML file using parse_sbml_exchange_bounds().
        Dictionary mapping exchange reaction IDs to (level, max_value) tuples.
        Example: {"EX_glc__D_e": (1, 500), "EX_o2_e_i": (1, 2000)}
        These are used if a column is not in custom_bounds.
    custom_bounds : Optional[Dict[str, Tuple[int, int]]]
        Custom bounds for specific exchange reactions (highest priority).
        Dictionary mapping exchange reaction IDs to (level, max_value) tuples.
        Example: {"EX_glc__D_e": (2, 500), "EX_o2_e_i": (1, 2000)}
        Takes precedence over sbml_bounds and defaults.

    Returns
    -------
    pd.DataFrame
        Configuration dataframe with columns:
        - "name": Row identifier ("level", "max_value", "ratio_drawing")
        - Exchange reaction columns: Values for each reaction

    Examples
    --------
    >>> # Create exchange matrix
    >>> matrix = create_supplement_exchange_matrix(growth_df, mapping)
    >>>
    >>> # Scenario 1: Generate bounds template with defaults
    >>> bounds = create_exchange_bounds_template(matrix)
    >>> print(bounds)
             name  EX_glc__D_e  EX_o2_e_i  ...
    0       level            1          1  ...
    1   max_value         1000       1000  ...
    2  ratio_drawing         0          0  ...
    >>>
    >>> # Scenario 2: With custom defaults
    >>> bounds = create_exchange_bounds_template(
    ...     matrix,
    ...     default_level=2,
    ...     default_max_value=500
    ... )
    >>>
    >>> # Scenario 3: Use SBML bounds
    >>> sbml_bounds = parse_sbml_exchange_bounds("iML1515.xml")
    >>> bounds = create_exchange_bounds_template(
    ...     matrix,
    ...     sbml_bounds=sbml_bounds
    ... )
    >>>
    >>> # Scenario 4: SBML bounds + custom overrides
    >>> sbml_bounds = parse_sbml_exchange_bounds("iML1515.xml")
    >>> custom = {"EX_glc__D_e": (5, 200)}
    >>> bounds = create_exchange_bounds_template(
    ...     matrix,
    ...     sbml_bounds=sbml_bounds,
    ...     custom_bounds=custom
    ... )

    Notes
    -----
    Precedence order (highest to lowest):
    1. custom_bounds - User-specified overrides
    2. sbml_bounds - Bounds from SBML file
    3. default_level and default_max_value - Fallback defaults

    Default values:
    - level: 1 (enabled)
    - max_value: 1000 (upper bound)
    - ratio_drawing: 0 (no ratio constraint)

    You can modify these values based on your metabolic model requirements.
    Custom bounds take precedence over SBML bounds, which take precedence over defaults.
    """
    # Get all exchange columns (exclude growth rate column)
    exchange_cols = [col for col in exchange_matrix.columns if col != growth_rate_column]

    # Create the template dataframe
    template_data: Dict[str, List[Union[str, int]]] = {
        "name": ["level", "max_value", "ratio_drawing"]
    }

    # Add columns for each exchange reaction with values
    for col in exchange_cols:
        # Precedence: custom_bounds > sbml_bounds > defaults
        if custom_bounds and col in custom_bounds:
            # Highest priority: custom bounds
            level_val, max_val = custom_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        elif sbml_bounds and col in sbml_bounds:
            # Second priority: SBML bounds
            level_val, max_val = sbml_bounds[col]
            template_data[col] = [level_val, max_val, 0]
        else:
            # Lowest priority: default values
            template_data[col] = [default_level, default_max_value, 0]

    return pd.DataFrame(template_data)


def parse_sbml_exchange_bounds(
    sbml_path: Union[str, Path],
    default_level: int = 1,
) -> Dict[str, Tuple[int, int]]:
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
        model = cobra.io.read_sbml_model(str(sbml_path))

        bounds_map = {}
        for reaction in model.reactions:
            if reaction.id.startswith("EX_"):
                # Get upper bound (max_value)
                upper_bound = reaction.upper_bound

                # Cap very large or infinite bounds at 1000
                if upper_bound > 10000 or upper_bound == float('inf'):
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
    sbml_path: Union[str, Path],
    default_level: int = 1,
) -> Dict[str, Tuple[int, int]]:
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
                            if value > 10000 or value == float('inf'):
                                upper_bound = 1000
                            else:
                                upper_bound = int(value)
                        except (ValueError, TypeError):
                            upper_bound = 1000

            bounds_map[rxn_id] = (default_level, upper_bound)

    return bounds_map


def load_default_iml1515_mapping() -> Dict[str, str]:
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


def load_minimal_media_exchanges() -> List[str]:
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
        "EX_pi_e_i",      # Phosphate
        "EX_co2_e_i",     # Carbon dioxide
        "EX_fe3_e_i",     # Iron (III)
        "EX_h_e_i",       # Proton
        "EX_mn2_e_i",     # Manganese
        "EX_fe2_e_i",     # Iron (II)
        "EX_zn2_e_i",     # Zinc
        "EX_mg2_e_i",     # Magnesium
        "EX_ca2_e_i",     # Calcium
        "EX_ni2_e_i",     # Nickel
        "EX_cu2_e_i",     # Copper
        "EX_sel_e_i",     # Selenium
        "EX_cobalt2_e_i", # Cobalt
        "EX_h2o_e_i",     # Water
        "EX_mobd_e_i",    # Molybdate
        "EX_so4_e_i",     # Sulfate
        "EX_nh4_e_i",     # Ammonium
        "EX_k_e_i",       # Potassium
        "EX_na1_e_i",     # Sodium
        "EX_cl_e_i",      # Chloride
        "EX_o2_e_i",      # Oxygen
        "EX_tungs_e_i",   # Tungsten
        "EX_slnt_e_i",    # Selenite
    ]
