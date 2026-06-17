##############################################################
#  labUtils: Synthetic dataset generator for in-silico FBA   #
#                                                            #
#  Author: Roozbeh H. Pazuki - 2026                          #
#  License: MIT                                              #
#                                                            #
##############################################################

from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from labUtils.amn_mappings import MediumSource, build_mappings


@dataclass
class SyntheticDataset:
    """Container for in-silico generated FBA-ready data.

    Attributes
    ----------
    conditions_df
        One row per condition with `condition_id` + per-supplement concentration columns (g/L).
        Includes a `supplements` semicolon-joined column for parity with `growth_rates_df`.
    medium_df
        One row per condition keyed by `condition_id`; columns are supplement common names
        (g/L). Acts as the shared medium specification when multiple organisms share a pool.
    flux_tables
        Mapping `organism_name -> flux DataFrame` in the shape produced by
        `amn_mappings.build_supplement_flux_dataframe` (exchange columns + `mu_max`).
    mappings
        Mapping `organism_name -> mappings_df` used to build the corresponding flux table.
    """

    conditions_df: pd.DataFrame
    medium_df: pd.DataFrame
    flux_tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    mappings: dict[str, pd.DataFrame] = field(default_factory=dict)

    @property
    def flux_df(self) -> pd.DataFrame:
        """Convenience accessor for the single-organism case."""
        if len(self.flux_tables) != 1:
            raise ValueError(
                f"flux_df is only defined for single-organism datasets; "
                f"got {len(self.flux_tables)} organisms: {list(self.flux_tables)}"
            )
        return next(iter(self.flux_tables.values()))


# ---------------------------------------------------------------------------
# Enumeration
# ---------------------------------------------------------------------------


def enumerate_conditions(spec: dict[str, Any]) -> pd.DataFrame:
    """Generate condition rows from a per-supplement enumeration spec.

    Parameters
    ----------
    spec
        Dict with keys:

        - `mode`: one of `cartesian`, `custom`, `presence_absence`, `random`.
        - `supplements`: dict mapping supplement name to its sampling spec.
          Each supplement spec accepts:
            * `levels: [v1, v2, ...]` — explicit values (used by `cartesian`, `custom`).
            * `binary: true` — interpreted as `levels: [0.0, mass_when_on]` where
              `mass_when_on` defaults to `on_value` (default 1.0).
            * `range: {min, max, n}` — `n` evenly spaced values; used by `cartesian`
              if no `levels`, or sampled continuously by `random`.
            * `on_value: float` — used by `presence_absence` and as the "on" mass
              for `binary` levels (default 1.0).
        - `n_samples`: int, only for `random`.
        - `seed`: int, for `random` / `presence_absence` determinism (optional).
        - `max_active`: int, optional bound for `presence_absence`.

    Returns
    -------
    pd.DataFrame
        One column per supplement (concentration in g/L); rows are conditions.
        An extra `condition_id` column (zero-padded) and a `supplements`
        semicolon-joined column listing supplements with concentration > 0 are
        added to match the `growth_rates_df` schema downstream.
    """
    mode = spec.get("mode", "cartesian")
    supplements = spec.get("supplements", {}) or {}
    if not supplements:
        return _format_conditions(pd.DataFrame([{}]))

    if mode == "cartesian":
        per_supp = {name: _levels_for(supp_spec) for name, supp_spec in supplements.items()}
        rows = [dict(zip(per_supp.keys(), combo, strict=True)) for combo in itertools.product(*per_supp.values())]
        df = pd.DataFrame(rows, columns=list(per_supp.keys()))

    elif mode == "custom":
        per_supp = {name: list(supp_spec.get("levels", [0.0])) for name, supp_spec in supplements.items()}
        rows = [dict(zip(per_supp.keys(), combo, strict=True)) for combo in itertools.product(*per_supp.values())]
        df = pd.DataFrame(rows, columns=list(per_supp.keys()))

    elif mode == "presence_absence":
        names = list(supplements.keys())
        on_values = {name: float(supplements[name].get("on_value", 1.0)) for name in names}
        max_active = spec.get("max_active", len(names))
        patterns: list[tuple[int, ...]] = []
        for k in range(0, min(max_active, len(names)) + 1):
            for combo in itertools.combinations(range(len(names)), k):
                pattern = [0] * len(names)
                for idx in combo:
                    pattern[idx] = 1
                patterns.append(tuple(pattern))
        rows = [{names[i]: on_values[names[i]] * bit for i, bit in enumerate(pattern)} for pattern in patterns]
        df = pd.DataFrame(rows, columns=names)

    elif mode == "random":
        n = int(spec.get("n_samples", 0))
        seed = spec.get("seed")
        names = list(supplements.keys())
        ranges = []
        for name in names:
            supp_spec = supplements[name]
            rng_spec = supp_spec.get("range")
            if rng_spec is None:
                # Fall back to levels min/max if provided
                lvls = supp_spec.get("levels")
                if lvls:
                    ranges.append((float(min(lvls)), float(max(lvls))))
                else:
                    ranges.append((0.0, float(supp_spec.get("on_value", 1.0))))
            else:
                ranges.append((float(rng_spec["min"]), float(rng_spec["max"])))
        samples = _lhs_or_uniform(n, len(names), seed)
        scaled = np.empty_like(samples)
        for i, (lo, hi) in enumerate(ranges):
            scaled[:, i] = lo + (hi - lo) * samples[:, i]
        df = pd.DataFrame(scaled, columns=names)

    else:
        raise ValueError(f"Unknown enumeration mode: '{mode}'. Expected one of cartesian, custom, presence_absence, random.")

    return _format_conditions(df)


def _levels_for(supp_spec: dict[str, Any]) -> list[float]:
    if "levels" in supp_spec:
        return [float(v) for v in supp_spec["levels"]]
    if supp_spec.get("binary", False):
        return [0.0, float(supp_spec.get("on_value", 1.0))]
    rng_spec = supp_spec.get("range")
    if rng_spec is not None:
        n = int(rng_spec.get("n", 2))
        return list(np.linspace(float(rng_spec["min"]), float(rng_spec["max"]), n))
    return [0.0]


def _lhs_or_uniform(n: int, d: int, seed: int | None) -> np.ndarray:
    """Latin hypercube samples in [0,1]^d via scipy if available, else uniform random."""
    if n <= 0 or d <= 0:
        return np.zeros((max(n, 0), d))
    try:
        from scipy.stats import qmc

        sampler = qmc.LatinHypercube(d=d, seed=seed)
        return sampler.random(n=n)
    except ImportError:
        rng = np.random.default_rng(seed)
        return rng.random(size=(n, d))


def _format_conditions(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index(drop=True).copy()
    df["condition_id"] = [f"cond_{i:06d}" for i in range(len(df))]
    supplement_cols = [c for c in df.columns if c != "condition_id"]

    def _active(row: pd.Series) -> str:
        return ";".join(c for c in supplement_cols if float(row[c]) > 0.0)

    df["supplements"] = df.apply(_active, axis=1) if supplement_cols else ""
    ordered = ["condition_id"] + supplement_cols + ["supplements"]
    return df[ordered]


# ---------------------------------------------------------------------------
# Flux dataframe construction
# ---------------------------------------------------------------------------


def _recover_mw(mass_per_litre: float, mmol: float) -> float | None:
    """Back-compute molecular weight (g/mol) from build_mappings' stored mmol.

    `build_mappings` stores `mmol = (mass_per_litre / MW) * 1000`. Invert that here
    to avoid duplicate PubChem lookups in the synthetic pipeline.
    """
    if mass_per_litre and mmol and mass_per_litre > 0 and mmol > 0:
        return (mass_per_litre / mmol) * 1000.0
    return None


def _flux_from_concentration(mass_per_litre: float, molecular_weight: float | None, max_time: float, od_at_gr: float) -> float:
    """Convert concentration (g/L) to uptake flux (mmol/gCDW/h).

    Matches the formula in `amn_mappings.build_supplement_flux_dataframe`
    (mmol_concentration / (max_time * od_at_gr)).
    """
    if not molecular_weight or molecular_weight <= 0 or max_time <= 0 or od_at_gr <= 0:
        return 0.0
    mmol = (mass_per_litre / molecular_weight) * 1000.0
    return mmol / (max_time * od_at_gr)


def build_flux_dataframe(
    conditions_df: pd.DataFrame,
    mappings_df: pd.DataFrame,
    kinetics_defaults: dict[str, float] | None = None,
    exchange_suffix: str | None = None,
    growth_rate_column: str = "mu_max",
    molecular_weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Build a flux dataframe in the shape of `build_supplement_flux_dataframe`.

    The output has one column per valid exchange (SUPPLEMENT/MEDIUM/FIXED source)
    plus `mu_max`, matching what `fba_tools.solve_fba` expects.

    Parameters
    ----------
    conditions_df
        Output of `enumerate_conditions`. Each non-meta column is a supplement
        common name; values are concentrations in g/L.
    mappings_df
        Output of `amn_mappings.build_mappings`. Provides the
        common-name → exchange-reaction mapping, fixed-source bounds, and stored
        mmol/mass for back-computing molecular weights.
    kinetics_defaults
        Defaults for the synthetic growth kinetics:
        `max_time` (h), `mv_mu_max_value` (OD at max growth rate). Both default to
        24.0 and 0.5 respectively.
    exchange_suffix
        Optional suffix appended to every exchange column name (mirrors
        `build_supplement_flux_dataframe`).
    growth_rate_column
        Name of the growth-rate column emitted in the output. Defaults to `mu_max`.
    molecular_weights
        Optional explicit MW (g/mol) per supplement common name. Takes precedence
        over MW back-derived from the mapping's `mmol_concentration` field. Use
        when `pubchempy` is unavailable or you want deterministic MW values.
    """
    mw_overrides = {k.lower(): float(v) for k, v in (molecular_weights or {}).items()}
    kin = {"max_time": 24.0, "mv_mu_max_value": 0.5, **(kinetics_defaults or {})}
    max_time = float(kin["max_time"])
    od_at_gr = float(kin["mv_mu_max_value"])

    valid_sources = {MediumSource.SUPPLEMENT.value, MediumSource.MEDIUM.value, MediumSource.FIXED.value}
    valid_mappings = mappings_df.loc[mappings_df["source"].isin(valid_sources)]

    fixed_mappings = valid_mappings.loc[valid_mappings["source"] == MediumSource.FIXED.value]
    medium_mappings = valid_mappings.loc[valid_mappings["source"] == MediumSource.MEDIUM.value]
    supp_mappings = valid_mappings.loc[valid_mappings["source"] == MediumSource.SUPPLEMENT.value]

    # Pre-compute per-supplement MW (override > back-derived from mapping).
    # mmol_concentration is also stored so that entries with a direct mmol_per_liter /
    # mol_per_liter override (where mass_per_litre may be 0 and MW cannot be derived)
    # can still produce a non-zero flux.
    supp_meta: dict[str, dict[str, Any]] = {}
    for _, mrow in supp_mappings.iterrows():
        name = str(mrow["name"]).lower()
        mw = mw_overrides.get(name) or _recover_mw(
            float(mrow.get("mass_per_litre", 0.0)),
            float(mrow.get("mmol_concentration", 0.0)),
        )
        supp_meta[name] = {
            "exchange_reaction": mrow["exchange_reaction"],
            "molecular_weight": mw,
            "mmol_concentration": float(mrow.get("mmol_concentration", 0.0)),
        }
        if mrow["iupac_name"]:
            supp_meta[str(mrow["iupac_name"]).lower()] = supp_meta[name]

    all_exchanges = list(valid_mappings["exchange_reaction"])
    # Build initial template with precedence FIXED -> MEDIUM -> SUPPLEMENT (supplement
    # is only set per-row when its concentration > 0; otherwise the medium / fixed
    # value carries over).
    template: dict[str, float] = dict.fromkeys(all_exchanges, 0.0)
    for _, mrow in fixed_mappings.iterrows():
        template[mrow["exchange_reaction"]] = float(mrow.get("flux_upper_bound", 0.0))
    for _, mrow in medium_mappings.iterrows():
        # Medium flux uses pre-set mmol from mapping (mass_per_litre × MW captured in
        # `mmol_concentration` by `build_mappings`).
        mmol = float(mrow.get("mmol_concentration", 0.0))
        flux = mmol / (max_time * od_at_gr) if max_time > 0 and od_at_gr > 0 else 0.0
        template[mrow["exchange_reaction"]] = flux

    supplement_cols = [c for c in conditions_df.columns if c not in {"condition_id", "supplements"}]

    rows: list[dict[str, float]] = []
    for _, cond in conditions_df.iterrows():
        new_row = template.copy()
        for supp_name in supplement_cols:
            mass = float(cond.get(supp_name, 0.0))
            if mass <= 0.0:
                continue
            meta = supp_meta.get(supp_name.lower())
            if meta is None:
                # Not in mapping as a supplement source — silently skip; this can be a
                # supplement column that maps onto a MEDIUM or FIXED exchange and is
                # therefore already accounted for via the template.
                continue
            if meta.get("molecular_weight") is None and meta.get("mmol_concentration", 0.0) > 0.0:
                # MW unavailable but a direct mmol override was set (mmol_per_liter /
                # mol_per_liter in the mapping).  Use the fixed mmol for any active condition.
                mmol = meta["mmol_concentration"]
                flux = mmol / (max_time * od_at_gr) if max_time > 0 and od_at_gr > 0 else 0.0
            else:
                flux = _flux_from_concentration(mass, meta["molecular_weight"], max_time, od_at_gr)
            new_row[meta["exchange_reaction"]] = flux
        rows.append(new_row)

    flux_df = pd.DataFrame(rows, columns=all_exchanges)
    flux_df = flux_df.loc[:, ~flux_df.columns.duplicated()].copy()

    exch_cols = sorted(flux_df.columns)
    if exchange_suffix:
        flux_df = flux_df.rename(columns={c: f"{c}{exchange_suffix}" for c in exch_cols})
        exch_cols = [f"{c}{exchange_suffix}" for c in exch_cols]

    flux_df[growth_rate_column] = 0.0
    return flux_df[exch_cols + [growth_rate_column]]


# ---------------------------------------------------------------------------
# Phenotype attachment
# ---------------------------------------------------------------------------


_FORMULA_WHITELIST = {
    "abs", "min", "max", "sqrt", "log", "exp", "pow",
    "sin", "cos", "tan", "where", "and", "or", "not",
    "True", "False",
}


def attach_phenotype(
    flux_df: pd.DataFrame,
    conditions_df: pd.DataFrame,
    mode: str = "empty",
    *,
    model: Any = None,
    medium_columns: list[str] | None = None,
    formula: str | None = None,
    noise_std: float = 0.0,
    seed: int | None = None,
    growth_rate_column: str = "mu_max",
) -> pd.DataFrame:
    """Populate the growth-rate column of `flux_df` according to `mode`.

    Modes:

    - `empty` — leave `mu_max` as 0.0 (the generator default).
    - `fba` — call `fba_tools.solve_fba` on a copy of `flux_df` and copy the
      resulting `fba_growth` column into `mu_max`. Requires `model` and
      `medium_columns`.
    - `formula` — evaluate `formula` via `pandas.DataFrame.eval` against
      `conditions_df`. Allowed names are restricted to supplement column names
      and a small math whitelist.
    """
    out = flux_df.copy()

    if mode == "empty":
        return out

    if mode == "fba":
        if model is None or medium_columns is None:
            raise ValueError("phenotype mode='fba' requires `model` and `medium_columns`.")
        from labUtils.fba.fba_tools import solve_fba

        solved = solve_fba(out, model, list(medium_columns), solution_column_name="_phenotype_growth")
        out[growth_rate_column] = solved["_phenotype_growth"].values
        if noise_std > 0:
            rng = np.random.default_rng(seed)
            out[growth_rate_column] = out[growth_rate_column] + rng.normal(0.0, noise_std, size=len(out))
        return out

    if mode == "formula":
        if not formula:
            raise ValueError("phenotype mode='formula' requires a non-empty `formula`.")
        _validate_formula(formula, allowed_names=set(conditions_df.columns) | _FORMULA_WHITELIST)
        result = conditions_df.eval(formula, engine="python")
        if np.isscalar(result):
            out[growth_rate_column] = float(result)
        else:
            out[growth_rate_column] = pd.Series(result).astype(float).values
        if noise_std > 0:
            rng = np.random.default_rng(seed)
            out[growth_rate_column] = out[growth_rate_column] + rng.normal(0.0, noise_std, size=len(out))
        return out

    raise ValueError(f"Unknown phenotype mode: '{mode}'. Expected one of empty, fba, formula.")


def _validate_formula(formula: str, allowed_names: set[str]) -> None:
    import ast

    tree = ast.parse(formula, mode="eval")
    # Check disallowed names first (more informative error for things like __import__).
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            raise ValueError("Attribute access is not permitted in formula.")
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and node.id not in allowed_names:
            raise ValueError(f"Formula references disallowed name: '{node.id}'.")
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if not isinstance(func, ast.Name) or func.id not in _FORMULA_WHITELIST:
                raise ValueError("Function call in formula is restricted to math whitelist.")


# ---------------------------------------------------------------------------
# YAML-driven orchestration
# ---------------------------------------------------------------------------


def _load_config(config: dict | str | Path) -> dict:
    if isinstance(config, (str, Path)):
        with open(config, encoding="utf-8") as f:
            data = yaml.safe_load(f)
    else:
        data = config
    if not isinstance(data, dict):
        raise ValueError("Synthetic config must be a mapping at top level.")
    if "synthetic" in data:
        data = data["synthetic"]
    return data


def _build_organism_mappings(
    organism_cfg: dict[str, Any],
    enumerated_supplements: list[str],
) -> pd.DataFrame:
    """Construct a `mappings_df` for one organism.

    The synthetic pipeline allows two ways:

    1. Provide `supplement_to_exchange_map` directly in the organism config.
    2. Reuse an existing yaml mapping via `exchange_mapping` (file in `yamls/`).

    Either way, `amn_mappings.build_mappings` is called so the schema is identical
    to the rest of the pipeline.
    """
    base_map = organism_cfg.get("supplement_to_exchange_map", {}) or {}
    custom_mapping = organism_cfg.get("custom_mapping")
    mapping_file = organism_cfg.get("exchange_mapping")

    # Provide a stub growth_rates_df so build_mappings has supplement context.
    stub = pd.DataFrame({
        "well": [f"w{i}" for i in range(len(enumerated_supplements) or 1)],
        "supplements": enumerated_supplements + [""] * max(0, 1 - len(enumerated_supplements)),
        "mu_max": [0.0] * (len(enumerated_supplements) or 1),
    })
    return build_mappings(
        stub,
        supplement_to_exchange_map=base_map,
        custom_mapping_file=mapping_file,
        custom_mapping=custom_mapping,
    )


def generate_dataset(config: dict | str | Path) -> SyntheticDataset:
    """Generate a SyntheticDataset from a YAML or dict spec.

    See `src/labUtils/yamls/synthetic_pipeline.yaml` for the schema.
    """
    cfg = _load_config(config)
    enumeration_spec = cfg.get("enumeration") or {}
    conditions_df = enumerate_conditions(enumeration_spec)

    supplement_cols = [c for c in conditions_df.columns if c not in {"condition_id", "supplements"}]
    medium_df = conditions_df[["condition_id", *supplement_cols]].copy()

    kinetics_defaults = cfg.get("growth_kinetics_defaults") or {}
    phenotype_cfg = cfg.get("phenotype") or {"mode": "empty"}

    organisms_cfg = cfg.get("organisms") or {}
    if not organisms_cfg:
        raise ValueError("Synthetic config must declare at least one organism under `organisms`.")

    flux_tables: dict[str, pd.DataFrame] = {}
    mappings_by_org: dict[str, pd.DataFrame] = {}
    for org_name, org_cfg in organisms_cfg.items():
        mappings_df = _build_organism_mappings(org_cfg, supplement_cols)
        flux_df = build_flux_dataframe(
            conditions_df,
            mappings_df,
            kinetics_defaults=kinetics_defaults,
            exchange_suffix=org_cfg.get("exchange_suffix"),
            molecular_weights=org_cfg.get("molecular_weights"),
        )

        pheno_mode = phenotype_cfg.get("mode", "empty")
        model = None
        medium_columns = None
        if pheno_mode == "fba":
            sbml = org_cfg.get("sbml")
            if not sbml:
                raise ValueError(f"phenotype mode='fba' requires `sbml` for organism '{org_name}'.")
            from labUtils.fba.fba_tools import load_model

            model = load_model(str(sbml))
            medium_columns = [c for c in flux_df.columns if c != "mu_max"]
        flux_df = attach_phenotype(
            flux_df,
            conditions_df,
            mode=pheno_mode,
            model=model,
            medium_columns=medium_columns,
            formula=phenotype_cfg.get("formula"),
            noise_std=float(phenotype_cfg.get("noise_std", 0.0)),
            seed=phenotype_cfg.get("seed"),
        )

        flux_tables[org_name] = flux_df
        mappings_by_org[org_name] = mappings_df

    return SyntheticDataset(
        conditions_df=conditions_df,
        medium_df=medium_df,
        flux_tables=flux_tables,
        mappings=mappings_by_org,
    )


def write_dataset(dataset: SyntheticDataset, out_dir: str | Path) -> dict[str, Path]:
    """Write all dataset components to `out_dir`. Returns the written paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    cond_path = out_dir / "conditions.csv"
    dataset.conditions_df.to_csv(cond_path, index=False)
    paths["conditions"] = cond_path
    med_path = out_dir / "medium.csv"
    dataset.medium_df.to_csv(med_path, index=False)
    paths["medium"] = med_path
    for org, df in dataset.flux_tables.items():
        p = out_dir / f"df_flux_{org}.csv"
        df.to_csv(p, index=False)
        paths[org] = p
    return paths
