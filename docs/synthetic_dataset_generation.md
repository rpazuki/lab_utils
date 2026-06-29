# Synthetic Dataset Generation

The `labUtils.synthetic` module generates in-silico FBA-ready datasets without
any lab measurements.  It produces flux dataframes in exactly the same shape as
`amn_mappings.build_supplement_flux_dataframe`, so the output slots directly
into `fba_tools.solve_fba` and any downstream sampling code.

---

## Concepts

| Term | Meaning |
|------|---------|
| **Condition** | One experimental medium: a vector of supplement concentrations (g/L). |
| **Enumeration** | The strategy used to generate all conditions (Cartesian product, random sampling, etc.). |
| **Flux table** | A DataFrame with one exchange-reaction column per supplement/medium/fixed source plus a `mu_max` column. One per organism. |
| **Phenotype** | How `mu_max` is populated: left empty, solved via FBA, or computed from a formula. |

---

## Outputs — `SyntheticDataset`

`generate_dataset()` returns a `SyntheticDataset` dataclass:

| Attribute | Description |
|-----------|-------------|
| `conditions_df` | One row per condition; columns are `condition_id`, per-supplement concentrations (g/L), and `supplements` (semicolon-joined active supplement names). |
| `medium_df` | Same as `conditions_df` minus the `supplements` text column; used as the shared pool specification in co-culture setups. |
| `flux_tables` | `dict[organism_name → flux DataFrame]` with exchange-reaction columns and `mu_max`. |
| `mappings` | `dict[organism_name → mappings_df]` (the `build_mappings` output used to build each flux table). |
| `flux_df` | Convenience shortcut for the single-organism case (raises if more than one organism). |

`write_dataset(dataset, out_dir)` writes `conditions.csv`, `medium.csv`, and
`df_flux_<organism>.csv` to `out_dir` and returns the written paths.

---

## YAML Configuration

All behaviour is driven by a single YAML file (or an equivalent Python dict).
The top-level key `synthetic:` is optional — the generator accepts either a
bare mapping or one nested under `synthetic:`.

### Full annotated schema

```yaml
synthetic:

  # ── Organisms ────────────────────────────────────────────────────────────────
  organisms:
    ecoli:                           # arbitrary organism key
      # Path to SBML model — only required when phenotype.mode == 'fba'.
      sbml: /path/to/iML1515.xml

      # Reuse an existing exchange-mapping YAML from src/labUtils/yamls/.
      exchange_mapping: custom_exchange_mapping.yaml

      # Suffix appended to every exchange column (e.g. EX_glc__D_e_i).
      exchange_suffix: _i

      # common-name → exchange-reaction-id. IDs must match your SBML exactly.
      supplement_to_exchange_map:
        glucose: EX_glc__D_e
        uracil:  EX_ura_e
        glycine: EX_gly_e

      # Explicit molecular weights (g/mol).
      # If omitted, MW is back-derived from the mapping (requires pubchempy).
      # Set these for deterministic, offline runs.
      molecular_weights:
        glucose: 180.156
        uracil:  112.0868
        glycine: 75.067

  # ── Enumeration ──────────────────────────────────────────────────────────────
  enumeration:
    # Mode — one of: cartesian | custom | presence_absence | random
    mode: cartesian
    n_samples: 200    # random / LHS only
    seed: 42          # deterministic sampling for random / presence_absence
    max_active: 2     # optional cap on how many supplements are active at once
                      # (presence_absence only)

    supplements:
      glucose:
        levels: [0.0, 2.0, 4.0]              # explicit g/L values
      uracil:
        binary: true                          # shorthand for levels: [0, on_value]
        on_value: 0.0224                      # g/L when "on"
        # Optional direct concentration override used in flux conversion.
        # Precedence: mmol_per_liter > mol_per_liter > g/L fields above.
        mmol_per_liter: 0.2                   # mmol/L (fixed for active conditions)
        # mol_per_liter: 0.0002               # mol/L (converted to mmol/L)
      glycine:
        range: { min: 0.0, max: 0.5, n: 5 }  # n evenly-spaced values

  # ── Phenotype ─────────────────────────────────────────────────────────────────
  phenotype:
    # Mode — one of: empty | fba | formula
    mode: empty
    noise_std: 0.0   # Gaussian noise σ added to mu_max (0 = no noise)
    seed: 42         # noise RNG seed
    # formula: "0.6 * (glucose > 0) + 0.1 * uracil"  # formula mode only

  # ── Kinetics defaults ────────────────────────────────────────────────────────
  growth_kinetics_defaults:
    max_time: 24.0        # hours (experiment duration)
    mv_mu_max_value: 0.5  # OD600 at maximum growth rate

  # ── Output ───────────────────────────────────────────────────────────────────
  output:
    dir: ./synthetic_out
```

### Enumeration modes

| Mode | What it generates |
|------|-------------------|
| `cartesian` | Full Cartesian product of all supplement levels. 3 supplements × 3 levels each → 27 conditions. |
| `custom` | Same as `cartesian` but levels are taken directly from each supplement's `levels` list. |
| `presence_absence` | Every binary on/off combination up to `max_active` active supplements at once. |
| `random` | `n_samples` points sampled via Latin Hypercube (requires `scipy`) or uniform random over each supplement's `range`. |

### Supplement spec keys

| Key | Used by | Meaning |
|-----|---------|---------|
| `levels` | `cartesian`, `custom` | Explicit concentration values (g/L). |
| `binary` | `cartesian` | Shorthand for `levels: [0.0, on_value]`. |
| `on_value` | `presence_absence`, `binary` | Concentration when supplement is "on" (default 1.0 g/L). |
| `range` | `cartesian`, `random` | `{min, max, n}` — linspace for `cartesian`, continuous range for `random`. |
| `mmol_per_liter` | all modes | Direct mmol/L override for supplement flux conversion; overrides g/L conversion when active. |
| `mol_per_liter` | all modes | Direct mol/L override (converted to mmol/L); lower precedence than `mmol_per_liter`. |

Override precedence for supplement concentration conversion is:
`mmol_per_liter` > `mol_per_liter` > g/L (`levels`, `on_value`, `range`).
When a direct mmol/mol override is set, it is applied as a fixed mmol/L value for each active condition.

### Phenotype modes

| Mode | Requirement | Behaviour |
|------|-------------|-----------|
| `empty` | — | `mu_max` remains 0.0 everywhere. |
| `fba` | `sbml` path per organism | Calls `solve_fba` on each flux row; copies `fba_growth` into `mu_max`. |
| `formula` | `formula` string | Evaluates a Python expression against `conditions_df` via `pandas.eval`. |

**Formula safety**: only supplement column names and a math whitelist
(`abs`, `min`, `max`, `sqrt`, `log`, `exp`, `pow`, `sin`, `cos`, `tan`,
`where`, `and`, `or`, `not`) are allowed — no attribute access, no arbitrary
function calls.

---

## Python API

### Minimal single-organism example

```python
from labUtils.synthetic import generate_dataset, write_dataset

dataset = generate_dataset("src/labUtils/yamls/synthetic_pipeline.yaml")

# Access DataFrames
print(dataset.conditions_df.head())
print(dataset.flux_df.head())           # single-organism shortcut

# Persist to disk
paths = write_dataset(dataset, "./out/synthetic")
# Writes: conditions.csv, medium.csv, df_flux_ecoli.csv
```

### Inline dict config (no YAML file)

```python
from labUtils.synthetic import generate_dataset

cfg = {
    "organisms": {
        "ecoli": {
            "supplement_to_exchange_map": {
                "glucose": "EX_glc__D_e",
                "glycine": "EX_gly_e",
            },
            "molecular_weights": {"glucose": 180.156, "glycine": 75.067},
        }
    },
    "enumeration": {
        "mode": "cartesian",
        "supplements": {
            "glucose": {"levels": [0.0, 2.0, 4.0]},
            "glycine": {"binary": True, "on_value": 0.75},
        },
    },
    "phenotype": {"mode": "empty"},
    "growth_kinetics_defaults": {"max_time": 24.0, "mv_mu_max_value": 0.5},
}

dataset = generate_dataset(cfg)
print(dataset.conditions_df)
# condition_id  glucose  glycine  supplements
# cond_000000      0.0     0.00
# cond_000001      0.0     0.75  glycine
# cond_000002      2.0     0.00  glucose
# ...
```

### Formula phenotype

```python
cfg = {
    "organisms": {"ecoli": {"supplement_to_exchange_map": {"glucose": "EX_glc__D_e"}}},
    "enumeration": {
        "mode": "random",
        "n_samples": 50,
        "seed": 0,
        "supplements": {"glucose": {"range": {"min": 0.0, "max": 5.0}}},
    },
    "phenotype": {
        "mode": "formula",
        "formula": "0.5 * (1 - exp(-glucose))",
        "noise_std": 0.01,
        "seed": 99,
    },
}

dataset = generate_dataset(cfg)
print(dataset.flux_df[["EX_glc__D_e", "mu_max"]].head())
```

### FBA phenotype (requires COBRApy + SBML)

```python
cfg = {
    "organisms": {
        "ecoli": {
            "sbml": "/path/to/iML1515.xml",
            "supplement_to_exchange_map": {"glucose": "EX_glc__D_e"},
        }
    },
    "enumeration": {
        "mode": "presence_absence",
        "supplements": {"glucose": {"on_value": 4.0}},
    },
    "phenotype": {"mode": "fba"},
}

dataset = generate_dataset(cfg)
print(dataset.flux_df)  # mu_max column now contains FBA-predicted growth rates
```

### Co-culture / community setup

```python
cfg = {
    "organisms": {
        "ecoli":     {"supplement_to_exchange_map": {"glucose": "EX_glc__D_e"}, "exchange_suffix": "_ec"},
        "bsubtilis": {"supplement_to_exchange_map": {"glucose": "R_EX_glc__D_e"}, "exchange_suffix": "_bs"},
    },
    "enumeration": {
        "mode": "cartesian",
        "supplements": {"glucose": {"levels": [0.0, 2.0, 4.0]}},
    },
    "phenotype": {"mode": "empty"},
}

dataset = generate_dataset(cfg)
print(list(dataset.flux_tables))   # ['ecoli', 'bsubtilis']
print(dataset.flux_tables["ecoli"].columns.tolist())
```

---

## Running as a script

`synthetic.py` can be run from the command line via Python's `-m` flag.  A
small `__main__` block is expected in the file; if it is not present, use the
snippet below directly.

### Using the bundled YAML

```bash
python -c "
from labUtils.synthetic import generate_dataset, write_dataset
ds = generate_dataset('src/labUtils/yamls/synthetic_pipeline.yaml')
write_dataset(ds, './synthetic_out')
print('Wrote', list(write_dataset.__doc__ or 'done'))
"
```

### One-liner for quick exploration

```bash
python - <<'EOF'
from labUtils.synthetic import generate_dataset, write_dataset
import yaml, pathlib

cfg = pathlib.Path("src/labUtils/yamls/synthetic_pipeline.yaml")
ds = generate_dataset(cfg)
paths = write_dataset(ds, "./synthetic_out")
for k, p in paths.items():
    print(f"  {k:20s} → {p}")
EOF
```

### Expected output structure

```
synthetic_out/
├── conditions.csv       # condition_id, glucose, uracil, glycine, supplements
├── medium.csv           # condition_id, glucose, uracil, glycine
└── df_flux_ecoli.csv    # EX_glc__D_e_i, EX_ura_e_i, EX_gly_e_i, mu_max
```

---

## Low-level functions

These are used internally by `generate_dataset` but are also importable for
custom pipelines.

| Function | Purpose |
|----------|---------|
| `enumerate_conditions(spec)` | Produces a `conditions_df` from an enumeration spec dict. |
| `build_flux_dataframe(conditions_df, mappings_df, ...)` | Converts g/L concentrations to mmol/gCDW/h fluxes using the exchange mapping. |
| `attach_phenotype(flux_df, conditions_df, mode, ...)` | Fills the `mu_max` column according to the chosen phenotype mode. |
| `write_dataset(dataset, out_dir)` | Writes all CSV outputs and returns `{key: Path}`. |

---

## Dependencies

| Package | Required for |
|---------|-------------|
| `numpy` | Array math in flux conversion and random sampling. |
| `pandas` | All DataFrame operations. |
| `pyyaml` | Parsing YAML config files. |
| `scipy` | Latin Hypercube sampling in `random` mode (falls back to uniform if absent). |
| `cobra` (COBRApy) | `fba` phenotype mode only. |
| `pubchempy` | MW lookup when `molecular_weights` are not set explicitly. |

All packages except `cobra` and `pubchempy` are standard lab_utils dependencies.
