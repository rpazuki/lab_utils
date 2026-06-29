# Exchange Mapping Convention

Reference for `custom_exchange_mapping.yaml` and how concentration values flow through the pipeline into FBA exchange fluxes.

---

## Overview

The exchange mapping connects biological supplement names to SBML exchange reaction IDs and carries the concentration of each compound in the growth medium. It is the central configuration that drives:

- `build_mappings()` — assembles a `mappings_df` DataFrame used throughout the pipeline
- `build_supplement_flux_dataframe()` — converts concentrations to uptake fluxes (mmol/gCDW/h)
- `build_flux_dataframe()` in `synthetic.py` — does the same for in-silico generated conditions

---

## YAML entry format

```yaml
custom_mapping:
    entry_name:                  # key used to match supplement names from experiment metadata
        exchange_name: EX_xxx_e  # SBML exchange reaction ID
        iupac_name: ""           # optional; used for PubChem lookup and fuzzy matching
        other_names: []          # optional aliases
        pubchem_name: ""         # preferred name for PubChem MW lookup (overrides iupac/name)
        pubchem_id: 0            # PubChem CID; used when name-based lookup is unreliable

        # --- Concentration (pick one) ---
        mass_per_litre: 0.0      # g/L; MW fetched from PubChem to derive mmol
        mmol_per_liter: 0.0      # mmol/L; bypasses PubChem lookup entirely
        mol_per_liter:  0.0      # mol/L; bypasses PubChem lookup entirely

        flux_upper_bound: 0.0    # explicit upper bound passed to the FBA solver
        source: unstated         # unstated | supplement | medium | fixed
        record_origin: smbl      # smbl | updated
        is_compound: false       # true if this entry dissolves into multiple ions/molecules
        compounds:               # list of sub-component entry names (only when is_compound: true)
            - sodium
            - chloride
```

---

## Concentration fields

### Precedence (highest → lowest)

| Field | Unit | PubChem call? | Notes |
|---|---|---|---|
| `mmol_per_liter` | mmol/L | No | Stored directly as `mmol_concentration` |
| `mol_per_liter` | mol/L | No | Multiplied by 1000, stored as `mmol_concentration` |
| `mass_per_litre` | g/L | Yes | Divided by MW (g/mol) × 1000 → `mmol_concentration` |

If multiple fields are present, only the highest-priority one takes effect. If none is present, `mmol_concentration` is 0.

`mass_per_litre` is always stored in the `mappings_df` regardless of which concentration field is used — it is kept for reference and for MW back-derivation in `synthetic.py`.

### When to use each field

**`mass_per_litre`** — the natural choice when you know the stock concentration in g/L (the most common lab measurement). Requires an internet-accessible PubChem lookup at build time.

**`mmol_per_liter`** — use when:
- You already know the molar concentration precisely (e.g. from a protocol)
- PubChem lookup would be unreliable or unavailable
- The compound is not in PubChem (e.g. custom conjugates)

**`mol_per_liter`** — same as `mmol_per_liter` but for very small concentrations expressed in mol/L.

### Compound entries (`is_compound: true`)

Compound entries (salts, hydrates, etc.) always use `mass_per_litre`. The total mass is split across sub-components proportionally by molecular weight ratio. `mmol_per_liter` and `mol_per_liter` are **not supported** for compound entries — there is no sensible way to distribute a fixed molar concentration across dissociated ions.

---

## Source field

Controls which precedence tier the entry belongs to and how flux is assigned.

| Value | Meaning | Typical origin |
|---|---|---|
| `unstated` | Default; no flux contribution unless overridden | SBML base map |
| `supplement` | Varies per experiment condition; read from metadata | Experimental design |
| `medium` | Present in every condition at a fixed concentration | Defined growth medium (e.g. M9 salts) |
| `fixed` | Fixed flux, not derived from concentration | Gases, ions held constant by model bounds |

Precedence during flux assignment: **supplement > medium > fixed > unstated**.

---

## How `mmol_concentration` becomes a flux

The flux (mmol/gCDW/h) used by the FBA solver is:

```
flux = mmol_concentration / (max_time × OD_at_max_growth_rate)
```

where:
- `max_time` — duration of the growth experiment (hours)
- `OD_at_max_growth_rate` — OD600 at the point of maximum growth rate (proxy for biomass)

This applies to both `medium` and `supplement` sources. `fixed` sources use `flux_upper_bound` directly.

---

## Mapping layers and precedence

`build_mappings()` builds `mappings_df` in four layers, each overriding the previous:

1. **SBML base** — exchange reactions parsed from the organism model (`supplement_to_exchange_map`)
2. **YAML file** — `custom_exchange_mapping.yaml` (or a custom file path)
3. **`custom_mapping` parameter** — passed directly at call time; highest file-level precedence
4. **Experiment metadata** — unique supplements found in the `growth_rates_df` supplement column; source is set to `supplement`

---

## Synthetic pipeline (`synthetic.py`)

`build_flux_dataframe()` in `synthetic.py` follows the same concentration logic:

- **Medium entries**: `mmol_concentration` from `mappings_df` is used directly (same formula above). All conditions share this fixed concentration.
- **Fixed entries**: `flux_upper_bound` from `mappings_df` is used as a constant template value.
- **Supplement entries**: concentration varies per condition (from `conditions_df` in g/L).
  - If MW can be back-derived from `mass_per_litre / mmol_concentration`, the per-condition g/L value is converted to mmol at runtime.
  - If MW cannot be derived (e.g. only `mmol_per_liter` was specified, `mass_per_litre` is 0), the fixed `mmol_concentration` from the mapping is used for all conditions where that supplement is active. This means every active condition gets the same flux, which is the only sensible interpretation when no per-unit-mass conversion is available.

In addition, synthetic enumeration supplement specs now support direct concentration overrides:

- `mmol_per_liter` — direct mmol/L for that supplement in synthetic conversion.
- `mol_per_liter` — direct mol/L (converted internally to mmol/L).

Precedence for synthetic supplement conversion is the same convention as mappings:
`mmol_per_liter` > `mol_per_liter` > g/L (`levels`, `on_value`, `range`).

---

## Example: uracil at known molar concentration

```yaml
uracil_exact:
    exchange_name: EX_ura_e
    mmol_per_liter: 0.200      # 0.2 mmol/L — no PubChem call needed
    source: supplement
```

## Example: glucose in M9 medium

```yaml
glucose:
    exchange_name: EX_glc__D_e
    mass_per_litre: 4.0        # 4 g/L; MW (~180 g/mol) fetched from PubChem → ~22.2 mmol/L
    source: medium
```

## Example: NaCl as a compound

```yaml
nacl:
    mass_per_litre: 0.5
    source: medium
    is_compound: true
    compounds:
        - sodium    # Na⁺  (MW ~23 g/mol → ~39% of total mass)
        - chloride  # Cl⁻  (MW ~35.5 g/mol → ~61% of total mass)
```

The 0.5 g/L is split by MW ratio into a sodium contribution and a chloride contribution, each added to those ions' existing `mass_per_litre` before the mmol calculation.
