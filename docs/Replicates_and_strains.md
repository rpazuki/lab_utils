# Replicates and Strains

The two replicate-statistics functions are:

```python
labUtils.media_bot.calculate_replicate_statistics_by_well
labUtils.media_bot.calculate_replicate_statistics_by_custom
```

## 1) Quick Comparison

- `calculate_replicate_statistics_by_well(...)`
    - Uses one global rule for the whole dataset.
    - You set one `direction` and one `sample_size`.

- `calculate_replicate_statistics_by_custom(...)`
    - Uses per-strain (or per-strain-pattern) rules.
    - `strain_pattern` decides how original strain labels are normalized.
    - `custom_rules` decides which normalized strains are processed and how wells are grouped.

## 2) Custom Rule Shape

```python
class CustomReplicateRule(TypedDict, total=False):
        direction: str      # alphabetical/alpha or numerical/num
        pattern: str        # regex to match strain_group values
        sample_size: int    # fixed replicate block size when provided
```

Defaults used by the function:

- `direction` defaults to `"alphabetical"` if missing.
- `sample_size` defaults to `3` for fixed-size grouping mode.

## 3) How `strain_pattern` Affects Matching

`strain_pattern` is applied first to `df_parsed["strain"]` and creates an internal `strain_group` column.

- If `strain_pattern is None`:
    - `strain_group = strain` (original values preserved).

- If `strain_pattern` is a regex string:
    - The function extracts the first regex capture from each strain value.
    - Example with default `r"[A-Za-z]+"`:
        - `"strainA1" -> "strainA"`
        - `"strainA2" -> "strainA"`
    - Rules in `custom_rules` are matched against this extracted `strain_group`, not the original raw `strain` string.

Important consequence:

- If extraction fails for a row, that row gets `NaN` in `strain_group` and will not match rules.
- Only strains matched by rules are processed; unmatched strains are skipped.

## 4) How `custom_rules` Are Interpreted

Each entry in `custom_rules` can be either:

- Exact-name rule (no `pattern` key):
    - Rule key is treated as an exact `strain_group` value.

- Pattern rule (has `pattern` key):
    - `pattern` regex is matched against all unique `strain_group` values.

Validation behavior:

- `pattern` must be a string.
- Invalid regex raises `ValueError`.
- If one strain matches multiple rules, `ValueError` is raised.
- `direction` must be `alphabetical/alpha/rows/numerical/num/columns`.
- In fixed-size mode, `sample_size` must be `int >= 1`.

## 5) Core Interaction: `pattern` vs `direction` vs `sample_size`

This is the most important part.

### Scenario A: Exact rule (no `pattern`)

Rule example:

```python
"control": {"direction": "alphabetical", "sample_size": 3}
```

Behavior:

- Matches only `strain_group == "control"`.
- Uses fixed-size grouping by `sample_size`.
- If `sample_size` is omitted, default is `3`.
- If well count is not divisible by `sample_size`, raises `ValueError`.
- Output `strain` column value is the rule key (`"control"`).

### Scenario B: Pattern rule with `sample_size` provided

Rule example:

```python
"EGMB_rule": {"pattern": r"EGMB_\d+", "direction": "alphabetical", "sample_size": 2}
```

Behavior:

- Finds all matching strains (for example `EGMB_61`, `EGMB_62`).
- Each matching strain is processed independently.
- Fixed-size grouping is used inside each strain (`sample_size=2`).
- Divisibility check is applied per strain.
- Output `strain` values are the matched strain names (`EGMB_61`, `EGMB_62`), not the rule key.

### Scenario C: Pattern rule without `sample_size`

Rule example:

```python
"EGMB": {"pattern": r"EGMB_\d+", "direction": "alphabetical"}
```

Behavior:

- All matching strains are combined into one logical output strain (the rule key, here `"EGMB"`).
- Grouping changes from fixed-size to full-direction groups:
    - `alphabetical`: group by full row (`well_row`).
    - `numerical`: group by full column (`well_col`).
- Replicate counts can vary between groups.
- No divisibility check against sample size is needed in this mode.

This mode is useful when plate layout is directional but replicate count differs by row/column.

### Scenario D: Pattern rule with explicit `sample_size=None`

Behavior is the same as Scenario C (full-direction grouping), because the code checks:

- `pattern` exists, and
- `sample_size is None`.

## 6) Complete Combination Matrix

| `strain_pattern` | Rule Type | `sample_size` in rule | Matching Target | Grouping Method | Output `strain` value |
|---|---|---|---|---|---|
| `None` | Exact | provided/missing | Original `strain` | Fixed-size blocks | Rule key |
| `None` | Pattern | provided (or defaulted by explicit int) | Original `strain` | Fixed-size blocks per matched strain | Matched strain name |
| `None` | Pattern | missing or `None` | Original `strain` | Full row/column groups | Rule key |
| Regex string | Exact | provided/missing | Extracted `strain_group` | Fixed-size blocks | Rule key |
| Regex string | Pattern | provided | Extracted `strain_group` | Fixed-size blocks per matched strain | Matched extracted strain |
| Regex string | Pattern | missing or `None` | Extracted `strain_group` | Full row/column groups | Rule key |

## 7) Output Semantics

The function returns one row per group and timepoint with columns such as:

- `strain`
- `group_id`
- `wells`, `well_rows`, `well_cols`
- `{value_column_name}_mean`, `{value_column_name}_std`
- `n_replicates`

Notes:

- `n_replicates` is based on non-null values in the value column.
- `group_id` format depends on direction:
    - Alphabetical: `A_1-3` style
    - Numerical: `ABC_1` style
- Time sorting uses first available in: `time_min`, `time_h`, `time_label`.

## 8) Practical Guidance

- Use `strain_pattern=None` when your raw strain labels are already rule-ready.
- Use a regex `strain_pattern` when you need to normalize names before rule matching.
- Use pattern rule + `sample_size` when each matched strain should remain separate.
- Use pattern rule without `sample_size` when matched strains should be merged and grouped by full rows/columns.
- Make patterns mutually exclusive to avoid conflict errors.

## 9) Minimal Examples

Exact rules only:

```python
custom_rules = {
        "control": {"direction": "alphabetical", "sample_size": 3},
        "blank": {"direction": "alphabetical", "sample_size": 1},
}
stats = calculate_replicate_statistics_by_custom(df, strain_pattern=None, custom_rules=custom_rules)
```

Pattern rules, per-strain fixed grouping:

```python
custom_rules = {
        "experimental": {"pattern": r"st_\d+_r\d+", "direction": "numerical", "sample_size": 4},
}
stats = calculate_replicate_statistics_by_custom(df, strain_pattern=None, custom_rules=custom_rules)
```

Pattern rules, merged variable replicate groups:

```python
custom_rules = {
        "EGMB": {"pattern": r"EGMB_\d+", "direction": "alphabetical"},
}
stats = calculate_replicate_statistics_by_custom(df, strain_pattern=None, custom_rules=custom_rules)
```

## 10) Actual Examples

### Alfie
1. The `strain_pattern=([A-Za-z]+)` is a *regex* that group starins by thier initial charachters only, for example, "strainA1" and "strainA2" are grouped as "strainA".
2. The custum rule does not have pattern with `direction=alphabetical` and `sample_size=3`, each rule's name directly matches with the strain group. As a result, each group **must** have `sample_size=3` replicates.

### Anthony
1. The `strain_pattern=([A-Za-z0-9]+)(?:_\d+)?_(R\d+)` is a *regex*. The middle optional `(?:_\d+)?` is non-capturing group and both `([A-Za-z0-9]+)` and `_(R\d+)` are captured, and are joined to gether. For example, "EGMB_61_R1", "EGMB_61_R1" and "EGMB_61_R1" are grouped as "EGMB_R1", and "MG1655_R12", MG1655_R1" are grouped as  "MG1655_R12", MG1655_R1".
2. The custum rule is `pattern=EGMB`, or `pattern=MG1655`, `direction=numerical` and `sample_size=None`, so, each strain group first simplifies to its strain names, then, they groups along columns (e.g., since `direction=numerical`, one group can be a1, b2, c1, d1.) and replicate size (or group size for columns) can be varied (`sample_size=None`).
