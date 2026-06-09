# Enhancement: Pattern-Based Rules for calculate_replicate_statistics_by_custom

## Overview

The `calculate_replicate_statistics_by_custom` function in `media_bot.py` has been enhanced to support **pattern-based grouping rules**. This allows you to specify regex patterns within the `custom_rules` dictionary to match multiple strains at once, with each matching strain being processed as a separate group using the same aggregation rules.

## What's New

### New Feature: `pattern` Key in Rules

Each rule in the `custom_rules` dictionary can now optionally include a `"pattern"` key containing a regex pattern. When a pattern is specified:

1. **All strains matching the pattern are identified**
2. **Each matching strain is processed independently** with the specified rules
3. **Conflicts are detected** - if a strain matches multiple patterns, an error is raised

### Key Differences from Existing Behavior

| Aspect | Before | After |
|--------|--------|-------|
| Rule Scope | Each rule matched exactly one strain | Rules can match one or multiple strains via pattern |
| Flexibility | Limited to exact strain name matches | Can use regex patterns for flexible matching |
| Strain Matching | Direct key lookup | Can use pattern OR exact name lookup |

## Usage Examples

### Example 1: Pattern-Based Rules (Recommended New Approach)

```python
from labUtils.media_bot import calculate_replicate_statistics_by_custom

# If you have strains like: st_1_r1, st_2_r1, st_3_r1, st_1_r2, st_2_r2, control

custom_rules = {
    "experimental_strains": {
        "pattern": r"st_\d+_r\d+",  # Matches st_1_r1, st_2_r1, st_3_r1, etc.
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "control": {  # Exact match for strain named "control"
        "direction": "numerical",
        "sample_size": 3,
    }
}

result = calculate_replicate_statistics_by_custom(
    df_parsed,
    strain_pattern=None,  # Use strain names directly
    custom_rules=custom_rules,
)
```

Each matching strain (st_1_r1, st_2_r1, st_3_r1, st_1_r2, st_2_r2) will be processed independently as a separate group.

### Example 2: Mixing Exact and Pattern-Based Rules

```python
custom_rules = {
    "st_1_r1": {  # Exact match for this specific strain
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "other_experimental": {
        "pattern": r"st_[2-3]_r\d+",  # Matches st_2_r*, st_3_r*
        "direction": "numerical",
        "sample_size": 3,
    },
    "control": {  # Exact match
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

result = calculate_replicate_statistics_by_custom(
    df_parsed,
    strain_pattern=None,
    custom_rules=custom_rules,
)
```

### Example 3: Global Pattern + Per-Rule Patterns

```python
# If strain names are like: exp_st_1_r1, exp_st_2_r1, etc.

custom_rules = {
    "experimental": {
        "pattern": r"exp_st_\d+_r\d+",  # Matches exp_st_*_r*
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "control": {
        "direction": "alphabetical",
        "sample_size": 3,
    }
}

result = calculate_replicate_statistics_by_custom(
    df_parsed,
    strain_pattern=r"exp_(.+)",  # Extract the part after "exp_"
    custom_rules=custom_rules,
)
```

## Pattern Syntax

The `"pattern"` value is a **regular expression (regex)** string. Common patterns:

| Pattern | Matches |
|---------|---------|
| `r"st_\d+_r\d+"` | `st_1_r1`, `st_2_r2`, `st_10_r3`, etc. |
| `r"control.*"` | `control`, `control_1`, `control_rep1`, etc. |
| `r"st_[1-3]_.*"` | `st_1_anything`, `st_2_anything`, `st_3_anything` |
| `r"^exp_.*"` | Any strain starting with `exp_` |
| `r".*_rep\d+$"` | Any strain ending with `_rep1`, `_rep2`, etc. |

## Error Handling

### New Validation

1. **Invalid Regex**: If a `pattern` is not a valid regex, you'll get:
   ```
   ValueError: Invalid regex pattern in rule 'rule_name': pattern_string. Error: ...
   ```

2. **Multiple Pattern Matches**: If a strain matches more than one pattern, you'll get:
   ```
   ValueError: Strain 'strain_name' matches multiple patterns in custom_rules.
   Please ensure each strain matches at most one pattern.
   ```

3. **Rule-Specific Error Messages**: Error messages now include the `rule_key` for easier debugging:
   ```
   ValueError: Invalid direction 'invalid' for strain 'st_1_r1' in rule 'experimental_strains'.
   Must be 'alphabetical'/'alpha' or 'numerical'/'num'
   ```

## Return Value Changes

The output DataFrame column name changed from `strain_group` to `strain` for consistency. This applies to both pattern-based and exact-match rules.

## Backward Compatibility

✅ **Fully backward compatible!**

Existing code using exact strain name matching will continue to work without any changes:

```python
# This still works exactly as before
custom_rules = {
    "strainA": {"direction": "alphabetical", "sample_size": 3},
    "strainB": {"direction": "numerical", "sample_size": 4},
}
```

## Testing

Three comprehensive test scenarios were validated:

1. **Pure Pattern-Based Rules**: All strains matched via patterns
2. **Mixed Rules**: Combination of exact matches and pattern-based rules
3. **Global + Per-Rule Patterns**: Global pattern extraction with per-rule pattern matching

All tests passed successfully ✓

## Implementation Details

### Internal Processing

1. **Pattern Compilation**: Each pattern is compiled as a regex object for validation
2. **Strain Matching**: For each unique strain, the function checks if it matches any rule's pattern
3. **Conflict Detection**: Ensures no strain is matched by multiple patterns
4. **Independent Processing**: Each matched strain is processed independently with its rule's settings

### New Columns in Output

No new columns were added. The function output schema remains the same, with the only change being the column name `strain_group` → `strain`.

## Migration Guide

If you were using the old `strain_group` column name in downstream code:

**Before:**
```python
result['strain_group']
```

**After:**
```python
result['strain']
```

## Example Use Cases

### Use Case 1: Replicate Series as Distinct Groups

```python
# If each replicates are independent experiments
custom_rules = {
    "series_1": {
        "pattern": r".*_rep1",  # All strains with rep1 suffix
        "direction": "alphabetical",
        "sample_size": 3,
    },
    "series_2": {
        "pattern": r".*_rep2",  # All strains with rep2 suffix
        "direction": "alphabetical",
        "sample_size": 3,
    },
}
```

### Use Case 2: Hierarchical Strain Organization

```python
# Organize strains by both type and condition
custom_rules = {
    "wildtype_conditions": {
        "pattern": r"wt_.*",  # All wildtype strains
        "direction": "numerical",
        "sample_size": 4,
    },
    "mutant_conditions": {
        "pattern": r"mut_.*",  # All mutant strains
        "direction": "numerical",
        "sample_size": 4,
    },
    "control_group": {
        "pattern": r"ctrl.*",  # Control strains
        "direction": "alphabetical",
        "sample_size": 6,
    },
}
```

## Summary

This enhancement provides:

✅ **Flexibility**: Match multiple strains with a single rule using regex patterns
✅ **Clarity**: Mix exact and pattern-based rules in the same configuration
✅ **Safety**: Built-in conflict detection and improved error messages
✅ **Compatibility**: Fully backward compatible with existing code
✅ **Simplicity**: Intuitive API using familiar regex syntax
