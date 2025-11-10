# calculate_replicate_statistics() Output Format

## Function Returns

The function returns a DataFrame with **only averaged values**, not the original data.

## Output Columns

| Column Name | Description | Example |
|-------------|-------------|---------|
| `group_id` | Identifier for the replicate group | "A_1-3" (alphabetical) or "ABC_1" (numerical) |
| `wells` | Comma-separated list of wells in group | "A1,A2,A3" |
| `well_rows` | Comma-separated list of well rows | "A,A,A" |
| `well_cols` | Comma-separated list of well columns | "1,2,3" |
| `od_mean` | Mean OD value across replicates | 0.089 |
| `od_std` | Standard deviation of OD | 0.0082 |
| `n_replicates` | Number of non-null values used | 3 |
| Time columns | Preserved from input (time_h, time_min, etc.) | 0.0, 15.0, 30.0, ... |
| Metadata columns | First well's metadata (strain, media_type, etc.) | Preserved if present |

## Example Usage

```python
from pathlib import Path
from labUtils.media_bot import parse, calculate_replicate_statistics

# Parse raw data
df = parse(raw_data_path, meta_data_path)

# Calculate statistics for 3 wells along rows
stats = calculate_replicate_statistics(
    df,
    direction="alphabetical",  # or "alpha"
    sample_size=3,
    ddof=1  # unbiased std (divide by N-1)
)

# Output columns
print(stats.columns)
# ['group_id', 'wells', 'well_rows', 'well_cols', 'od_mean', 'od_std',
#  'n_replicates', 'time_h', 'time_min', 'strain', 'media_type', ...]
```

## Sample Output

### Alphabetical Direction (sample_size=3)
```
   group_id      wells  well_rows well_cols  od_mean   od_std  n_replicates  time_h
0     A_1-3  A1,A2,A3      A,A,A     1,2,3   0.0890  0.00819             3     0.0
1     A_1-3  A1,A2,A3      A,A,A     1,2,3   0.0993  0.00252             3     0.25
2     B_1-3  B1,B2,B3      B,B,B     1,2,3   0.0920  0.00436             3     0.0
```

### Numerical Direction (sample_size=3)
```
   group_id      wells  well_rows well_cols  od_mean   od_std  n_replicates  time_h
0     ABC_1  A1,B1,C1    A,B,C       1,1,1   0.0837  0.00550             3     0.0
1     ABC_1  A1,B1,C1    A,B,C       1,1,1   0.0963  0.00400             3     0.25
2     ABC_2  A2,B2,C2    A,B,C       2,2,2   0.0937  0.00404             3     0.0
```

## Key Features

1. **Only aggregated data** - Original OD values are not included
2. **Column names** - `od_mean` and `od_std` (not od600_mean/od600_std)
3. **Well information** - Three columns (`wells`, `well_rows`, `well_cols`) provide complete grouping details
4. **Flexible std calculation** - Use `ddof=1` for sample std or `ddof=0` for population std
