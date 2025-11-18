# Ruff Migration Summary

## Changes Made

Successfully migrated from **Pylint + Black + isort** to **Ruff** (all-in-one tool).

### 1. Updated `.vscode/settings.json`

- **Removed**: Pylint, Black formatter, isort
- **Added**: Ruff as the default formatter and linter
- **Configured**: Format on save and auto-fix on save

### 2. Created `ruff.toml`

Configuration file with:
- Line length: 120 (matching your previous Black config)
- Target: Python 3.13
- Enabled rules: pycodestyle, Pyflakes, isort, pep8-naming, pyupgrade, bugbear, comprehensions, simplify
- Ignored rules that match your previous Pylint configuration

### 3. Installed Ruff

- Installed `ruff` package in your Python environment
- Version: 0.14.5

## Benefits of Ruff

1. **Speed**: 10-100x faster than existing tools
2. **All-in-one**: Replaces Black, isort, Pylint, and more
3. **Auto-fix**: Can automatically fix most issues
4. **Modern**: Suggests Python 3.10+ type hints

## Usage

### In VS Code
- Formatting and linting happen automatically on save
- Issues are underlined in the editor

### Command Line

```powershell
# Check for issues
ruff check .

# Fix auto-fixable issues
ruff check --fix .

# Format code
ruff format .

# Check specific file
ruff check src/labUtils/amn_mappings.py
```

## Current Status

Ruff found 51 suggestions in `amn_mappings.py`:
- 46 auto-fixable (mostly type hint modernization)
- 2 require manual review (simplify if-else blocks)

Most suggestions are about using Python 3.10+ type hints:
- `Dict[str, str]` → `dict[str, str]`
- `List[str]` → `list[str]`
- `Optional[str]` → `str | None`
- `Union[str, Path]` → `str | Path`

## Next Steps (Optional)

1. **Auto-fix suggestions**:
   ```powershell
   ruff check --fix src/
   ```

2. **Format all code**:
   ```powershell
   ruff format src/
   ```

3. **Update dev-requirements.txt**:
   - Remove: `black`, `isort`, `pylint`
   - Add: `ruff`

## Configuration Files

- **Active**: `ruff.toml` (project root)
- **Deprecated**: Can remove `.pylintrc` if it exists
