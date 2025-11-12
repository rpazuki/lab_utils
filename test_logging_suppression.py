"""Test that COBRApy INFO logging is suppressed"""

import tempfile
from pathlib import Path

from labUtils.amn_mappings import (parse_sbml_exchange_bounds,
                                   parse_sbml_exchanges)

print("Testing COBRApy logging suppression...")
print("=" * 70)

model_path = Path(tempfile.gettempdir()) / "iML1515.xml"

if model_path.exists():
    print("\n1. Testing parse_sbml_exchanges (should NOT see 'cobra.core.model - INFO'):")
    mapping = parse_sbml_exchanges(model_path)
    print(f"   ✓ Parsed {len(mapping)} metabolite mappings")

    print("\n2. Testing parse_sbml_exchange_bounds (should NOT see 'cobra.core.model - INFO'):")
    bounds = parse_sbml_exchange_bounds(model_path)
    print(f"   ✓ Parsed {len(bounds)} exchange bounds")

    print("\n" + "=" * 70)
    print("✓ If you don't see 'cobra.core.model - INFO' messages above,")
    print("  the logging suppression is working correctly!")
    print("=" * 70)
else:
    print(f"\niML1515.xml not found at {model_path}")
    print("Please run the test_sbml_parsing.py script first to download it.")
