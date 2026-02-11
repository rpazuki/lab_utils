import json
import requests
from .fba_tools import load_model

def normalize_stoichiometry(reaction):
    """Return a frozenset of metabolite:coefficient pairs for comparison."""
    return frozenset((met.id, coeff) for met, coeff in reaction.metabolites.items())

def build_reaction_mapping_ijo_to_iml(model_ijo1366, model_iml1515):
    """
    Map reactions from iJO1366 to iML1515 based on:
    1. Exact ID match
    2. Stoichiometry match
    Returns: dict {ijo1366_id: iml1515_id}
    """
    mapping = {}
    iml_stoich_map = {normalize_stoichiometry(r): r.id for r in model_iml1515.reactions}

    for r_ijo in model_ijo1366.reactions:
        # Case 1: Exact ID match
        if r_ijo.id in model_iml1515.reactions:
            mapping[r_ijo.id] = r_ijo.id
            continue

        # Case 2: Stoichiometry match
        stoich = normalize_stoichiometry(r_ijo)
        if stoich in iml_stoich_map:
            mapping[r_ijo.id] = iml_stoich_map[stoich]

    return mapping

def create_iml1515_map(source_map_name, iml1515_model_path, ijo1366_model_path, output_path):
    """
    Create a new Escher map for iML1515 by adapting an iJO1366 map.

    Args:
        source_map_name: Name of the iJO1366 map (e.g., 'iJO1366.Nucleotide metabolism')
        iml1515_model_path: Path to iML1515 model
        ijo1366_model_path: Path to iJO1366 model
        output_path: Where to save the new iML1515 map
    """
    print(f"Loading models...")
    model_iml = load_model(iml1515_model_path)
    model_ijo = load_model(ijo1366_model_path)
    print(f"  iML1515: {len(model_iml.reactions)} reactions")
    print(f"  iJO1366: {len(model_ijo.reactions)} reactions")

    print(f"\nDownloading {source_map_name} map...")
    # Download the map directly from Escher's GitHub
    map_url = f"https://escher.github.io/1-0-0/6/maps/Escherichia%20coli/{source_map_name.replace(' ', '%20')}.json"
    print(f"  URL: {map_url}")

    response = requests.get(map_url)
    if response.status_code != 200:
        raise ValueError(f"Failed to download map. Status code: {response.status_code}")

    map_data = response.json()
    print(f"  Map downloaded successfully")
    print(f"  Map type: {type(map_data)}")

    print(f"\nBuilding reaction mapping from iJO1366 to iML1515...")
    reaction_mapping = build_reaction_mapping_ijo_to_iml(model_ijo, model_iml)
    print(f"  Found {len(reaction_mapping)} mapped reactions")

    print(f"\nTranslating map reaction IDs...")

    translated_count = 0
    not_found_count = 0

    # Find the reactions dictionary in the map structure
    map_object = None
    if isinstance(map_data, list):
        print(f"  Map is a list with {len(map_data)} items")
        # Try to find the map object in the list
        for i, item in enumerate(map_data):
            if isinstance(item, dict):
                print(f"    Item {i}: dict with keys {list(item.keys())[:5]}")
                if 'reactions' in item:
                    map_object = item
                    print(f"    Found reactions in item {i}")
                    break
    elif isinstance(map_data, dict):
        print(f"  Map is a dict with keys: {list(map_data.keys())[:10]}")
        map_object = map_data

    if map_object and 'reactions' in map_object:
        reactions_dict = map_object['reactions']
        print(f"  Found {len(reactions_dict)} reactions in map")

        for rxn_key, rxn_data in reactions_dict.items():
            if isinstance(rxn_data, dict) and 'bigg_id' in rxn_data:
                old_id = rxn_data['bigg_id']
                if old_id in reaction_mapping:
                    new_id = reaction_mapping[old_id]
                    rxn_data['bigg_id'] = new_id
                    translated_count += 1
                else:
                    not_found_count += 1

        print(f"  Translated: {translated_count} reactions")
        print(f"  Not found in iML1515: {not_found_count} reactions")
    else:
        print("  ERROR: Could not find reactions in map structure!")
        raise ValueError("Map structure doesn't contain reactions")

    print(f"\nSaving new map to {output_path}...")
    with open(output_path, 'w') as f:
        json.dump(map_data, f, indent=2)

    print(f"\n✓ Success! Created iML1515 map at: {output_path}")
    return output_path
