import json
import requests
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
from escher import Builder
from .fba_tools import load_model



def plot_fba_predictions(df_data,
                         gr_column,
                         per_strain: bool = True,
                         replication: str = "replicates",
                         strain: str = "purB",
                         experiment: str = "mediabotJLF1",
                         title=""):
    fig = plt.figure(figsize=(10,6))
    ax1 = fig.subplots(1,1)

    # Get filtered data for plotting
    y = np.array([p for p, status in df_data[["fba_growth", "fba_status"]].itertuples(index=False)
                if not np.isnan(p) and status == 'optimal'])
    x = np.array([v for p, status, v in df_data[["fba_growth", "fba_status", gr_column]].itertuples(index=False)
                if not np.isnan(p) and status == 'optimal'])

    # Get supplement groups for color coding
    groups = np.array([g for p, status, g in df_data[["fba_growth", "fba_status", "group"]].itertuples(index=False)
                        if not np.isnan(p) and status == 'optimal']).astype(str)
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)

    # Use a colormap with distinct colors
    if n_groups <= 10:
        colors = plt.cm.tab10(np.linspace(0, 1, 10))[:n_groups] # type: ignore
    elif n_groups <= 20:
        colors = plt.cm.tab20(np.linspace(0, 1, 20))[:n_groups] # type: ignore
    else:
        colors = plt.cm.gist_rainbow(np.linspace(0, 1, n_groups)) # type: ignore

    # Create color map for supplements
    group_colors = {str(supp): colors[i] for i, supp in enumerate(unique_groups)}

    # Plot each supplement group with its own color
    for supp in unique_groups:
        supp_mask = groups == supp

        ax1.scatter(x[supp_mask], y[supp_mask],
                    color=group_colors[supp],
                    label=supp,
                    alpha=0.7,
                    s=50,
                    edgecolors='black',
                    linewidths=0.5)
    ax1.plot([min(x), max(x)], [min(x), max(x)], color='gray', linestyle='--', label='y=x')

    # Add legend
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8,
                title='Supplement', framealpha=0.9)
    m, b = np.polyfit(x, y, 1)
    ax1.plot(x, m*x + b, color='red', label=f'y={m:.2f}x+{b:.2f}')

    # Calculate both R² of y = x, correlation, and R2_score metrics
    y_pred = x
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    correlation_matrix = np.corrcoef(x, y)
    correlation_xy = correlation_matrix[0,1]
    r_squared_corr = correlation_xy**2
    r2_uncalibrated = r2_score(x, y, multioutput='variance_weighted')

    # Calculate R² with calibrated predictions
    # y_calibrated = m * x + b
    # r2_calibrated = r2_score(y, y_calibrated)

    # print(f"Calibration parameters: m={m:.4f}, b={b:.4f}")
    # print(f"Correlation²: {r_squared_corr:.4f}")
    # print(f"R² (uncalibrated, y vs x): {r2_uncalibrated:.4f}")
    # print(f"R² (calibrated, y vs mx+b): {r2_calibrated:.4f}")
    # print(f"\nNote: R²(calibrated) ≈ Correlation² as expected!")

    ax1.set_xlabel("Experimental growth rates")
    ax1.set_ylabel("FBA predicted growth rates")
    if per_strain:
        subtitle= (f"Flux converted - Strain:{strain}, "
                    f"\n Data points n:{len(x)}")
    else:
        subtitle= (f"Flux converted - experiment:{experiment},  "
                    f"\n Data points n:{len(x)}")
    ax1.set_title(f"{title} {subtitle}"
                f"(Corr²={r_squared_corr:.4f}, R²={r2:.4f}, R² scored={r2_uncalibrated:.4f}) \n"
                f"Calibration: y={m:.3f}x+{b:.3f} "
                )
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    groups = df_data['group'].astype(str).values
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    n_groups = len(unique_groups)
    n_cols = min(3, n_groups)
    n_rows = int(np.ceil(n_groups / n_cols))

    fig2, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_groups == 1:
        axes = np.array([axes])
    axes = axes.flatten() if n_groups > 1 else axes

    for idx, supp in enumerate(unique_groups):
        ax = axes[idx] if n_groups > 1 else axes[0]

        # Get data for this supplement group
        supp_mask = groups == supp
        x_supp = x[supp_mask]
        y_supp = y[supp_mask]

        # Get the individual supplements within this group for coloring
        supp_supplements = df_data['supplements_unified'].values[supp_mask]
        unique_supps_in_group = np.unique(supp_supplements)
        n_supps_in_group = len(unique_supps_in_group)

        # Create color map for individual supplements within this group
        if n_supps_in_group <= 10:
            supp_colors_local = plt.cm.tab10(np.linspace(0, 1, 10))[:n_supps_in_group] # type: ignore
        elif n_supps_in_group <= 20:
            supp_colors_local = plt.cm.tab20(np.linspace(0, 1, 20))[:n_supps_in_group] # type: ignore
        else:
            supp_colors_local = plt.cm.gist_rainbow(np.linspace(0, 1, n_supps_in_group)) # type: ignore

        supp_color_map = {s: supp_colors_local[i] for i, s in enumerate(unique_supps_in_group)}

        # Plot scatter with individual supplement colors
        for individual_supp in unique_supps_in_group:
            individual_mask = supp_supplements == individual_supp
            ax.scatter(x_supp[individual_mask], y_supp[individual_mask],
                      color=supp_color_map[individual_supp],
                      alpha=0.7,
                      s=50,
                      edgecolors='black',
                      linewidths=0.5,
                      label=individual_supp)

        # Fit line and calculate R²
        if len(x_supp) > 1:
            m_supp, b_supp = np.polyfit(x_supp, y_supp, 1)
            ax.plot(x_supp, m_supp*x_supp + b_supp,
                   color='red',
                   linewidth=2,
                   label=f'y={m_supp:.2f}x+{b_supp:.2f}')

            # Calculate R²
            correlation_matrix_supp = np.corrcoef(x_supp, y_supp)
            correlation_xy_supp = correlation_matrix_supp[0,1]
            r_squared_supp = correlation_xy_supp**2

            ax.set_title(f'{supp}\nR²={r_squared_supp:.4f}, n={len(x_supp)}',
                        fontsize=10, fontweight='bold')
        else:
            ax.set_title(f'{supp}\nn={len(x_supp)} (insufficient data)',
                        fontsize=10)

        ax.set_xlabel("Experimental growth rates")
        ax.set_ylabel("FBA predicted growth rates")
        ax.grid(True, alpha=0.3)
        ax.plot([min(x_supp), max(x_supp)], [min(x_supp), max(x_supp)], color='gray', linestyle='--', label='y=x')
        # if len(x_supp) > 1:
        #     ax.legend(fontsize=6, ncol=2, loc='best')
    average_type = "All data points" if replication == "no_replicates" else "OD average" if replication == "replicates" else "Growth rate average"
    # Hide unused subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].axis('off')

    fig2.suptitle(f'Individual Supplement Group Analysis - {strain if per_strain else experiment}, for {average_type}',
                 fontsize=14, fontweight='bold')
    fig2.tight_layout()
    plt.show()

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

def _build_escher_map_and_builder_from_reactions(
    reaction_objects,
    map_name,
    map_id,
    map_description,
    map_json_path,
    model,
    reaction_data,
    metabolite_data,
    hide_secondary_metabolites=True,
    hide_all_labels=False,
    reaction_scale_preset="GaBuRd",
):
    import json
    import os

    output_dir = os.path.dirname(map_json_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    nodes = {}
    reactions = {}
    metabolite_node_ids = {}
    node_id_counter = 1
    reaction_id_counter = 1
    segment_id_counter = 1

    def new_node_id():
        nonlocal node_id_counter
        node_id = str(node_id_counter)
        node_id_counter += 1
        return node_id

    def get_or_create_metabolite_node(metabolite, x, y):
        if metabolite.id in metabolite_node_ids:
            return metabolite_node_ids[metabolite.id]
        node_id = new_node_id()
        nodes[node_id] = {
            "node_type": "metabolite",
            "x": x,
            "y": y,
            "bigg_id": metabolite.id,
            "name": metabolite.name,
            "label_x": x + 20,
            "label_y": y + 20,
            "node_is_primary": False,
        }
        metabolite_node_ids[metabolite.id] = node_id
        return node_id

    for row_index, rxn in enumerate(reaction_objects):
        y_base = 200 + (row_index * 220)

        left_marker_id = new_node_id()
        mid_marker_id = new_node_id()
        right_marker_id = new_node_id()

        nodes[left_marker_id] = {"node_type": "multimarker", "x": 320, "y": y_base}
        nodes[mid_marker_id] = {"node_type": "midmarker", "x": 420, "y": y_base}
        nodes[right_marker_id] = {"node_type": "multimarker", "x": 520, "y": y_base}

        segments = {}

        segments[str(segment_id_counter)] = {
            "from_node_id": left_marker_id,
            "to_node_id": mid_marker_id,
            "b1": None,
            "b2": None,
        }
        segment_id_counter += 1

        segments[str(segment_id_counter)] = {
            "from_node_id": mid_marker_id,
            "to_node_id": right_marker_id,
            "b1": None,
            "b2": None,
        }
        segment_id_counter += 1

        reactants = [(met, coeff) for met, coeff in rxn.metabolites.items() if coeff < 0]
        products = [(met, coeff) for met, coeff in rxn.metabolites.items() if coeff > 0]

        for met_index, (metabolite, _) in enumerate(reactants):
            met_y = y_base - 60 + (met_index * 40)
            met_node_id = get_or_create_metabolite_node(metabolite, x=180, y=met_y)
            segments[str(segment_id_counter)] = {
                "from_node_id": met_node_id,
                "to_node_id": left_marker_id,
                "b1": None,
                "b2": None,
            }
            segment_id_counter += 1

        for met_index, (metabolite, _) in enumerate(products):
            met_y = y_base - 60 + (met_index * 40)
            met_node_id = get_or_create_metabolite_node(metabolite, x=660, y=met_y)
            segments[str(segment_id_counter)] = {
                "from_node_id": right_marker_id,
                "to_node_id": met_node_id,
                "b1": None,
                "b2": None,
            }
            segment_id_counter += 1

        reactions[str(reaction_id_counter)] = {
            "name": rxn.name,
            "bigg_id": rxn.id,
            "reversibility": bool(rxn.reversibility),
            "label_x": 420,
            "label_y": y_base - 30,
            "gene_reaction_rule": rxn.gene_reaction_rule or "",
            "genes": [{"bigg_id": gene.id, "name": gene.name} for gene in rxn.genes],
            "metabolites": [
                {"bigg_id": met.id, "coefficient": coeff}
                for met, coeff in rxn.metabolites.items()
            ],
            "segments": segments,
        }
        reaction_id_counter += 1

    for rxn in reaction_objects:
        for met, coeff in rxn.metabolites.items():
            if met.id in metabolite_node_ids and coeff != 0:
                node_id = metabolite_node_ids[met.id]
                if not nodes[node_id].get("node_is_primary", False):
                    nodes[node_id]["node_is_primary"] = True
                    break

    canvas_height = max(600, 260 + (len(reaction_objects) * 220))

    map_data = [
        {
            "map_name": map_name,
            "map_id": map_id,
            "map_description": map_description,
            "homepage": "https://escher.github.io",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        },
        {
            "reactions": reactions,
            "nodes": nodes,
            "text_labels": {},
            "canvas": {"x": 0, "y": 0, "width": 900, "height": canvas_height},
        },
    ]

    with open(map_json_path, "w") as f:
        json.dump(map_data, f, indent=2)

    builder = Builder(
        map_json=map_json_path,
        model=model,
        reaction_data=reaction_data,
        metabolite_data=metabolite_data,
    )

    builder.hide_secondary_metabolites = hide_secondary_metabolites
    builder.hide_all_labels = hide_all_labels
    builder.reaction_scale_preset = reaction_scale_preset

    return builder, {
        "reaction_count": len(reaction_objects),
        "metabolite_count": len(metabolite_node_ids),
        "map_json_path": map_json_path,
    }

def build_reaction_escher_builder(
    reaction_name,
    df_fluxes,
    df_shadow_prices,
    model,
    solution_index=0,
    search_depth=0,
    expansion_mode="both",
    nonzero_only=False,
    flux_threshold=1e-9,
    keep_seed_reaction=True,
    map_json_path=None,
    hide_secondary_metabolites=True,
    hide_all_labels=False,
    reaction_scale_preset="GaBuRd",
    verbose=True,
):
    """Build and return an Escher Builder for a reaction neighborhood.

    Neighborhood definition:
    - depth 0: only the seed reaction
    - depth 1: reactions sharing selected metabolites with the seed
    - depth 2+: repeat the same expansion from newly found reactions

    expansion_mode controls which metabolites are used for expansion:
    - both: use reactants and products
    - reactants_only: use only metabolites with negative coefficients
    - products_only: use only metabolites with positive coefficients

    nonzero_only controls plotting filter based on flux magnitude:
    - False: include all discovered reactions
    - True: include only reactions with |flux| > flux_threshold
    """

    if solution_index >= len(df_fluxes):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for fluxes length {len(df_fluxes)}"
        )
    if solution_index >= len(df_shadow_prices):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for shadow prices length {len(df_shadow_prices)}"
        )
    if not isinstance(search_depth, int) or search_depth < 0:
        raise ValueError("search_depth must be a non-negative integer")
    if flux_threshold < 0:
        raise ValueError("flux_threshold must be non-negative")

    valid_expansion_modes = {"both", "reactants_only", "products_only"}
    if expansion_mode not in valid_expansion_modes:
        raise ValueError(
            f"expansion_mode must be one of {sorted(valid_expansion_modes)}; got '{expansion_mode}'"
        )

    seed_reaction = model.reactions.get_by_id(reaction_name)
    flux_row = df_fluxes.iloc[solution_index]

    def metabolites_for_expansion(reaction_obj):
        if expansion_mode == "both":
            return [met for met, coeff in reaction_obj.metabolites.items() if coeff != 0]
        if expansion_mode == "reactants_only":
            return [met for met, coeff in reaction_obj.metabolites.items() if coeff < 0]
        return [met for met, coeff in reaction_obj.metabolites.items() if coeff > 0]

    visited_reaction_ids = {seed_reaction.id}
    reaction_depth = {seed_reaction.id: 0}
    frontier = {seed_reaction.id}

    for depth in range(1, search_depth + 1):
        next_frontier = set()
        for rxn_id in frontier:
            rxn_obj = model.reactions.get_by_id(rxn_id)
            for metabolite in metabolites_for_expansion(rxn_obj):
                for neighbor_reaction in metabolite.reactions:
                    if neighbor_reaction.id not in visited_reaction_ids:
                        visited_reaction_ids.add(neighbor_reaction.id)
                        reaction_depth[neighbor_reaction.id] = depth
                        next_frontier.add(neighbor_reaction.id)
        frontier = next_frontier
        if not frontier:
            break

    def get_flux_value(reaction_id):
        val = flux_row.get(reaction_id, 0.0)
        try:
            return float(val)
        except (TypeError, ValueError):
            return 0.0

    filtered_reaction_ids = set(visited_reaction_ids)
    if nonzero_only:
        filtered_reaction_ids = {
            rid
            for rid in visited_reaction_ids
            if abs(get_flux_value(rid)) > flux_threshold
        }
        if keep_seed_reaction and seed_reaction.id in visited_reaction_ids:
            filtered_reaction_ids.add(seed_reaction.id)

    reaction_objects = [model.reactions.get_by_id(rxn_id) for rxn_id in sorted(filtered_reaction_ids)]

    if not reaction_objects:
        raise ValueError(
            "No reactions left to plot after applying nonzero flux filter. "
            "Lower flux_threshold, disable nonzero_only, or change solution_index."
        )

    if map_json_path is None:
        map_json_path = f"iML1515.{reaction_name}.depth{search_depth}.{expansion_mode}.json"

    map_name = f"iML1515.{reaction_name}.depth{search_depth}.{expansion_mode}"
    map_id = f"iML1515_{reaction_name}_d{search_depth}_{expansion_mode}"
    map_description = (
        f"Auto-generated Escher map for {reaction_name} with "
        f"search_depth={search_depth} and expansion_mode={expansion_mode}"
    )

    builder, map_stats = _build_escher_map_and_builder_from_reactions(
        reaction_objects=reaction_objects,
        map_name=map_name,
        map_id=map_id,
        map_description=map_description,
        map_json_path=map_json_path,
        model=model,
        reaction_data=df_fluxes.iloc[solution_index],
        metabolite_data=df_shadow_prices.iloc[solution_index],
        hide_secondary_metabolites=hide_secondary_metabolites,
        hide_all_labels=hide_all_labels,
        reaction_scale_preset=reaction_scale_preset,
    )

    if verbose:
        depth_counts = {}
        for rxn_id, dep in reaction_depth.items():
            depth_counts[dep] = depth_counts.get(dep, 0) + 1

        print("=" * 60)
        print(f"Seed reaction: {seed_reaction.id}")
        print(f"Requested search_depth: {search_depth}")
        print(f"Expansion mode: {expansion_mode}")
        print(f"Nonzero only: {nonzero_only} (threshold={flux_threshold})")
        print(f"Discovered reactions (topology): {len(visited_reaction_ids)}")
        print(f"Plotted reactions (after filter): {map_stats['reaction_count']}")
        print(f"Discovered metabolites (unique): {map_stats['metabolite_count']}")
        print(f"Map output: {map_stats['map_json_path']}")
        print("Reactions per depth:")
        for dep in sorted(depth_counts):
            print(f"  depth {dep}: {depth_counts[dep]}")
        print("=" * 60)

    return builder

def build_path_between_reactions_escher_builder(
    start_reaction_name,
    end_reaction_name,
    df_fluxes,
    df_shadow_prices,
    model,
    solution_index=0,
    nonzero_only=False,
    flux_threshold=1e-9,
    keep_endpoints=True,
    map_json_path=None,
    hide_secondary_metabolites=True,
    hide_all_labels=False,
    reaction_scale_preset="GaBuRd",
    verbose=True,
):
    from collections import deque

    if solution_index >= len(df_fluxes):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for fluxes length {len(df_fluxes)}"
        )
    if solution_index >= len(df_shadow_prices):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for shadow prices length {len(df_shadow_prices)}"
        )
    if flux_threshold < 0:
        raise ValueError("flux_threshold must be non-negative")

    start_reaction = model.reactions.get_by_id(start_reaction_name)
    end_reaction = model.reactions.get_by_id(end_reaction_name)
    flux_row = df_fluxes.iloc[solution_index]

    def get_flux_value(reaction_id):
        val = flux_row.get(reaction_id, 0.0)
        try:
            return float(val)
        except (TypeError, ValueError):
            return 0.0

    def reaction_is_allowed(reaction_id):
        if not nonzero_only:
            return True
        if keep_endpoints and reaction_id in {start_reaction.id, end_reaction.id}:
            return True
        return abs(get_flux_value(reaction_id)) > flux_threshold

    queue = deque([start_reaction.id])
    visited = {start_reaction.id}
    parent = {start_reaction.id: None}

    found = False
    while queue:
        current_id = queue.popleft()
        if current_id == end_reaction.id:
            found = True
            break

        current_rxn = model.reactions.get_by_id(current_id)
        neighbors = set()
        for met in current_rxn.metabolites:
            for neighbor in met.reactions:
                if neighbor.id != current_id:
                    neighbors.add(neighbor.id)

        for neighbor_id in sorted(neighbors):
            if neighbor_id in visited:
                continue
            if not reaction_is_allowed(neighbor_id):
                continue
            visited.add(neighbor_id)
            parent[neighbor_id] = current_id
            queue.append(neighbor_id)

    if not found:
        raise ValueError(
            f"No path found between {start_reaction.id} and {end_reaction.id} "
            f"under current filter settings (nonzero_only={nonzero_only}, flux_threshold={flux_threshold})."
        )

    reaction_path_ids = []
    node = end_reaction.id
    while node is not None:
        reaction_path_ids.append(node)
        node = parent[node]
    reaction_path_ids.reverse()

    reaction_objects = [model.reactions.get_by_id(rid) for rid in reaction_path_ids]

    if map_json_path is None:
        map_json_path = f"iML1515.path.{start_reaction.id}_to_{end_reaction.id}.json"

    map_name = f"iML1515.path.{start_reaction.id}_to_{end_reaction.id}"
    map_id = f"iML1515_path_{start_reaction.id}_to_{end_reaction.id}"
    map_description = f"Auto-generated Escher path map from {start_reaction.id} to {end_reaction.id}"

    builder, map_stats = _build_escher_map_and_builder_from_reactions(
        reaction_objects=reaction_objects,
        map_name=map_name,
        map_id=map_id,
        map_description=map_description,
        map_json_path=map_json_path,
        model=model,
        reaction_data=df_fluxes.iloc[solution_index],
        metabolite_data=df_shadow_prices.iloc[solution_index],
        hide_secondary_metabolites=hide_secondary_metabolites,
        hide_all_labels=hide_all_labels,
        reaction_scale_preset=reaction_scale_preset,
    )

    if verbose:
        print("=" * 60)
        print(f"Path search: {start_reaction.id} -> {end_reaction.id}")
        print(f"Path length (reactions): {len(reaction_path_ids)}")
        print("Reaction path:")
        print("  " + " -> ".join(reaction_path_ids))
        print(f"Nonzero only: {nonzero_only} (threshold={flux_threshold})")
        print(f"Discovered metabolites (map unique): {map_stats['metabolite_count']}")
        print(f"Map output: {map_stats['map_json_path']}")
        print("=" * 60)

    return builder

def _resolve_metabolite_or_exchange(entity_id, model):
    """Resolve an input ID to a cobra Metabolite.

    Accepts:
    - metabolite ID (e.g., glc__D_e)
    - exchange/sink/demand reaction ID with a single metabolite (e.g., EX_glc__D_e)
    """
    if entity_id in model.metabolites:
        return model.metabolites.get_by_id(entity_id)

    if entity_id in model.reactions:
        rxn = model.reactions.get_by_id(entity_id)
        mets = list(rxn.metabolites.keys())
        if len(mets) == 1:
            return mets[0]
        raise ValueError(
            f"Reaction '{entity_id}' does not map to a single metabolite; "
            "please pass a metabolite ID or a single-metabolite exchange/sink/demand reaction ID."
        )

    raise KeyError(f"'{entity_id}' not found as metabolite or reaction in model")


def _normalize_excluded_metabolites(exclude_metabolites, model):
    excluded_met_ids = set()
    if not exclude_metabolites:
        return excluded_met_ids

    for item in exclude_metabolites:
        met = _resolve_metabolite_or_exchange(item, model)
        excluded_met_ids.add(met.id)

    return excluded_met_ids


def build_metabolite_escher_builder(
    metabolite_name,
    df_fluxes,
    df_shadow_prices,
    model,
    solution_index=0,
    search_depth=0,
    nonzero_only=False,
    flux_threshold=1e-9,
    keep_seed_metabolite=True,
    exclude_metabolites=None,
    map_json_path=None,
    hide_secondary_metabolites=True,
    hide_all_labels=False,
    reaction_scale_preset="GaBuRd",
    verbose=True,
):
    """Build and return an Escher Builder around a metabolite neighborhood.

    Depth logic (metabolite-reaction bipartite expansion):
    - depth 0: reactions directly connected to seed metabolite
    - depth 1+: expand via metabolites of discovered reactions

    exclude_metabolites lets you skip currency-like metabolites during expansion.
    Items may be metabolite IDs or exchange IDs.
    """

    if solution_index >= len(df_fluxes):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for fluxes length {len(df_fluxes)}"
        )
    if solution_index >= len(df_shadow_prices):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for shadow prices length {len(df_shadow_prices)}"
        )
    if not isinstance(search_depth, int) or search_depth < 0:
        raise ValueError("search_depth must be a non-negative integer")
    if flux_threshold < 0:
        raise ValueError("flux_threshold must be non-negative")

    seed_metabolite = _resolve_metabolite_or_exchange(metabolite_name, model)
    excluded_met_ids = _normalize_excluded_metabolites(exclude_metabolites, model)
    excluded_met_ids.discard(seed_metabolite.id)

    flux_row = df_fluxes.iloc[solution_index]

    def get_flux_value(reaction_id):
        val = flux_row.get(reaction_id, 0.0)
        try:
            return float(val)
        except (TypeError, ValueError):
            return 0.0

    reaction_ids = {rxn.id for rxn in seed_metabolite.reactions}
    visited_metabolite_ids = {seed_metabolite.id}
    frontier_metabolite_ids = {seed_metabolite.id}
    reaction_depth = {}
    for rid in reaction_ids:
        reaction_depth[rid] = 0

    for depth in range(1, search_depth + 1):
        next_frontier_metabolites = set()
        for met_id in frontier_metabolite_ids:
            if met_id in excluded_met_ids:
                continue
            met_obj = model.metabolites.get_by_id(met_id)
            for rxn in met_obj.reactions:
                reaction_ids.add(rxn.id)
                if rxn.id not in reaction_depth:
                    reaction_depth[rxn.id] = depth
                for met2 in rxn.metabolites:
                    if met2.id in excluded_met_ids:
                        continue
                    if met2.id not in visited_metabolite_ids:
                        visited_metabolite_ids.add(met2.id)
                        next_frontier_metabolites.add(met2.id)
        frontier_metabolite_ids = next_frontier_metabolites
        if not frontier_metabolite_ids:
            break

    filtered_reaction_ids = set(reaction_ids)
    if nonzero_only:
        filtered_reaction_ids = {
            rid
            for rid in reaction_ids
            if abs(get_flux_value(rid)) > flux_threshold
        }
        if keep_seed_metabolite:
            seed_rxn_ids = {rxn.id for rxn in seed_metabolite.reactions}
            filtered_reaction_ids.update(seed_rxn_ids)

    reaction_objects = [model.reactions.get_by_id(rid) for rid in sorted(filtered_reaction_ids)]
    if not reaction_objects:
        raise ValueError(
            "No reactions left to plot after applying nonzero flux filter. "
            "Lower flux_threshold, disable nonzero_only, or change solution_index."
        )

    if map_json_path is None:
        map_json_path = f"iML1515.met.{seed_metabolite.id}.depth{search_depth}.json"

    map_name = f"iML1515.met.{seed_metabolite.id}.depth{search_depth}"
    map_id = f"iML1515_met_{seed_metabolite.id}_d{search_depth}"
    map_description = (
        f"Auto-generated Escher map for metabolite {seed_metabolite.id} with search_depth={search_depth}"
    )

    builder, map_stats = _build_escher_map_and_builder_from_reactions(
        reaction_objects=reaction_objects,
        map_name=map_name,
        map_id=map_id,
        map_description=map_description,
        map_json_path=map_json_path,
        model=model,
        reaction_data=df_fluxes.iloc[solution_index],
        metabolite_data=df_shadow_prices.iloc[solution_index],
        hide_secondary_metabolites=hide_secondary_metabolites,
        hide_all_labels=hide_all_labels,
        reaction_scale_preset=reaction_scale_preset,
    )

    if verbose:
        depth_counts = {}
        for rid, dep in reaction_depth.items():
            depth_counts[dep] = depth_counts.get(dep, 0) + 1
        print("=" * 60)
        print(f"Seed metabolite: {seed_metabolite.id}")
        print(f"Requested search_depth: {search_depth}")
        print(f"Nonzero only: {nonzero_only} (threshold={flux_threshold})")
        print(f"Excluded metabolites: {len(excluded_met_ids)}")
        print(f"Discovered reactions (topology): {len(reaction_ids)}")
        print(f"Plotted reactions (after filter): {map_stats['reaction_count']}")
        print(f"Discovered metabolites (topology): {len(visited_metabolite_ids)}")
        print(f"Discovered metabolites (map unique): {map_stats['metabolite_count']}")
        print(f"Map output: {map_stats['map_json_path']}")
        print("Reactions per depth:")
        for dep in sorted(depth_counts):
            print(f"  depth {dep}: {depth_counts[dep]}")
        print("=" * 60)

    return builder


def build_path_between_metabolites_escher_builder(
    start_metabolite_name,
    end_metabolite_name,
    df_fluxes,
    df_shadow_prices,
    model,
    solution_index=0,
    nonzero_only=False,
    flux_threshold=1e-9,
    keep_endpoints=True,
    exclude_metabolites=None,
    map_json_path=None,
    hide_secondary_metabolites=True,
    hide_all_labels=False,
    reaction_scale_preset="GaBuRd",
    verbose=True,
):
    from collections import deque

    if solution_index >= len(df_fluxes):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for fluxes length {len(df_fluxes)}"
        )
    if solution_index >= len(df_shadow_prices):
        raise IndexError(
            f"solution_index {solution_index} is out of bounds for shadow prices length {len(df_shadow_prices)}"
        )
    if flux_threshold < 0:
        raise ValueError("flux_threshold must be non-negative")

    start_met = _resolve_metabolite_or_exchange(start_metabolite_name, model)
    end_met = _resolve_metabolite_or_exchange(end_metabolite_name, model)
    endpoint_met_ids = {start_met.id, end_met.id}

    excluded_met_ids = _normalize_excluded_metabolites(exclude_metabolites, model)
    excluded_met_ids.difference_update(endpoint_met_ids)

    flux_row = df_fluxes.iloc[solution_index]

    def get_flux_value(reaction_id):
        val = flux_row.get(reaction_id, 0.0)
        try:
            return float(val)
        except (TypeError, ValueError):
            return 0.0

    def reaction_is_allowed(reaction_id):
        if not nonzero_only:
            return True
        if keep_endpoints:
            endpoint_rxn_ids = {rxn.id for rxn in start_met.reactions}.union(
                {rxn.id for rxn in end_met.reactions}
            )
            if reaction_id in endpoint_rxn_ids:
                return True
        return abs(get_flux_value(reaction_id)) > flux_threshold

    start_node = ("met", start_met.id)
    queue = deque([start_node])
    visited_nodes = {start_node}
    parent = {start_node: None}
    found_end_node = None

    while queue:
        node_type, node_id = queue.popleft()

        if node_type == "met" and node_id == end_met.id:
            found_end_node = (node_type, node_id)
            break

        if node_type == "met":
            if node_id in excluded_met_ids:
                continue
            met_obj = model.metabolites.get_by_id(node_id)
            for rxn in sorted(met_obj.reactions, key=lambda r: r.id):
                if not reaction_is_allowed(rxn.id):
                    continue
                nxt = ("rxn", rxn.id)
                if nxt not in visited_nodes:
                    visited_nodes.add(nxt) # type: ignore
                    parent[nxt] = (node_type, node_id)  # type: ignore
                    queue.append(nxt)  # type: ignore
        else:
            rxn_obj = model.reactions.get_by_id(node_id)
            for met in sorted(rxn_obj.metabolites, key=lambda m: m.id):
                if met.id in excluded_met_ids:
                    continue
                nxt = ("met", met.id)
                if nxt not in visited_nodes:
                    visited_nodes.add(nxt)
                    parent[nxt] = (node_type, node_id)  # type: ignore
                    queue.append(nxt)

    if found_end_node is None:
        raise ValueError(
            f"No path found between {start_met.id} and {end_met.id} "
            f"under current filter settings (nonzero_only={nonzero_only}, flux_threshold={flux_threshold}, excluded_mets={len(excluded_met_ids)})."
        )

    path_nodes = []
    cur = found_end_node
    while cur is not None:
        path_nodes.append(cur)
        cur = parent[cur]
    path_nodes.reverse()

    reaction_path_ids = [node_id for node_type, node_id in path_nodes if node_type == "rxn"]
    reaction_objects = [model.reactions.get_by_id(rid) for rid in reaction_path_ids]

    if not reaction_objects:
        raise ValueError("No reactions found on the metabolite path.")

    if map_json_path is None:
        map_json_path = f"iML1515.metpath.{start_met.id}_to_{end_met.id}.json"

    map_name = f"iML1515.metpath.{start_met.id}_to_{end_met.id}"
    map_id = f"iML1515_metpath_{start_met.id}_to_{end_met.id}"
    map_description = f"Auto-generated Escher metabolite path map from {start_met.id} to {end_met.id}"

    builder, map_stats = _build_escher_map_and_builder_from_reactions(
        reaction_objects=reaction_objects,
        map_name=map_name,
        map_id=map_id,
        map_description=map_description,
        map_json_path=map_json_path,
        model=model,
        reaction_data=df_fluxes.iloc[solution_index],
        metabolite_data=df_shadow_prices.iloc[solution_index],
        hide_secondary_metabolites=hide_secondary_metabolites,
        hide_all_labels=hide_all_labels,
        reaction_scale_preset=reaction_scale_preset,
    )

    if verbose:
        print("=" * 60)
        print(f"Metabolite path search: {start_met.id} -> {end_met.id}")
        print(f"Path length (reactions): {len(reaction_path_ids)}")
        print("Reaction path:")
        print("  " + " -> ".join(reaction_path_ids))
        print(f"Nonzero only: {nonzero_only} (threshold={flux_threshold})")
        print(f"Excluded metabolites: {len(excluded_met_ids)}")
        print(f"Discovered metabolites (map unique): {map_stats['metabolite_count']}")
        print(f"Map output: {map_stats['map_json_path']}")
        print("=" * 60)

    return builder
