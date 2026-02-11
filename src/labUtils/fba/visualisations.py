import json
import requests
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
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
    ax1.set_title(f"{title}{subtitle}"
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
