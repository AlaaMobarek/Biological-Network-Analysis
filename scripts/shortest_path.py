import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# Import graph loading function
from load_graph import load_graph


# =============================================================================
# STEP 1: Find ALL Acyclic Shortest Paths Using Yen's Algorithm
# =============================================================================

def find_all_shortest_paths(G, source, target):
    """
    Find ALL simple (acyclic) paths that share the minimum total cost.
    """

    # Validate source and target existence
    if source not in G:
        raise ValueError(f"Source '{source}' not found in the network.")
    if target not in G:
        raise ValueError(f"Target '{target}' not found in the network.")

    # Use Dijkstra to compute the minimum path cost
    try:
        min_cost = nx.shortest_path_length(
            G,
            source,
            target,
            weight='cost'
        )

    except nx.NetworkXNoPath:
        raise nx.NetworkXNoPath(
            f"No path between '{source}' and '{target}'."
        )

    print(
        f"[Dijkstra]  Minimum internal cost "
        f"(sum of 1-conf) = {min_cost:.6f}"
    )

    # Collect all shortest simple paths with the same minimum cost
    shortest_paths = []
    tolerance = 1e-9

    for path in nx.shortest_simple_paths(
        G,
        source,
        target,
        weight='cost'
    ):

        path_cost = sum(
            G[path[i]][path[i + 1]]['cost']
            for i in range(len(path) - 1)
        )

        if abs(path_cost - min_cost) < tolerance:
            shortest_paths.append(path)
        else:
            break

    # Compute total confidence score
    total_score = sum(
        G[shortest_paths[0][i]][shortest_paths[0][i + 1]]['confidence']
        for i in range(len(shortest_paths[0]) - 1)
    )

    print(f"[Yen's]     {len(shortest_paths)} path(s) found.")
    print(
        f"[INFO]      Total Path Score "
        f"(sum of confidence) = {total_score:.6f}"
    )

    return shortest_paths, min_cost, total_score


# =============================================================================
# STEP 2: Write Results to Text File
# =============================================================================

def write_paths_to_file(
    G,
    paths,
    min_cost,
    total_score,
    source,
    target,
    output_file
):
    """
    Save all shortest paths and statistics into a text report.
    """

    # Create output directory if it does not exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w', encoding='utf-8') as f:

        # Report header
        f.write("=" * 70 + "\n")
        f.write("  ACYCLIC SHORTEST PATH(S) ANALYSIS\n")
        f.write("=" * 70 + "\n")

        f.write(f"  Source Protein                        : {source}\n")
        f.write(f"  Target Protein                        : {target}\n")
        f.write(f"  Total Shortest Paths Found            : {len(paths)}\n")
        f.write(f"\n")

        f.write(
            f"  Total Path Score "
            f"(sum of confidence)  : {total_score:.6f}\n"
        )

        f.write(
            f"  [Internal cost  "
            f"(sum of 1-conf)       : {min_cost:.6f}]\n"
        )

        f.write("=" * 70 + "\n")

        # Write every shortest path
        for idx, path in enumerate(paths, start=1):

            f.write("\n\n")
            f.write("*" * 70 + "\n")
            f.write(f"  PATH {idx} of {len(paths)}\n")
            f.write("*" * 70 + "\n")

            f.write(
                f"\n  Sequence "
                f"({len(path)} nodes, {len(path)-1} edges):\n\n"
            )

            f.write(f"      {' -> '.join(path)}\n")
            f.write("\n")

            path_score = 0.0

            f.write(
                f"  {'Edge':<45} "
                f"{'Weight (Confidence)':>20}\n"
            )

            f.write(
                f"  {'-'*45} "
                f"{'-'*20}\n"
            )

            # Write edge-by-edge confidence
            for i in range(len(path) - 1):

                u, v = path[i], path[i + 1]
                conf = G[u][v]['confidence']

                f.write(
                    f"  {f'{u} -> {v}':<45} "
                    f"{conf:>20.6f}\n"
                )

                path_score += conf

            avg_conf = path_score / (len(path) - 1)

            f.write("\n")

            f.write(
                f"  {'Total Path Score (sum of confidence)':<45} "
                f"{path_score:>20.6f}\n"
            )

            f.write(
                f"  {'Average Confidence per Edge':<45} "
                f"{avg_conf:>20.6f}\n"
            )

        # Report footer
        f.write("\n" + "=" * 70 + "\n")
        f.write("  END OF REPORT\n")
        f.write("=" * 70 + "\n")


# =============================================================================
# STEP 3: Draw the Sub-Network
# =============================================================================

def draw_subnetwork(
    G,
    paths,
    source,
    target,
    total_score,
    output_image
):
    """
    Draw the shortest-path subnetwork and save it as an image.
    """

    # Build subgraph from all shortest paths
    sub_G = nx.DiGraph()

    for path in paths:

        for i in range(len(path) - 1):

            u, v = path[i], path[i + 1]

            sub_G.add_edge(
                u,
                v,
                confidence=G[u][v]['confidence'],
                cost=G[u][v]['cost']
            )

    # Count intermediate nodes
    n_mid = len([
        n for n in sub_G.nodes()
        if n != source and n != target
    ])

    # Generate node positions
    pos = nx.spring_layout(
        sub_G,
        k=2.5,
        iterations=200,
        seed=42
    )

    # Define node appearance
    node_colors = []
    node_sizes = []

    for n in sub_G.nodes():

        if n == source:
            node_colors.append('#00FF7F')
            node_sizes.append(3500)

        elif n == target:
            node_colors.append('#FF4500')
            node_sizes.append(3500)

        else:
            node_colors.append('#00BFFF')
            node_sizes.append(2500)

    # Edge confidence values
    confidences = [
        sub_G[u][v]['confidence']
        for u, v in sub_G.edges()
    ]

    edge_colors = plt.cm.plasma(np.array(confidences))

    # Edge labels
    edge_labels = {
        (u, v): f"{sub_G[u][v]['confidence']:.2f}"
        for u, v in sub_G.edges()
    }

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 12))

    fig.patch.set_facecolor('#121212')
    ax.set_facecolor('#121212')

    # Draw edges
    nx.draw_networkx_edges(
        sub_G,
        pos,
        ax=ax,
        edge_color=edge_colors,
        edge_cmap=plt.cm.plasma,
        arrows=True,
        arrowsize=30,
        arrowstyle='-|>',
        width=1.2,
        alpha=0.9,
        connectionstyle='arc3,rad=0.2',
        min_source_margin=25,
        min_target_margin=25,
    )

    # Draw nodes
    nx.draw_networkx_nodes(
        sub_G,
        pos,
        ax=ax,
        node_color=node_colors,
        node_size=node_sizes,
        alpha=0.9,
        linewidths=3,
        edgecolors='#FFFFFF'
    )

    # Draw node labels
    nx.draw_networkx_labels(
        sub_G,
        pos,
        ax=ax,
        font_size=11,
        font_color='black',
        font_weight='bold'
    )

    # Draw edge labels
    nx.draw_networkx_edge_labels(
        sub_G,
        pos,
        edge_labels=edge_labels,
        ax=ax,
        font_size=9,
        font_color='#FFFFFF',
        bbox=dict(
            boxstyle='round,pad=0.3',
            fc='#2E2E2E',
            alpha=0.8,
            ec='none'
        ),
        label_pos=0.5,
    )

    # Legend
    ax.legend(
        handles=[
            mpatches.Patch(
                color='#00FF7F',
                label=f'Source: {source}'
            ),

            mpatches.Patch(
                color='#FF4500',
                label=f'Target: {target}'
            ),

            mpatches.Patch(
                color='#00BFFF',
                label=f'Intermediate nodes ({n_mid})'
            ),
        ],

        loc='upper left',
        facecolor='#1E1E1E',
        edgecolor='gray',
        labelcolor='white',
        fontsize=12
    )

    # Colorbar for confidence values
    sm = plt.cm.ScalarMappable(
        cmap='plasma',
        norm=plt.Normalize(0, 1)
    )

    sm.set_array([])

    cbar = plt.colorbar(
        sm,
        ax=ax,
        shrink=0.5,
        pad=0.02
    )

    cbar.set_label(
        'Interaction Confidence',
        color='white',
        fontsize=12,
        fontweight='bold'
    )

    cbar.ax.yaxis.set_tick_params(color='white')

    plt.setp(
        plt.getp(cbar.ax.axes, 'yticklabels'),
        color='white'
    )

    # Figure title
    ax.set_title(
        f"Protein-Protein Interaction: "
        f"Acyclic Shortest Paths Network\n"
        f"{source}  →  {target}  |  "
        f"{len(paths)} path(s)  |  "
        f"Total Path Score: {total_score:.4f}",

        color='white',
        fontsize=16,
        fontweight='bold',
        pad=25
    )

    ax.axis('off')

    plt.tight_layout()

    # Create output folder
    os.makedirs(
        os.path.dirname(output_image),
        exist_ok=True
    )

    # Save figure
    plt.savefig(
        output_image,
        dpi=300,
        bbox_inches='tight',
        facecolor=fig.get_facecolor()
    )

    # Close figure to free memory
    plt.close(fig)


# =============================================================================
# STEP 4: Complete Analysis Pipeline
# =============================================================================

def run_shortest_paths_pipeline(
    G,
    source,
    target,
    output_txt,
    output_img
):
    """
    Run the complete shortest-path analysis pipeline:

    1. Find shortest paths
    2. Save text report
    3. Draw and save subnetwork
    """

    print(f"\n[INFO] Analyzing: {source} -> {target}")

    try:

        # Find shortest paths
        paths, min_cost, total_score = find_all_shortest_paths(
            G,
            source,
            target
        )

        # Save text report
        write_paths_to_file(
            G,
            paths,
            min_cost,
            total_score,
            source,
            target,
            output_txt
        )

        print(f"[INFO] Text report saved -> {output_txt}")

        # Draw network figure
        draw_subnetwork(
            G,
            paths,
            source,
            target,
            total_score,
            output_img
        )

        print(f"[INFO] Figure saved -> {output_img}")

    except (ValueError, nx.NetworkXNoPath) as e:
        print(f"[ERROR] {e}")