import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# استدعاء دالة تحميل الجراف من الملف الخارجي
from load_graph import load_graph

# =============================================================================
# STEP 1: Find ALL Acyclic Shortest Paths Using Yen's Algorithm
# =============================================================================

def find_all_shortest_paths(G, source, target):
    """
    Find ALL simple (acyclic) paths that share the minimum total cost.
    """
    if source not in G:
        raise ValueError(f"Source '{source}' not found in the network.")
    if target not in G:
        raise ValueError(f"Target '{target}' not found in the network.")

    # Step 1: Dijkstra -> get minimum cost using internal cost attribute
    try:
        min_cost = nx.shortest_path_length(G, source, target, weight='cost')
    except nx.NetworkXNoPath:
        raise nx.NetworkXNoPath(f"No path between '{source}' and '{target}'.")

    print(f"[Dijkstra]  Minimum internal cost (sum of 1-conf) = {min_cost:.6f}")

    # Step 2: Yen's -> collect all paths that equal min_cost
    shortest_paths = []
    tolerance = 1e-9

    for path in nx.shortest_simple_paths(G, source, target, weight='cost'):
        path_cost = sum(
            G[path[i]][path[i + 1]]['cost']
            for i in range(len(path) - 1)
        )
        if abs(path_cost - min_cost) < tolerance:
            shortest_paths.append(path)
        else:
            break

    # Total path score = sum of confidence (same for all tied paths)
    total_score = sum(
        G[shortest_paths[0][i]][shortest_paths[0][i + 1]]['confidence']
        for i in range(len(shortest_paths[0]) - 1)
    )

    print(f"[Yen's]     {len(shortest_paths)} path(s) found.")
    print(f"[INFO]      Total Path Score (sum of confidence) = {total_score:.6f}")
    return shortest_paths, min_cost, total_score


# =============================================================================
# STEP 2: Write Results to Text File
# =============================================================================

def write_paths_to_file(G, paths, min_cost, total_score, source, target, output_file):
    """
    Save every shortest path to a text file.
    Reported values use CONFIDENCE as the edge weight.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("  ACYCLIC SHORTEST PATH(S) ANALYSIS\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Source Protein                        : {source}\n")
        f.write(f"  Target Protein                        : {target}\n")
        f.write(f"  Total Shortest Paths Found            : {len(paths)}\n")
        f.write(f"\n")
        f.write(f"  Total Path Score (sum of confidence)  : {total_score:.6f}\n")
        f.write(f"  [Internal cost  (sum of 1-conf)       : {min_cost:.6f}]\n")
        f.write("=" * 70 + "\n")

        for idx, path in enumerate(paths, start=1):
            f.write("\n\n")
            f.write("*" * 70 + "\n")
            f.write(f"  PATH {idx} of {len(paths)}\n")
            f.write("*" * 70 + "\n")
            f.write(f"\n  Sequence ({len(path)} nodes, {len(path)-1} edges):\n\n")
            f.write(f"      {' -> '.join(path)}\n")
            f.write("\n")

            path_score = 0.0
            f.write(f"  {'Edge':<45} {'Weight (Confidence)':>20}\n")
            f.write(f"  {'-'*45} {'-'*20}\n")

            for i in range(len(path) - 1):
                u, v  = path[i], path[i + 1]
                conf  = G[u][v]['confidence']
                f.write(f"  {f'{u} -> {v}':<45} {conf:>20.6f}\n")
                path_score += conf

            avg_conf = path_score / (len(path) - 1)
            f.write("\n")
            f.write(f"  {'Total Path Score (sum of confidence)':<45} {path_score:>20.6f}\n")
            f.write(f"  {'Average Confidence per Edge':<45} {avg_conf:>20.6f}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("  END OF REPORT\n")
        f.write("=" * 70 + "\n")


# =============================================================================
# STEP 3: Draw the Sub-Network (Thinner Edges & Prominent Arrows)
# =============================================================================

def draw_subnetwork(G, paths, source, target, total_score, output_image):
    """
    Draw the sub-network with thinner edges and prominent arrows, 
    and save it without showing the plot window.
    """
    sub_G = nx.DiGraph()
    for path in paths:
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            sub_G.add_edge(u, v,
                           confidence=G[u][v]['confidence'],
                           cost=G[u][v]['cost'])

    n_mid = len([n for n in sub_G.nodes() if n != source and n != target])

    pos = nx.spring_layout(sub_G, k=2.5, iterations=200, seed=42)

    node_colors, node_sizes = [], []
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

    confidences = [sub_G[u][v]['confidence'] for u, v in sub_G.edges()]
    edge_colors = plt.cm.plasma(np.array(confidences))

    edge_labels = {
        (u, v): f"{sub_G[u][v]['confidence']:.2f}"
        for u, v in sub_G.edges()
    }

    fig, ax = plt.subplots(figsize=(16, 12))
    fig.patch.set_facecolor('#121212')
    ax.set_facecolor('#121212')

    nx.draw_networkx_edges(
        sub_G, pos, ax=ax,
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

    nx.draw_networkx_nodes(
        sub_G, pos, ax=ax,
        node_color=node_colors, node_size=node_sizes,
        alpha=0.9, linewidths=3, edgecolors='#FFFFFF'
    )

    nx.draw_networkx_labels(
        sub_G, pos, ax=ax,
        font_size=11, font_color='black', font_weight='bold'
    )

    nx.draw_networkx_edge_labels(
        sub_G, pos, edge_labels=edge_labels, ax=ax,
        font_size=9,
        font_color='#FFFFFF',
        bbox=dict(boxstyle='round,pad=0.3', fc='#2E2E2E', alpha=0.8, ec='none'),
        label_pos=0.5,
    )

    ax.legend(handles=[
        mpatches.Patch(color='#00FF7F', label=f'Source: {source}'),
        mpatches.Patch(color='#FF4500', label=f'Target: {target}'),
        mpatches.Patch(color='#00BFFF', label=f'Intermediate nodes ({n_mid})'),
    ], loc='upper left', facecolor='#1E1E1E',
       edgecolor='gray', labelcolor='white', fontsize=12)

    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('Interaction Confidence', color='white', fontsize=12, fontweight='bold')
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

    ax.set_title(
        f"Protein-Protein Interaction: Acyclic Shortest Paths Network\n"
        f"{source}  →  {target}  |  "
        f"{len(paths)} path(s)  |  Total Path Score: {total_score:.4f}",
        color='white', fontsize=16, fontweight='bold', pad=25
    )

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    
    # تم إغلاق الـ plot لتفادي استهلاك الذاكرة وعدم عرض النافذة
    plt.close(fig)


# =============================================================================
# STEP 4: The Pipeline Function
# =============================================================================

def run_shortest_paths_pipeline(G, source, target, output_txt, output_img):
    """
    دالة بتجمع بايب لاين تحليل المسارات الأقصر كلها:
    1- إيجاد المسارات.
    2- كتابة التقرير النصي.
    3- رسم وحفظ الشبكة.
    """
    print(f"\n[INFO] Analyzing: {source} -> {target}")
    try:
        paths, min_cost, total_score = find_all_shortest_paths(G, source, target)
        
        write_paths_to_file(G, paths, min_cost, total_score, source, target, output_txt)
        print(f"[INFO] Text report saved -> {output_txt}")

        draw_subnetwork(G, paths, source, target, total_score, output_img)
        print(f"[INFO] Figure saved -> {output_img}")

    except (ValueError, nx.NetworkXNoPath) as e:
        print(f"[ERROR] {e}")


# =============================================================================
# MAIN
# =============================================================================

def main():

    INTERACTOME_FILE = "C:\\Users\\nadah\\Documents\\final_bioinfo\\Biological-Network-Analysis\\data\\PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    OUTPUT_TXT       = "shortest_paths_result.txt"
    OUTPUT_IMG       = "shortest_paths_subnetwork.png"

    SOURCE_PROTEIN = "P04637"   # TP53
    TARGET_PROTEIN = "P00533"   # EGFR

    # 1. Load Graph using the external module
    if not os.path.exists(INTERACTOME_FILE):
        print(f"[ERROR] File not found: {INTERACTOME_FILE}")
        return
        
    print(f"[INFO] Loading graph from {INTERACTOME_FILE}...")
    G = load_graph(INTERACTOME_FILE)

    # 2. Execute the Pipeline
    run_shortest_paths_pipeline(
        G=G,
        source=SOURCE_PROTEIN,
        target=TARGET_PROTEIN,
        output_txt=OUTPUT_TXT,
        output_img=OUTPUT_IMG
    )

    print("\n[DONE] Pipeline execution finished.")


if __name__ == "__main__":
    main()
# import networkx as nx
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import numpy as np
# import os


# # =============================================================================
# # STEP 2: Find ALL Acyclic Shortest Paths Using Yen's Algorithm
# # =============================================================================

# def find_all_shortest_paths(G, source, target):
#     """
#     Find ALL simple (acyclic) paths that share the minimum total cost.

#     Algorithm
#     ---------
#     Step 1 - Dijkstra's Algorithm:
#         Uses edge attribute 'cost' (= 1 - confidence) internally.
#         -> gives us min_cost (threshold for filtering)

#     Step 2 - Yen's K-Shortest Paths:
#         nx.shortest_simple_paths() with weight='cost'.
#         Returns paths ordered from lowest to highest cost.
#         We collect all paths whose cost == min_cost and break
#         immediately when a path exceeds it.

#     Returns
#     -------
#     shortest_paths : list of lists  - all paths with cost == min_cost
#     min_cost       : float          - the minimum total cost (sum of 1-conf)
#     total_score    : float          - total path score (sum of confidence)
#     """

#     if source not in G:
#         raise ValueError(f"Source '{source}' not found in the network.")
#     if target not in G:
#         raise ValueError(f"Target '{target}' not found in the network.")

#     # Step 1: Dijkstra -> get minimum cost using internal cost attribute
#     try:
#         min_cost = nx.shortest_path_length(G, source, target, weight='cost')
#     except nx.NetworkXNoPath:
#         raise nx.NetworkXNoPath(f"No path between '{source}' and '{target}'.")

#     print(f"[Dijkstra]  Minimum internal cost (sum of 1-conf) = {min_cost:.6f}")

#     # Step 2: Yen's -> collect all paths that equal min_cost
#     shortest_paths = []
#     tolerance = 1e-9

#     for path in nx.shortest_simple_paths(G, source, target, weight='cost'):
#         path_cost = sum(
#             G[path[i]][path[i + 1]]['cost']
#             for i in range(len(path) - 1)
#         )
#         if abs(path_cost - min_cost) < tolerance:
#             shortest_paths.append(path)
#         else:
#             break

#     # Total path score = sum of confidence (same for all tied paths)
#     total_score = sum(
#         G[shortest_paths[0][i]][shortest_paths[0][i + 1]]['confidence']
#         for i in range(len(shortest_paths[0]) - 1)
#     )

#     print(f"[Yen's]     {len(shortest_paths)} path(s) found.")
#     print(f"[INFO]      Total Path Score (sum of confidence) = {total_score:.6f}")
#     return shortest_paths, min_cost, total_score


# # =============================================================================
# # STEP 3: Write Results to Text File
# # =============================================================================

# def write_paths_to_file(G, paths, min_cost, total_score, source, target, output_file):
#     """
#     Save every shortest path to a text file.
#     Reported values use CONFIDENCE as the edge weight (as given in the file).
#     Paths are separated clearly for readability.
#     """
#     with open(output_file, 'w', encoding='utf-8') as f:
#         f.write("=" * 70 + "\n")
#         f.write("  ACYCLIC SHORTEST PATH(S) ANALYSIS\n")
#         f.write("=" * 70 + "\n")
#         f.write(f"  Source Protein                        : {source}\n")
#         f.write(f"  Target Protein                        : {target}\n")
#         f.write(f"  Total Shortest Paths Found            : {len(paths)}\n")
#         f.write(f"\n")
#         f.write(f"  Total Path Score (sum of confidence)  : {total_score:.6f}\n")
#         f.write(f"  [Internal cost  (sum of 1-conf)       : {min_cost:.6f}]\n")
#         f.write("=" * 70 + "\n")

#         for idx, path in enumerate(paths, start=1):
#             # Extra blank lines between paths for readability
#             f.write("\n\n")
#             f.write("*" * 70 + "\n")
#             f.write(f"  PATH {idx} of {len(paths)}\n")
#             f.write("*" * 70 + "\n")
#             f.write(f"\n  Sequence ({len(path)} nodes, {len(path)-1} edges):\n\n")
#             f.write(f"      {' -> '.join(path)}\n")
#             f.write("\n")

#             path_score = 0.0

#             f.write(f"  {'Edge':<45} {'Weight (Confidence)':>20}\n")
#             f.write(f"  {'-'*45} {'-'*20}\n")

#             for i in range(len(path) - 1):
#                 u, v  = path[i], path[i + 1]
#                 conf  = G[u][v]['confidence']
#                 f.write(f"  {f'{u} -> {v}':<45} {conf:>20.6f}\n")
#                 path_score += conf

#             avg_conf = path_score / (len(path) - 1)
#             f.write("\n")
#             f.write(f"  {'Total Path Score (sum of confidence)':<45} {path_score:>20.6f}\n")
#             f.write(f"  {'Average Confidence per Edge':<45} {avg_conf:>20.6f}\n")

#         f.write("\n" + "=" * 70 + "\n")
#         f.write("  END OF REPORT\n")
#         f.write("=" * 70 + "\n")

#     print(f"[INFO] Text results saved -> {output_file}")


# # =============================================================================
# # STEP 4: Draw the Sub-Network (Thinner Edges & Prominent Arrows)
# # =============================================================================

# def draw_subnetwork(G, paths, source, target, total_score, output_image):
#     """
#     Draw the sub-network with thinner edges and prominent arrows.
#     """
#     # Build sub-graph
#     sub_G = nx.DiGraph()
#     for path in paths:
#         for i in range(len(path) - 1):
#             u, v = path[i], path[i + 1]
#             sub_G.add_edge(u, v,
#                            confidence=G[u][v]['confidence'],
#                            cost=G[u][v]['cost'])

#     n_mid = len([n for n in sub_G.nodes() if n != source and n != target])

#     # ── Positions (Force-Directed with spacing) ──────────────────────────────
#     pos = nx.spring_layout(sub_G, k=2.5, iterations=200, seed=42)

#     # ── Node colors and sizes ─────────────────────────────────────────────────
#     node_colors, node_sizes = [], []
#     for n in sub_G.nodes():
#         if n == source:
#             node_colors.append('#00FF7F')
#             node_sizes.append(3500)
#         elif n == target:
#             node_colors.append('#FF4500')
#             node_sizes.append(3500)
#         else:
#             node_colors.append('#00BFFF')
#             node_sizes.append(2500)

#     # ── Edge colors by confidence ─────────────────────────────────────────────
#     confidences = [sub_G[u][v]['confidence'] for u, v in sub_G.edges()]
#     edge_colors = plt.cm.plasma(np.array(confidences))

#     # ── Edge labels: show confidence value on every edge ──────────────────────
#     edge_labels = {
#         (u, v): f"{sub_G[u][v]['confidence']:.2f}"
#         for u, v in sub_G.edges()
#     }

#     # ── Figure ────────────────────────────────────────────────────────────────
#     fig, ax = plt.subplots(figsize=(16, 12))
#     fig.patch.set_facecolor('#121212')
#     ax.set_facecolor('#121212')

#     # Draw directed edges (Thinner lines, bigger arrows)
#     nx.draw_networkx_edges(
#         sub_G, pos, ax=ax,
#         edge_color=edge_colors,
#         edge_cmap=plt.cm.plasma,
#         arrows=True,
#         arrowsize=30,     
#         arrowstyle='-|>',
#         width=1.2,         
#         alpha=0.9,
#         connectionstyle='arc3,rad=0.2',
#         min_source_margin=25,
#         min_target_margin=25,
#     )

#     # Draw nodes
#     nx.draw_networkx_nodes(
#         sub_G, pos, ax=ax,
#         node_color=node_colors, node_size=node_sizes,
#         alpha=0.9, linewidths=3, edgecolors='#FFFFFF'
#     )

#     # Node labels
#     nx.draw_networkx_labels(
#         sub_G, pos, ax=ax,
#         font_size=11, font_color='black', font_weight='bold'
#     )

#     # Confidence label on every edge
#     nx.draw_networkx_edge_labels(
#         sub_G, pos, edge_labels=edge_labels, ax=ax,
#         font_size=9,
#         font_color='#FFFFFF',
#         bbox=dict(boxstyle='round,pad=0.3', fc='#2E2E2E', alpha=0.8, ec='none'),
#         label_pos=0.5,
#     )

#     # ── Legend ────────────────────────────────────────────────────────────────
#     ax.legend(handles=[
#         mpatches.Patch(color='#00FF7F', label=f'Source: {source}'),
#         mpatches.Patch(color='#FF4500', label=f'Target: {target}'),
#         mpatches.Patch(color='#00BFFF', label=f'Intermediate nodes ({n_mid})'),
#     ], loc='upper left', facecolor='#1E1E1E',
#        edgecolor='gray', labelcolor='white', fontsize=12)

#     # ── Colorbar ──────────────────────────────────────────────────────────────
#     sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(0, 1))
#     sm.set_array([])
#     cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
#     cbar.set_label('Interaction Confidence', color='white', fontsize=12, fontweight='bold')
#     cbar.ax.yaxis.set_tick_params(color='white')
#     plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

#     ax.set_title(
#         f"Protein-Protein Interaction: Acyclic Shortest Paths Network\n"
#         f"{source}  →  {target}  |  "
#         f"{len(paths)} path(s)  |  Total Path Score: {total_score:.4f}",
#         color='white', fontsize=16, fontweight='bold', pad=25
#     )

#     ax.axis('off')
#     plt.tight_layout()
#     plt.savefig(output_image, dpi=300, bbox_inches='tight',
#                 facecolor=fig.get_facecolor())
#     plt.show()
#     print(f"[INFO] Figure saved -> {output_image}")

# # =============================================================================
# # MAIN
# # =============================================================================

# def main():

#     INTERACTOME_FILE = "C:\\Users\\nadah\\Documents\\final_bioinfo\\Biological-Network-Analysis\\data\\PathLinker_2018_human-ppi-weighted-cap0_75.txt"
#     OUTPUT_TXT       = "shortest_paths_result.txt"
#     OUTPUT_IMG       = "shortest_paths_subnetwork.png"

#     # Change these to your two proteins of choice
#     SOURCE_PROTEIN = "P04637"   # TP53
#     TARGET_PROTEIN = "P00533"   # EGFR

#     # 1. Load
#     if not os.path.exists(INTERACTOME_FILE):
#         print(f"[ERROR] File not found: {INTERACTOME_FILE}")
#         return
#     G = load_interactome(INTERACTOME_FILE)

#     # 2. Find all shortest paths (Dijkstra + Yen's) using cost = 1-confidence
#     print(f"\n[INFO] Analyzing: {SOURCE_PROTEIN} -> {TARGET_PROTEIN}")
#     try:
#         paths, min_cost, total_score = find_all_shortest_paths(
#             G, SOURCE_PROTEIN, TARGET_PROTEIN)
#     except (ValueError, nx.NetworkXNoPath) as e:
#         print(f"[ERROR] {e}")
#         return

#     # 3. Write text file — reports confidence as the weight
#     write_paths_to_file(G, paths, min_cost, total_score,
#                         SOURCE_PROTEIN, TARGET_PROTEIN, OUTPUT_TXT)

#     # 4. Draw directed figure with edge weight labels
#     draw_subnetwork(G, paths, SOURCE_PROTEIN, TARGET_PROTEIN,
#                     total_score, OUTPUT_IMG)

#     print(f"\n[DONE]  {OUTPUT_TXT}  |  {OUTPUT_IMG}")


# if __name__ == "__main__":
#     main()




    
# # def draw_subnetwork(G, paths, source, target, total_score, output_image):
# #     """
# #     Draw the sub-network as a DIRECTED graph with an organic, professional layout
# #     similar to biological papers.
# #     """
# #     # Build sub-graph
# #     sub_G = nx.DiGraph()
# #     for path in paths:
# #         for i in range(len(path) - 1):
# #             u, v = path[i], path[i + 1]
# #             sub_G.add_edge(u, v,
# #                            confidence=G[u][v]['confidence'],
# #                            cost=G[u][v]['cost'])

# #     n_mid = len([n for n in sub_G.nodes() if n != source and n != target])

# #     # ── Positions (Force-Directed / Organic Layout) ───────────────────────────
# #     # Using spring_layout to simulate physical repulsion/attraction for a natural look
# #     pos = nx.spring_layout(sub_G, k=0.8, iterations=50, seed=42)

# #     # ── Node colors and sizes ─────────────────────────────────────────────────
# #     node_colors, node_sizes = [], []
# #     for n in sub_G.nodes():
# #         if n == source:
# #             node_colors.append('#2ecc71'); node_sizes.append(3500)
# #         elif n == target:
# #             node_colors.append('#e74c3c'); node_sizes.append(3500)
# #         else:
# #             node_colors.append('#3498db'); node_sizes.append(2500)

# #     # ── Edge colors by confidence ─────────────────────────────────────────────
# #     confidences = [sub_G[u][v]['confidence'] for u, v in sub_G.edges()]
# #     edge_colors = plt.cm.YlOrRd(np.array(confidences))

# #     # ── Edge labels: show confidence value on every edge ──────────────────────
# #     edge_labels = {
# #         (u, v): f"{sub_G[u][v]['confidence']:.2f}"
# #         for u, v in sub_G.edges()
# #     }

# #     # ── Figure ────────────────────────────────────────────────────────────────
# #     fig, ax = plt.subplots(figsize=(14, 10))
# #     fig.patch.set_facecolor('#1a1a2e')
# #     ax.set_facecolor('#1a1a2e')

# #     # Draw directed edges with clear arrowheads and a slight curve (organic look)
# #     nx.draw_networkx_edges(
# #         sub_G, pos, ax=ax,
# #         edge_color=edge_colors,
# #         arrows=True,
# #         arrowsize=25,
# #         arrowstyle='-|>',
# #         width=2.5,
# #         alpha=0.9,
# #         connectionstyle='arc3,rad=0.1', # Slight curve makes overlapping paths clearer
# #         min_source_margin=25,
# #         min_target_margin=25,
# #     )

# #     # Draw nodes
# #     nx.draw_networkx_nodes(
# #         sub_G, pos, ax=ax,
# #         node_color=node_colors, node_size=node_sizes,
# #         alpha=0.95, linewidths=2.5, edgecolors='white'
# #     )

# #     # Node labels
# #     nx.draw_networkx_labels(
# #         sub_G, pos, ax=ax,
# #         font_size=10, font_color='white', font_weight='bold'
# #     )

# #     # Confidence label on every edge
# #     nx.draw_networkx_edge_labels(
# #         sub_G, pos, edge_labels=edge_labels, ax=ax,
# #         font_size=8,
# #         font_color='#f0e68c',
# #         bbox=dict(boxstyle='round,pad=0.2', fc='#1a1a2e', alpha=0.7, ec='none'),
# #         label_pos=0.5,  # Centered on the curved edge
# #     )

# #     # ── Legend ────────────────────────────────────────────────────────────────
# #     ax.legend(handles=[
# #         mpatches.Patch(color='#2ecc71', label=f'Source: {source}'),
# #         mpatches.Patch(color='#e74c3c', label=f'Target: {target}'),
# #         mpatches.Patch(color='#3498db', label=f'Intermediate nodes ({n_mid})'),
# #     ], loc='upper left', facecolor='#16213e',
# #        edgecolor='white', labelcolor='white', fontsize=11)

# #     # ── Colorbar ──────────────────────────────────────────────────────────────
# #     sm = plt.cm.ScalarMappable(cmap='YlOrRd', norm=plt.Normalize(0, 1))
# #     sm.set_array([])
# #     cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
# #     cbar.set_label('Interaction Confidence', color='white', fontsize=11, fontweight='bold')
# #     cbar.ax.yaxis.set_tick_params(color='white')
# #     plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

# #     ax.set_title(
# #         f"Protein-Protein Interaction: Shortest Paths Network\n"
# #         f"{source}  →  {target}  |  "
# #         f"{len(paths)} path(s)  |  Total Path Score: {total_score:.4f}",
# #         color='white', fontsize=15, fontweight='bold', pad=25
# #     )

# #     ax.axis('off')
# #     plt.tight_layout()
# #     # Increased DPI for a high-quality "publication-ready" output
# #     plt.savefig(output_image, dpi=300, bbox_inches='tight',
# #                 facecolor=fig.get_facecolor())
# #     plt.show()
# #     print(f"[INFO] Figure saved -> {output_image}")





    
# # def draw_subnetwork(G, paths, source, target, total_score, output_image):
# #     """
# #     Draw the sub-network as a DIRECTED graph with 3-column layered layout:

# #         LEFT  (x=-1.0)         -> Source node
# #         MIDDLE (x~0, jittered) -> Intermediate nodes spread with x-jitter
# #         RIGHT (x=+1.0)         -> Target node

# #     Features:
# #         - Directed arrows on every edge
# #         - Confidence value label on every edge
# #         - Nodes spread (alternating x-jitter) so labels do not overlap
# #         - Edge color by confidence (YlOrRd colormap)
# #     """

# #     # Build sub-graph
# #     sub_G = nx.DiGraph()
# #     for path in paths:
# #         for i in range(len(path) - 1):
# #             u, v = path[i], path[i + 1]
# #             sub_G.add_edge(u, v,
# #                            confidence=G[u][v]['confidence'],
# #                            cost=G[u][v]['cost'])

# #     intermediates = [n for n in sub_G.nodes() if n != source and n != target]
# #     n_mid = len(intermediates)

# #     # ── Positions ─────────────────────────────────────────────────────────────
# #     pos = {}
# #     pos[source] = np.array([-1.0, 0.0])
# #     pos[target] = np.array([ 1.0, 0.0])

# #     if n_mid == 1:
# #         pos[intermediates[0]] = np.array([0.0, 0.0])
# #     else:
# #         y_positions = np.linspace(1.0, -1.0, n_mid)
# #         # Alternate left/right of centre so nodes and labels don't overlap
# #         x_jitter = [0.15 if i % 2 == 0 else -0.15 for i in range(n_mid)]
# #         for node, y, xj in zip(intermediates, y_positions, x_jitter):
# #             pos[node] = np.array([xj, y])

# #     # ── Node colors and sizes ─────────────────────────────────────────────────
# #     node_colors, node_sizes = [], []
# #     for n in sub_G.nodes():
# #         if n == source:
# #             node_colors.append('#2ecc71'); node_sizes.append(3500)
# #         elif n == target:
# #             node_colors.append('#e74c3c'); node_sizes.append(3500)
# #         else:
# #             node_colors.append('#3498db'); node_sizes.append(2200)

# #     # ── Edge colors by confidence ─────────────────────────────────────────────
# #     confidences = [sub_G[u][v]['confidence'] for u, v in sub_G.edges()]
# #     edge_colors = plt.cm.YlOrRd(np.array(confidences))

# #     # ── Edge labels: show confidence value on every edge ──────────────────────
# #     edge_labels = {
# #         (u, v): f"{sub_G[u][v]['confidence']:.2f}"
# #         for u, v in sub_G.edges()
# #     }

# #     # ── Figure ────────────────────────────────────────────────────────────────
# #     fig_h = max(10, n_mid * 0.6)
# #     fig, ax = plt.subplots(figsize=(18, fig_h))
# #     fig.patch.set_facecolor('#1a1a2e')
# #     ax.set_facecolor('#1a1a2e')

# #     # Draw directed edges with clear arrowheads
# #     nx.draw_networkx_edges(
# #         sub_G, pos, ax=ax,
# #         edge_color=edge_colors,
# #         arrows=True,
# #         arrowsize=20,
# #         arrowstyle='-|>',
# #         width=2.0,
# #         alpha=0.9,
# #         connectionstyle='arc3,rad=0.12',
# #         min_source_margin=25,
# #         min_target_margin=25,
# #     )

# #     # Draw nodes
# #     nx.draw_networkx_nodes(
# #         sub_G, pos, ax=ax,
# #         node_color=node_colors, node_size=node_sizes,
# #         alpha=0.95, linewidths=2, edgecolors='white'
# #     )

# #     # Node labels
# #     font_sz = max(6, 9 - n_mid // 12)
# #     nx.draw_networkx_labels(
# #         sub_G, pos, ax=ax,
# #         font_size=font_sz, font_color='white', font_weight='bold'
# #     )

# #     # Confidence label on every edge
# #     nx.draw_networkx_edge_labels(
# #         sub_G, pos, edge_labels=edge_labels, ax=ax,
# #         font_size=7,
# #         font_color='#f0e68c',
# #         bbox=dict(boxstyle='round,pad=0.15', fc='#1a1a2e', alpha=0.7),
# #         label_pos=0.35,   # closer to source so arrowhead stays visible
# #     )

# #     # ── Column headers ────────────────────────────────────────────────────────
# #     y_top = 1.12
# #     ax.text(-1.0, y_top, 'SOURCE',
# #             ha='center', va='bottom', color='#2ecc71',
# #             fontsize=12, fontweight='bold', transform=ax.transData)
# #     ax.text( 0.0, y_top, f'INTERMEDIATE PROTEINS ({n_mid})',
# #             ha='center', va='bottom', color='#3498db',
# #             fontsize=12, fontweight='bold', transform=ax.transData)
# #     ax.text( 1.0, y_top, 'TARGET',
# #             ha='center', va='bottom', color='#e74c3c',
# #             fontsize=12, fontweight='bold', transform=ax.transData)

# #     # Separator lines
# #     ax.axvline(x=-0.5, color='white', linestyle='--', alpha=0.2, linewidth=1)
# #     ax.axvline(x= 0.5, color='white', linestyle='--', alpha=0.2, linewidth=1)

# #     # ── Legend ────────────────────────────────────────────────────────────────
# #     ax.legend(handles=[
# #         mpatches.Patch(color='#2ecc71', label=f'Source: {source}'),
# #         mpatches.Patch(color='#e74c3c', label=f'Target: {target}'),
# #         mpatches.Patch(color='#3498db', label=f'Intermediate nodes ({n_mid})'),
# #     ], loc='lower left', facecolor='#16213e',
# #        edgecolor='white', labelcolor='white', fontsize=10)

# #     # ── Colorbar ──────────────────────────────────────────────────────────────
# #     sm = plt.cm.ScalarMappable(cmap='YlOrRd', norm=plt.Normalize(0, 1))
# #     sm.set_array([])
# #     cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
# #     cbar.set_label('Edge Weight (Confidence)', color='white', fontsize=10)
# #     cbar.ax.yaxis.set_tick_params(color='white')
# #     plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

# #     ax.set_title(
# #         f"Directed Shortest Path Sub-Network\n"
# #         f"{source}  ->  {target}  |  "
# #         f"{len(paths)} path(s)  |  Total Path Score = {total_score:.4f}",
# #         color='white', fontsize=13, fontweight='bold', pad=20
# #     )

# #     ax.axis('off')
# #     plt.tight_layout()
# #     plt.savefig(output_image, dpi=150, bbox_inches='tight',
# #                 facecolor=fig.get_facecolor())
# #     plt.show()
# #     print(f"[INFO] Figure saved -> {output_image}")
