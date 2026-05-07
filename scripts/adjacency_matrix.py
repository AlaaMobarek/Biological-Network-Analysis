"""
adjacency_matrix.py
───────────────────
Member 4 — Adjacency Matrix

Responsibilities:
    - Convert the directed weighted graph to an unweighted graph
    - Save it using the adjacency matrix method as a CSV
    - Optionally visualize the matrix as a heatmap (for a subgraph sample)

Usage (standalone):
    python scripts/adjacency_matrix.py

Usage (from main.py):
    from scripts.adjacency_matrix import save_adjacency_matrix
    save_adjacency_matrix(G)
"""

import os
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────
MATRIX_DIR = "Code/outputs/matrices"
FIG_DIR    = "Code/outputs/figures"
os.makedirs(MATRIX_DIR, exist_ok=True)
os.makedirs(FIG_DIR,    exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Convert to unweighted graph and save adjacency matrix
# ══════════════════════════════════════════════════════════════════════════════
def save_adjacency_matrix(G: nx.DiGraph, sample_size: int = 50):
    """
    Converts the directed weighted interactome to an unweighted directed graph,
    then saves it as an adjacency matrix CSV file.

    Because the full matrix (~18k × 18k) is too large to save as CSV,
    this function saves two versions:
        1. A sample adjacency matrix (top N hub proteins) as CSV + heatmap
        2. The full adjacency matrix as a sparse edge list

    Parameters
    ----------
    G           : nx.DiGraph — the full loaded interactome graph
    sample_size : int        — number of top hub proteins to include
                               in the sample matrix (default: 50)
    """
    print(f"\n[adjacency_matrix] Converting graph to unweighted...")

    # ── Step 1: Convert to unweighted DiGraph ─────────────────────────────────
    # Copy graph structure but remove weights (all edges become 1)
    G_unweighted = nx.DiGraph()
    for u, v in G.edges():
        G_unweighted.add_edge(u, v)   # no weight attribute = unweighted

    print(f"  Unweighted graph: {G_unweighted.number_of_nodes():,} nodes, "
          f"{G_unweighted.number_of_edges():,} edges")

    # ── Step 2: Save full graph as sparse edge list ────────────────────────────
    # The full adjacency matrix is ~18k×18k = 324M cells — too large for CSV.
    # The standard alternative is a sparse edge-list representation.
    edge_list_path = os.path.join(MATRIX_DIR, "adjacency_edgelist.txt")
    with open(edge_list_path, "w") as f:
        f.write("# Unweighted adjacency edge list\n")
        f.write("# source\ttarget\n")
        for u, v in G_unweighted.edges():
            f.write(f"{u}\t{v}\n")
    print(f"  Full edge list saved → {edge_list_path}")

    # ── Step 3: Build a SAMPLE adjacency matrix (top hub proteins) ────────────
    print(f"\n  Building sample adjacency matrix "
          f"(top {sample_size} hub proteins)...")

    # Pick top N nodes by total degree
    top_nodes = [n for n, _ in
                 sorted(G_unweighted.degree(),
                        key=lambda x: x[1], reverse=True)[:sample_size]]

    # Extract subgraph for these nodes only
    sub = G_unweighted.subgraph(top_nodes).copy()

    # Build adjacency matrix using NetworkX
    # nodelist ensures rows/columns are in the same consistent order
    adj_matrix = nx.to_numpy_array(sub, nodelist=top_nodes, dtype=int)

    # ── Step 4: Save sample matrix as CSV ─────────────────────────────────────
    df = pd.DataFrame(adj_matrix,
                      index=top_nodes,
                      columns=top_nodes)

    csv_path = os.path.join(MATRIX_DIR, "adjacency_matrix.csv")
    df.to_csv(csv_path)
    print(f"  Sample CSV saved  → {csv_path}")
    print(f"  Matrix shape      : {df.shape[0]} × {df.shape[1]}")
    print(f"  Non-zero entries  : {int(adj_matrix.sum())} "
          f"(interactions present)")

    # ── Step 5: Visualize matrix as heatmap ───────────────────────────────────
    _plot_adjacency_heatmap(adj_matrix, top_nodes, sample_size)

    return df


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Plot adjacency matrix heatmap
# ══════════════════════════════════════════════════════════════════════════════
def _plot_adjacency_heatmap(matrix: np.ndarray,
                             labels: list, sample_size: int):
    """
    Plots and saves a heatmap of the sample adjacency matrix.

    Parameters
    ----------
    matrix      : np.ndarray — the adjacency matrix (0s and 1s)
    labels      : list       — protein UniProt IDs (row/column labels)
    sample_size : int        — used in title
    """
    fig, ax = plt.subplots(figsize=(14, 12))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    # Plot heatmap: cyan for 1 (edge exists), dark for 0 (no edge)
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="auto",
                   interpolation="nearest")

    # Axis labels
    fontsize = max(4, int(80 / sample_size))   # scale font with matrix size
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=fontsize,
                        color="#e0e0e0")
    ax.set_yticklabels(labels, fontsize=fontsize, color="#e0e0e0")

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Edge present (1) / absent (0)",
                   color="#e0e0e0", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="#e0e0e0", labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#e0e0e0")

    # Gridlines for readability
    ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="#2e3250", linewidth=0.4)

    # Labels & title
    ax.set_xlabel("Target Protein (head)", color="#e0e0e0", fontsize=10)
    ax.set_ylabel("Source Protein (tail)", color="#e0e0e0", fontsize=10)
    ax.set_title(
        f"Unweighted Adjacency Matrix — Top {sample_size} Hub Proteins\n"
        f"(1 = interaction exists, 0 = no direct interaction)",
        color="#e0e0e0", fontsize=12, fontweight="bold", pad=12
    )

    plt.tight_layout()

    out_path = os.path.join(FIG_DIR, "adjacency_matrix_heatmap.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Heatmap saved     → {out_path}\n")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN — runs when executed directly
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    from load_graph import load_graph

    DATA_PATH = "Code/data/PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    G = load_graph(DATA_PATH)

    save_adjacency_matrix(G, sample_size=50)

    print("✅ adjacency_matrix.py complete!\n")