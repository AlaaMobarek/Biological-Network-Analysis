"""
neighbors.py
────────────
Member 3 — Neighbor Analysis

Responsibilities:
    - Given one protein, list all directly connected proteins
    - Report the degree of the protein
    - Report each neighbor with its interaction weight and direction
    - Save results to a text file

Usage (standalone):
    python scripts/neighbors.py

Usage (from main.py):
    from scripts.neighbors import get_neighbors
    get_neighbors(G, protein="P04637")
"""

import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────────
TXT_DIR = "Code/outputs/txt"
FIG_DIR = "Code/outputs/figures"
os.makedirs(TXT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Get all neighbors of a given protein
# ══════════════════════════════════════════════════════════════════════════════
def get_neighbors(G: nx.DiGraph, protein: str):
    """
    Lists all proteins directly connected to the given protein (in both
    directions), reports their interaction weights and directions, and
    saves results to a text file.

    Parameters
    ----------
    G       : nx.DiGraph — the loaded interactome graph
    protein : str        — UniProt ID of the query protein
    """
    print(f"\n[neighbors] Analyzing neighbors of: {protein}")

    # ── Validate node exists ───────────────────────────────────────────────────
    if protein not in G:
        print(f"  ⚠ Protein '{protein}' not found in graph.")
        return

    # ── Collect successors (protein → neighbor) ────────────────────────────────
    successors = []
    for neighbor in G.successors(protein):
        w = G[protein][neighbor]["weight"]
        m = G[protein][neighbor]["method"]
        successors.append((neighbor, w, m, "outgoing"))

    # ── Collect predecessors (neighbor → protein) ──────────────────────────────
    predecessors = []
    for neighbor in G.predecessors(protein):
        w = G[neighbor][protein]["weight"]
        m = G[neighbor][protein]["method"]
        predecessors.append((neighbor, w, m, "incoming"))

    # ── Degree stats ───────────────────────────────────────────────────────────
    total_degree  = G.degree(protein)
    in_degree     = G.in_degree(protein)
    out_degree    = G.out_degree(protein)

    # Sort each list by weight descending
    successors   = sorted(successors,   key=lambda x: x[1], reverse=True)
    predecessors = sorted(predecessors, key=lambda x: x[1], reverse=True)

    all_neighbors = successors + predecessors

    print(f"  Total degree  : {total_degree}  "
          f"(in={in_degree}, out={out_degree})")
    print(f"  Outgoing neighbors : {len(successors)}")
    print(f"  Incoming neighbors : {len(predecessors)}")

    # ── Save to text file ─────────────────────────────────────────────────────
    out_txt = os.path.join(TXT_DIR, "neighbors.txt")
    with open(out_txt, "w") as f:
        f.write(f"Neighbor Analysis\n")
        f.write(f"{'='*60}\n")
        f.write(f"Query protein  : {protein}\n")
        f.write(f"Total degree   : {total_degree}\n")
        f.write(f"  In-degree    : {in_degree}   "
                f"(proteins that interact WITH {protein})\n")
        f.write(f"  Out-degree   : {out_degree}  "
                f"(proteins that {protein} interacts WITH)\n")
        f.write(f"{'='*60}\n\n")

        # Outgoing neighbors
        f.write(f"OUTGOING Neighbors ({len(successors)}) "
                f"— {protein} → neighbor:\n")
        f.write(f"{'─'*60}\n")
        for neighbor, weight, method, _ in successors:
            f.write(f"  {neighbor:<12}  weight: {weight:.6f}  "
                    f"method: {method}\n")

        f.write(f"\nINCOMING Neighbors ({len(predecessors)}) "
                f"— neighbor → {protein}:\n")
        f.write(f"{'─'*60}\n")
        for neighbor, weight, method, _ in predecessors:
            f.write(f"  {neighbor:<12}  weight: {weight:.6f}  "
                    f"method: {method}\n")

    print(f"  Results saved → {out_txt}")

    # ── Visualize ─────────────────────────────────────────────────────────────
    _draw_neighbor_graph(G, protein, successors, predecessors)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Draw the neighbor sub-network
# ══════════════════════════════════════════════════════════════════════════════
def _draw_neighbor_graph(G: nx.DiGraph, protein: str,
                         successors: list, predecessors: list):
    """
    Draws and saves the ego sub-network centered on the query protein,
    color-coding outgoing and incoming neighbors separately.
    """
    # Build ego subgraph (radius=1)
    ego = nx.ego_graph(G, protein, radius=1)

    fig, ax = plt.subplots(figsize=(13, 9))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    pos = nx.spring_layout(ego, seed=42, k=0.8)

    # Identify node types
    succ_set = {n for n, _, _, _ in successors}
    pred_set = {n for n, _, _, _ in predecessors}
    both_set = succ_set & pred_set     # neighbors that are both in & out

    node_colors = []
    node_sizes  = []
    for n in ego.nodes():
        if n == protein:
            node_colors.append("#f39c12")   # orange = query protein
            node_sizes.append(1200)
        elif n in both_set:
            node_colors.append("#9b59b6")   # purple = bidirectional
            node_sizes.append(400)
        elif n in succ_set:
            node_colors.append("#2ecc71")   # green = outgoing
            node_sizes.append(350)
        else:
            node_colors.append("#e74c3c")   # red = incoming
            node_sizes.append(350)

    edge_weights = [ego[u][v]["weight"] for u, v in ego.edges()]
    edge_colors  = [plt.cm.YlOrRd(w) for w in edge_weights]

    nx.draw_networkx_nodes(ego, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.92, ax=ax)
    nx.draw_networkx_edges(ego, pos, edge_color=edge_colors,
                           width=1.5, arrows=True, arrowsize=12,
                           connectionstyle="arc3,rad=0.1", ax=ax)
    nx.draw_networkx_labels(ego, pos, font_size=6,
                            font_color="white", ax=ax)

    # Legend
    legend_patches = [
        mpatches.Patch(color="#f39c12", label=f"Query: {protein}"),
        mpatches.Patch(color="#2ecc71", label="Outgoing (→ neighbor)"),
        mpatches.Patch(color="#e74c3c", label="Incoming (neighbor →)"),
        mpatches.Patch(color="#9b59b6", label="Bidirectional"),
    ]
    ax.legend(handles=legend_patches, loc="upper left",
              facecolor="#1a1d27", edgecolor="#2e3250",
              labelcolor="#e0e0e0", fontsize=8)

    ax.set_title(
        f"Neighbor Network of {protein}\n"
        f"Total degree: {G.degree(protein)}  "
        f"(in={G.in_degree(protein)}, out={G.out_degree(protein)})",
        color="#e0e0e0", fontsize=12, fontweight="bold", pad=10
    )
    ax.axis("off")
    plt.tight_layout()

    out_path = os.path.join(FIG_DIR, "neighbors_subgraph.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Figure saved  → {out_path}\n")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN — runs when executed directly
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    from load_graph import load_graph

    DATA_PATH = "Code/data/PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    G = load_graph(DATA_PATH)

    # Change this to any UniProt ID in your dataset
    get_neighbors(G, protein="P04637")

    print("✅ neighbors.py complete!\n")