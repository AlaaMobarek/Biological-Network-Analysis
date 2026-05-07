"""
shortest_path.py
────────────────
Member 2 — Shortest Path Analysis

Responsibilities:
    - Given two proteins, find all acyclic shortest paths between them
    - Report total path score and each edge weight
    - Draw and save the sub-network formed by these paths

Usage (standalone):
    python scripts/shortest_path.py

Usage (from main.py):
    from scripts.shortest_path import find_shortest_paths
    find_shortest_paths(G, source="P04637", target="P00533")
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
# FUNCTION 1 — Find all acyclic shortest paths between two proteins
# ══════════════════════════════════════════════════════════════════════════════
def find_shortest_paths(G: nx.DiGraph, source: str, target: str):
    """
    Finds all acyclic shortest paths between two proteins in the interactome,
    reports path scores and edge weights, saves results to a text file,
    and draws the sub-network of these paths.

    Parameters
    ----------
    G      : nx.DiGraph  — the loaded interactome graph
    source : str         — UniProt ID of the source/start protein
    target : str         — UniProt ID of the target/end protein
    """
    print(f"\n[shortest_path] Finding paths: {source} → {target}")

    # ── Validate nodes exist in graph ─────────────────────────────────────────
    for node, label in [(source, "Source"), (target, "Target")]:
        if node not in G:
            print(f"  ⚠ {label} protein '{node}' not found in graph.")
            return

    # ── Find all simple (acyclic) shortest paths ───────────────────────────────
    try:
        # Get the shortest path length first
        shortest_len = nx.shortest_path_length(G, source=source, target=target, weight= "weight")

        # Get ALL simple paths of exactly that length (all shortest paths)
        all_paths = [
            p for p in nx.all_simple_paths(G, source=source, target=target,
                                           cutoff=shortest_len)
            if len(p) - 1 == shortest_len
        ]
    except nx.NetworkXNoPath:
        print(f"  ✗ No path exists between {source} and {target}.")
        return
    except nx.NodeNotFound as e:
        print(f"  ✗ Node error: {e}")
        return

    print(f"  Found {len(all_paths)} shortest path(s) of length {shortest_len}")

    # ── Compute path scores & edge weights ────────────────────────────────────
    path_results = []

    for path in all_paths:
        edges       = list(zip(path[:-1], path[1:]))           # pairs of nodes
        edge_weights = [G[u][v]["weight"] for u, v in edges]   # weight per edge
        methods      = [G[u][v]["method"] for u, v in edges]   # method per edge

        # Path score = product of all edge weights (combined confidence)
        path_score = float(np.prod(edge_weights))

        path_results.append({
            "path"        : path,
            "edges"       : edges,
            "weights"     : edge_weights,
            "methods"     : methods,
            "score"       : path_score,
        })

    # Sort paths by score descending (best path first)
    path_results.sort(key=lambda x: x["score"], reverse=True)

    # ── Save results to text file ─────────────────────────────────────────────
    out_txt = os.path.join(TXT_DIR, "shortest_paths.txt")
    with open(out_txt, "w") as f:
        f.write(f"Shortest Paths Analysis\n")
        f.write(f"{'='*60}\n")
        f.write(f"Source protein : {source}\n")
        f.write(f"Target protein : {target}\n")
        f.write(f"Path length    : {shortest_len} interaction(s)\n")
        f.write(f"Total paths    : {len(all_paths)}\n")
        f.write(f"{'='*60}\n\n")

        for i, result in enumerate(path_results, 1):
            f.write(f"Path {i}:\n")
            f.write(f"  Route      : {' → '.join(result['path'])}\n")
            f.write(f"  Path score : {result['score']:.6f}  "
                    f"(product of all edge weights)\n")
            f.write(f"  Edge details:\n")
            for (u, v), w, m in zip(result["edges"],
                                     result["weights"],
                                     result["methods"]):
                f.write(f"    {u} → {v}  |  weight: {w:.6f}  |  "
                        f"method: {m}\n")
            f.write("\n")

    print(f"  Results saved → {out_txt}")

    # ── Draw subnetwork ────────────────────────────────────────────────────────
    _draw_path_subnetwork(G, path_results, source, target)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Draw the sub-network of shortest paths
# ══════════════════════════════════════════════════════════════════════════════
def _draw_path_subnetwork(G: nx.DiGraph, path_results: list,
                          source: str, target: str):
    """
    Draws and saves the sub-network formed by all shortest paths.

    Parameters
    ----------
    G            : nx.DiGraph — full graph (for edge attribute lookup)
    path_results : list       — output of find_shortest_paths()
    source       : str        — source protein UniProt ID
    target       : str        — target protein UniProt ID
    """
    # ── Build subgraph from all path edges ────────────────────────────────────
    sub_edges = set()
    sub_nodes = set()
    for result in path_results:
        for edge in result["edges"]:
            sub_edges.add(edge)
        sub_nodes.update(result["path"])

    sub = G.edge_subgraph(sub_edges).copy()

    # ── Layout & styling ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 8))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    # Hierarchical left-to-right layout using shell layout
    pos = nx.spring_layout(sub, seed=7, k=1.2)

    # Node colors: source=green, target=red, intermediate=cyan
    node_colors = []
    for n in sub.nodes():
        if n == source:
            node_colors.append("#2ecc71")   # green
        elif n == target:
            node_colors.append("#e74c3c")   # red
        else:
            node_colors.append("#00d4ff")   # cyan

    node_sizes = [900 if n in (source, target) else 500
                  for n in sub.nodes()]

    # Edge colors by weight
    edge_weights = [sub[u][v]["weight"] for u, v in sub.edges()]
    edge_colors  = [plt.cm.YlOrRd(w) for w in edge_weights]

    # Draw
    nx.draw_networkx_nodes(sub, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.95, ax=ax)
    nx.drrkx_edges(sub, pos, edge_color=edge_colors,
                           width=2.5, arrows=True, arrowsize=18,
                           connectionstyle="arc3,rad=0.12", ax=ax)
    nx.draw_networkx_labels(sub, pos, font_size=8,
                            font_color="white", font_weight="bold", ax=ax)

    # Edge weight labels
    edge_labels = {(u, v): f"{sub[u][v]['weight']:.3f}"
                   for u, v in sub.edges()}
    nx.draw_networkx_edge_labels(sub, pos, edge_labels=edge_labels,
                                 font_size=6, font_color="#f0e68c", ax=ax)

    # Colorbar
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.YlOrRd,
        norm=plt.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Edge Weight (confidence)", color="#e0e0e0", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="#e0e0e0", labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#e0e0e0")

    # Legend
    legend_patches = [
        mpatches.Patch(color="#2ecc71", label=f"Source: {source}"),
        mpatches.Patch(color="#e74c3c", label=f"Target: {target}"),
        mpatches.Patch(color="#00d4ff", label="Intermediate protein"),
    ]
    ax.legend(handles=legend_patches, loc="upper left",
              facecolor="#1a1d27", edgecolor="#2e3250",
              labelcolor="#e0e0e0", fontsize=8)

    ax.set_title(
        f"Shortest Path Sub-Network: {source} → {target}\n"
        f"({len(path_results)} path(s), length = "
        f"{len(path_results[0]['path']) - 1} interactions)",
        color="#e0e0e0", fontsize=12, fontweight="bold", pad=10
    )
    ax.axis("off")
    plt.tight_layout()

    out_path = os.path.join(FIG_DIR, "shortest_path_subnetwork.png")
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

    # Change these two proteins to any UniProt IDs in your dataset
    find_shortest_paths(G, source="P04637", target="P00533")

    print("✅ shortest_path.py complete!\n")