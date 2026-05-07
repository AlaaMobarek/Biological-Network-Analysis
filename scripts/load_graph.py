# """
# load_graph.py
# ─────────────
# Member 1 — Graph Construction & Visualization

# Responsibilities:
#     - Load the full PathLinker human PPI interactome into a directed NetworkX graph
#     - Compute and print full network statistics
#     - Visualize a sample ego subgraph + full statistics dashboard
#     - Save the figure to outputs/figures/graph_overview.png

# Usage (standalone):
#     python scripts/load_graph.py

# Usage (from main.py):
#     from scripts.load_graph import load_graph
#     G = load_graph("data/interactome.txt")
# """

# import os
# import networkx as nx
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np


# # ── Constants ──────────────────────────────────────────────────────────────────
# DATA_PATH = "data/interactome.txt"
# FIG_DIR   = "outputs/figures"
# os.makedirs(FIG_DIR, exist_ok=True)


# # ══════════════════════════════════════════════════════════════════════════════
# # FUNCTION 1 — Load the interactome into a NetworkX DiGraph
# # ══════════════════════════════════════════════════════════════════════════════
# def load_graph(path: str) -> nx.DiGraph:
#     """
#     Reads the PathLinker interactome file and builds a directed weighted graph.

#     Parameters
#     ----------
#     path : str
#         Path to the interactome .txt file (tab-separated, 4 columns).

#     Returns
#     -------
#     G : nx.DiGraph
#         Directed graph where:
#             - nodes  = proteins (UniProt IDs)
#             - edges  = interactions with 'weight' and 'method' attributes
#     """
#     print(f"\n[load_graph] Reading file: {path}")

#     G = nx.DiGraph()  # Directed graph: tail → head

#     with open(path, "r") as f:
#         for line in f:
#             # Skip header / comment lines
#             if line.startswith("#") or line.strip() == "":
#                 continue

#             parts = line.strip().split("\t")

#             # Need at least 3 columns: tail, head, weight
#             if len(parts) < 3:
#                 continue

#             tail   = parts[0]                               # Source protein
#             head   = parts[1]                               # Target protein
#             weight = float(parts[2])                        # Confidence [0,1]
#             method = parts[3] if len(parts) > 3 else "unknown"  # Method(s)

#             # Add directed edge with attributes stored on the edge
#             G.add_edge(tail, head, weight=weight, method=method)

#     print(f"[load_graph] Graph loaded successfully.\n")
#     return G


# # ══════════════════════════════════════════════════════════════════════════════
# # FUNCTION 2 — Compute and print full network statistics
# # ══════════════════════════════════════════════════════════════════════════════
# def compute_statistics(G: nx.DiGraph) -> dict:
#     """
#     Computes comprehensive statistics for the full interactome graph.

#     Parameters
#     ----------
#     G : nx.DiGraph
#         The loaded directed graph.

#     Returns
#     -------
#     stats : dict
#         Dictionary containing all computed statistics.
#     """
#     print("[compute_statistics] Computing network statistics...")

#     # ── Basic counts ──────────────────────────────────────────────────────────
#     num_nodes = G.number_of_nodes()
#     num_edges = G.number_of_edges()

#     # ── Degree statistics ─────────────────────────────────────────────────────
#     degrees     = [d for _, d in G.degree()]       # total degree (in + out)
#     in_degrees  = [d for _, d in G.in_degree()]    # in-degree only
#     out_degrees = [d for _, d in G.out_degree()]   # out-degree only

#     avg_degree    = np.mean(degrees)
#     max_degree    = max(degrees)
#     min_degree    = min(degrees)
#     median_degree = np.median(degrees)

#     # ── Edge weight statistics ────────────────────────────────────────────────
#     weights    = [data["weight"] for _, _, data in G.edges(data=True)]
#     avg_weight = np.mean(weights)
#     max_weight = max(weights)
#     min_weight = min(weights)

#     # ── Connectivity ──────────────────────────────────────────────────────────
#     weakly_connected = nx.is_weakly_connected(G)
#     num_wcc          = nx.number_weakly_connected_components(G)
#     largest_wcc_size = len(max(nx.weakly_connected_components(G), key=len))

#     # ── Density ───────────────────────────────────────────────────────────────
#     density = nx.density(G)

#     # ── Top 10 hub proteins by total degree ───────────────────────────────────
#     top_proteins = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]

#     # ── Bundle everything into a dictionary ───────────────────────────────────
#     stats = {
#         "num_nodes"        : num_nodes,
#         "num_edges"        : num_edges,
#         "avg_degree"       : avg_degree,
#         "max_degree"       : max_degree,
#         "min_degree"       : min_degree,
#         "median_degree"    : median_degree,
#         "avg_weight"       : avg_weight,
#         "max_weight"       : max_weight,
#         "min_weight"       : min_weight,
#         "weakly_connected" : weakly_connected,
#         "num_wcc"          : num_wcc,
#         "largest_wcc_size" : largest_wcc_size,
#         "density"          : density,
#         "top_proteins"     : top_proteins,
#         "degrees"          : degrees,
#         "in_degrees"       : in_degrees,
#         "out_degrees"      : out_degrees,
#         "weights"          : weights,
#     }

#     # ── Pretty-print to console ───────────────────────────────────────────────
#     print(f"\n{'═'*52}")
#     print(f"   FULL NETWORK STATISTICS")
#     print(f"{'═'*52}")
#     print(f"  Nodes (proteins)            : {num_nodes:,}")
#     print(f"  Edges (interactions)        : {num_edges:,}")
#     print(f"  Graph density               : {density:.2e}")
#     print(f"  Is directed                 : Yes")
#     print(f"  Weakly connected            : {weakly_connected}")
#     print(f"  Weakly connected components : {num_wcc:,}")
#     print(f"  Largest component size      : {largest_wcc_size:,}")
#     print(f"{'─'*52}")
#     print(f"  Degree Statistics:")
#     print(f"    Average degree            : {avg_degree:.2f}")
#     print(f"    Median  degree            : {median_degree:.0f}")
#     print(f"    Max     degree            : {max_degree}")
#     print(f"    Min     degree            : {min_degree}")
#     print(f"{'─'*52}")
#     print(f"  Edge Weight Statistics:")
#     print(f"    Average weight            : {avg_weight:.4f}")
#     print(f"    Max     weight            : {max_weight:.4f}")
#     print(f"    Min     weight            : {min_weight:.4f}")
#     print(f"{'─'*52}")
#     print(f"  Top 5 Hub Proteins (by degree):")
#     for i, (node, deg) in enumerate(top_proteins[:5], 1):
#         print(f"    {i}. {node}  →  degree {deg}")
#     print(f"{'═'*52}\n")

#     return stats


# # ══════════════════════════════════════════════════════════════════════════════
# # FUNCTION 3 — Plot full network statistics dashboard + ego subgraph
# # ══════════════════════════════════════════════════════════════════════════════
# def plot_statistics(G: nx.DiGraph, stats: dict):
#     """
#     Creates a 7-panel statistics dashboard for the full interactome:
#         Panel 1 : Total degree distribution (log scale)
#         Panel 2 : In-degree distribution
#         Panel 3 : Out-degree distribution
#         Panel 4 : Edge weight distribution
#         Panel 5 : Top 10 hub proteins (bar chart)
#         Panel 6 : Network summary text box
#         Panel 7 : Ego subgraph of the top hub protein (full bottom row)

#     Saves figure to: outputs/figures/graph_overview.png

#     Parameters
#     ----------
#     G     : nx.DiGraph  — the loaded graph
#     stats : dict        — output of compute_statistics()
#     """
#     print("[plot_statistics] Generating statistics dashboard...")

#     # ── Figure & dark theme setup ─────────────────────────────────────────────
#     fig = plt.figure(figsize=(20, 14))
#     fig.patch.set_facecolor("#0f1117")

#     # 3-row, 3-column grid; bottom row spans all 3 columns for ego graph
#     gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.35)

#     # Color palette
#     ACCENT   = "#00d4ff"   # cyan
#     ACCENT2  = "#ff6b6b"   # coral/red  (used for seed node)
#     ACCENT3  = "#a8ff78"   # lime green
#     BG_PANEL = "#1a1d27"   # panel background
#     TEXT_COL = "#e0e0e0"   # axis text

#     def style_ax(ax, title):
#         """Apply consistent dark-theme styling to an axes panel."""
#         ax.set_facecolor(BG_PANEL)
#         ax.set_title(title, color=TEXT_COL, fontsize=11,
#                      fontweight="bold", pad=8)
#         ax.tick_params(colors=TEXT_COL, labelsize=8)
#         for spine in ax.spines.values():
#             spine.set_edgecolor("#2e3250")

#     # ── Panel 1 : Total degree distribution ───────────────────────────────────
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.hist(stats["degrees"], bins=60, color=ACCENT,
#              edgecolor="none", alpha=0.85, log=True)
#     ax1.set_xlabel("Degree", color=TEXT_COL, fontsize=8)
#     ax1.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
#     ax1.axvline(stats["avg_degree"], color=ACCENT2, linewidth=1.5,
#                 linestyle="--", label=f"Mean = {stats['avg_degree']:.1f}")
#     ax1.legend(fontsize=7, labelcolor=TEXT_COL,
#                facecolor=BG_PANEL, edgecolor="none")
#     style_ax(ax1, "Total Degree Distribution")

#     # ── Panel 2 : In-degree distribution ──────────────────────────────────────
#     ax2 = fig.add_subplot(gs[0, 1])
#     ax2.hist(stats["in_degrees"], bins=60, color=ACCENT3,
#              edgecolor="none", alpha=0.85, log=True)
#     ax2.set_xlabel("In-Degree", color=TEXT_COL, fontsize=8)
#     ax2.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
#     style_ax(ax2, "In-Degree Distribution")

#     # ── Panel 3 : Out-degree distribution ─────────────────────────────────────
#     ax3 = fig.add_subplot(gs[0, 2])
#     ax3.hist(stats["out_degrees"], bins=60, color=ACCENT2,
#              edgecolor="none", alpha=0.85, log=True)
#     ax3.set_xlabel("Out-Degree", color=TEXT_COL, fontsize=8)
#     ax3.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
#     style_ax(ax3, "Out-Degree Distribution")

#     # ── Panel 4 : Edge weight distribution ────────────────────────────────────
#     ax4 = fig.add_subplot(gs[1, 0])
#     ax4.hist(stats["weights"], bins=50, color="#f7971e",
#              edgecolor="none", alpha=0.85)
#     ax4.set_xlabel("Edge Weight (confidence score)", color=TEXT_COL, fontsize=8)
#     ax4.set_ylabel("Count", color=TEXT_COL, fontsize=8)
#     ax4.axvline(stats["avg_weight"], color=ACCENT, linewidth=1.5,
#                 linestyle="--", label=f"Mean = {stats['avg_weight']:.3f}")
#     ax4.legend(fontsize=7, labelcolor=TEXT_COL,
#                facecolor=BG_PANEL, edgecolor="none")
#     style_ax(ax4, "Edge Weight Distribution")

#     # ── Panel 5 : Top 10 hub proteins horizontal bar chart ────────────────────
#     ax5 = fig.add_subplot(gs[1, 1])
#     top10  = stats["top_proteins"]                      # already sorted top-10
#     labels = [n for n, _ in top10]
#     values = [d for _, d in top10]
#     colors = plt.cm.plasma(np.linspace(0.3, 0.9, len(labels)))
#     bars   = ax5.barh(labels[::-1], values[::-1],
#                       color=colors, edgecolor="none")
#     ax5.set_xlabel("Total Degree", color=TEXT_COL, fontsize=8)
#     for bar, val in zip(bars, values[::-1]):
#         ax5.text(bar.get_width() + 0.3,
#                  bar.get_y() + bar.get_height() / 2,
#                  str(val), va="center", color=TEXT_COL, fontsize=7)
#     style_ax(ax5, "Top 10 Hub Proteins")

#     # ── Panel 6 : Summary text box ────────────────────────────────────────────
#     ax6 = fig.add_subplot(gs[1, 2])
#     ax6.set_facecolor(BG_PANEL)
#     ax6.axis("off")
#     summary_text = (
#         f"NETWORK SUMMARY\n"
#         f"{'─'*30}\n"
#         f"Proteins (nodes)   : {stats['num_nodes']:,}\n"
#         f"Interactions       : {stats['num_edges']:,}\n"
#         f"Graph density      : {stats['density']:.2e}\n"
#         f"Avg degree         : {stats['avg_degree']:.2f}\n"
#         f"Median degree      : {stats['median_degree']:.0f}\n"
#         f"Max degree         : {stats['max_degree']}\n"
#         f"Avg edge weight    : {stats['avg_weight']:.4f}\n"
#         f"Max edge weight    : {stats['max_weight']:.4f}\n"
#         f"Min edge weight    : {stats['min_weight']:.4f}\n"
#         f"Connected comps    : {stats['num_wcc']:,}\n"
#         f"Largest comp size  : {stats['largest_wcc_size']:,}\n"
#         f"Directed graph     : Yes\n"
#     )
#     ax6.text(0.05, 0.97, summary_text,
#              transform=ax6.transAxes,
#              fontsize=9, verticalalignment="top",
#              fontfamily="monospace", color=ACCENT,
#              bbox=dict(facecolor="#0f1117", edgecolor=ACCENT,
#                        boxstyle="round,pad=0.6", linewidth=1))
#     ax6.set_title("Statistics Summary", color=TEXT_COL,
#                   fontsize=11, fontweight="bold", pad=8)

#     # ── Panel 7 : Ego subgraph — full bottom row ───────────────────────────────
#     ax7 = fig.add_subplot(gs[2, :])   # spans all 3 columns
#     ax7.set_facecolor(BG_PANEL)

#     # Use the highest-degree node as the ego seed
#     seed = stats["top_proteins"][0][0]
#     ego  = nx.ego_graph(G, seed, radius=1)

#     pos = nx.spring_layout(ego, seed=42, k=0.65)

#     # Color: seed = coral, neighbors = cyan
#     node_colors = [ACCENT2 if n == seed else ACCENT for n in ego.nodes()]
#     node_sizes  = [900    if n == seed else 200    for n in ego.nodes()]

#     # Edge color mapped to weight
#     edge_weights = [ego[u][v]["weight"] for u, v in ego.edges()]
#     edge_colors  = [plt.cm.YlOrRd(w) for w in edge_weights]

#     nx.draw_networkx_nodes(ego, pos, node_color=node_colors,
#                            node_size=node_sizes, alpha=0.9, ax=ax7)
#     nx.draw_networkx_edges(ego, pos, edge_color=edge_colors,
#                            width=1.2, arrows=True, arrowsize=10,
#                            connectionstyle="arc3,rad=0.08", ax=ax7)
#     nx.draw_networkx_labels(ego, pos, font_size=5.5,
#                             font_color="white", ax=ax7)

#     # Colorbar for edge weights
#     sm = plt.cm.ScalarMappable(
#         cmap=plt.cm.YlOrRd,
#         norm=plt.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
#     )
#     sm.set_array([])
#     cbar = plt.colorbar(sm, ax=ax7, fraction=0.015, pad=0.01)
#     cbar.set_label("Edge Weight (confidence)", color=TEXT_COL, fontsize=8)
#     cbar.ax.yaxis.set_tick_params(color=TEXT_COL, labelsize=7)
#     plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COL)

#     ax7.set_title(
#         f"Ego Subgraph of Top Hub Protein: {seed}  "
#         f"({ego.number_of_nodes()} proteins, {ego.number_of_edges()} interactions)",
#         color=TEXT_COL, fontsize=11, fontweight="bold", pad=8
#     )
#     ax7.axis("off")

#     # ── Main figure title ─────────────────────────────────────────────────────
#     fig.suptitle(
#         "PathLinker Human PPI Interactome — Full Network Analysis Dashboard",
#         fontsize=15, fontweight="bold", color=TEXT_COL, y=0.99
#     )

#     # ── Save figure ───────────────────────────────────────────────────────────
#     out_path = os.path.join(FIG_DIR, "graph_overview.png")
#     plt.savefig(out_path, dpi=150, bbox_inches="tight",
#                 facecolor=fig.get_facecolor())
#     plt.close()
#     print(f"[plot_statistics] Figure saved → {out_path}\n")


# # ══════════════════════════════════════════════════════════════════════════════
# # MAIN — runs when executed directly: python scripts/load_graph.py
# # ══════════════════════════════════════════════════════════════════════════════
# if __name__ == "__main__":
#     G     = load_graph(DATA_PATH)       # Step 1: load graph
#     stats = compute_statistics(G)       # Step 2: compute statistics
#     plot_statistics(G, stats)           # Step 3: plot dashboard + ego subgraph
#     print("✅ load_graph.py complete!\n")
"""
load_graph.py
─────────────
Member 1 — Graph Construction & Visualization

Saves TWO separate figures:
    1. graph_statistics.png  — 6-panel statistics dashboard (no network)
    2. graph_network.png     — clean ego subgraph of the top hub protein
                               (sampled to max 80 neighbors for readability)

Usage (standalone):
    python scripts/load_graph.py

Usage (from main.py):
    from scripts.load_graph import load_graph, compute_statistics, plot_statistics
    G     = load_graph("data/interactome.txt")
    stats = compute_statistics(G)
    plot_statistics(G, stats)
"""

import os
import random
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


# ── Constants ──────────────────────────────────────────────────────────────────
DATA_PATH = "data/interactome.txt"
FIG_DIR   = "outputs/figures"
os.makedirs(FIG_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Load the interactome
# ══════════════════════════════════════════════════════════════════════════════
def load_graph(path: str) -> nx.DiGraph:
    """
    Reads the PathLinker interactome and builds a directed weighted graph.

    Parameters
    ----------
    path : str — path to the tab-separated interactome file

    Returns
    -------
    G : nx.DiGraph
        Nodes  = proteins (UniProt IDs)
        Edges  = interactions with 'weight' (confidence) and 'method'
    """
    print(f"\n[load_graph] Reading: {path}")
    G = nx.DiGraph()

    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            tail   = parts[0]
            head   = parts[1]
            weight = float(parts[2])
            method = parts[3] if len(parts) > 3 else "unknown"
            G.add_edge(tail, head, weight=weight, method=method)

    print(f"[load_graph] Done.\n")
    return G


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Compute statistics
# ══════════════════════════════════════════════════════════════════════════════
def compute_statistics(G: nx.DiGraph) -> dict:
    """
    Computes comprehensive statistics for the full interactome.

    Parameters
    ----------
    G : nx.DiGraph

    Returns
    -------
    stats : dict of all computed values
    """
    print("[compute_statistics] Computing...")

    num_nodes   = G.number_of_nodes()
    num_edges   = G.number_of_edges()
    degrees     = [d for _, d in G.degree()]
    in_degrees  = [d for _, d in G.in_degree()]
    out_degrees = [d for _, d in G.out_degree()]
    weights     = [d["weight"] for _, _, d in G.edges(data=True)]

    avg_degree    = np.mean(degrees)
    median_degree = np.median(degrees)
    max_degree    = max(degrees)
    min_degree    = min(degrees)

    avg_weight = np.mean(weights)
    max_weight = max(weights)
    min_weight = min(weights)

    num_wcc          = nx.number_weakly_connected_components(G)
    largest_wcc_size = len(max(nx.weakly_connected_components(G), key=len))
    density          = nx.density(G)
    top_proteins     = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]

    stats = {
        "num_nodes": num_nodes, "num_edges": num_edges,
        "degrees": degrees, "in_degrees": in_degrees,
        "out_degrees": out_degrees, "weights": weights,
        "avg_degree": avg_degree, "median_degree": median_degree,
        "max_degree": max_degree, "min_degree": min_degree,
        "avg_weight": avg_weight, "max_weight": max_weight,
        "min_weight": min_weight, "num_wcc": num_wcc,
        "largest_wcc_size": largest_wcc_size, "density": density,
        "top_proteins": top_proteins,
    }

    print(f"\n{'═'*52}")
    print(f"   FULL NETWORK STATISTICS")
    print(f"{'═'*52}")
    print(f"  Nodes              : {num_nodes:,}")
    print(f"  Edges              : {num_edges:,}")
    print(f"  Density            : {density:.2e}")
    print(f"  Avg degree         : {avg_degree:.2f}")
    print(f"  Median degree      : {median_degree:.0f}")
    print(f"  Max degree         : {max_degree}")
    print(f"  Avg weight         : {avg_weight:.4f}")
    print(f"  WCC components     : {num_wcc:,}")
    print(f"  Largest comp size  : {largest_wcc_size:,}")
    print(f"{'─'*52}")
    print(f"  Top 5 hub proteins:")
    for i, (n, d) in enumerate(top_proteins[:5], 1):
        print(f"    {i}. {n}  →  degree {d}")
    print(f"{'═'*52}\n")

    return stats


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Plot statistics dashboard (NO network)
# ══════════════════════════════════════════════════════════════════════════════
def _plot_stats_dashboard(stats: dict):
    """
    Saves a clean 6-panel statistics dashboard to:
        outputs/figures/graph_statistics.png

    Panels:
        1 — Total degree distribution (log)
        2 — In-degree distribution (log)
        3 — Out-degree distribution (log)
        4 — Edge weight distribution
        5 — Top 10 hub proteins bar chart
        6 — Summary text box
    """
    print("[plot_statistics] Saving statistics dashboard...")

    fig = plt.figure(figsize=(18, 10))
    fig.patch.set_facecolor("#0f1117")
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)

    BG   = "#1a1d27"
    TXT  = "#e0e0e0"
    CYAN = "#00d4ff"
    RED  = "#ff6b6b"
    GRN  = "#a8ff78"

    def style(ax, title):
        ax.set_facecolor(BG)
        ax.set_title(title, color=TXT, fontsize=11, fontweight="bold", pad=8)
        ax.tick_params(colors=TXT, labelsize=8)
        ax.xaxis.label.set_color(TXT)
        ax.yaxis.label.set_color(TXT)
        for sp in ax.spines.values():
            sp.set_edgecolor("#2e3250")

    # ── Panel 1: Total degree ─────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(stats["degrees"], bins=60, color=CYAN,
             edgecolor="none", alpha=0.85, log=True)
    ax1.set_xlabel("Degree")
    ax1.set_ylabel("Count (log scale)")
    ax1.axvline(stats["avg_degree"], color=RED, lw=1.5, ls="--",
                label=f"Mean = {stats['avg_degree']:.1f}")
    ax1.legend(fontsize=7, labelcolor=TXT, facecolor=BG, edgecolor="none")
    style(ax1, "Total Degree Distribution")

    # ── Panel 2: In-degree ────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(stats["in_degrees"], bins=60, color=GRN,
             edgecolor="none", alpha=0.85, log=True)
    ax2.set_xlabel("In-Degree")
    ax2.set_ylabel("Count (log scale)")
    style(ax2, "In-Degree Distribution")

    # ── Panel 3: Out-degree ───────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(stats["out_degrees"], bins=60, color=RED,
             edgecolor="none", alpha=0.85, log=True)
    ax3.set_xlabel("Out-Degree")
    ax3.set_ylabel("Count (log scale)")
    style(ax3, "Out-Degree Distribution")

    # ── Panel 4: Edge weights ─────────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.hist(stats["weights"], bins=50, color="#f7971e",
             edgecolor="none", alpha=0.85)
    ax4.set_xlabel("Edge Weight (confidence score)")
    ax4.set_ylabel("Count")
    ax4.axvline(stats["avg_weight"], color=CYAN, lw=1.5, ls="--",
                label=f"Mean = {stats['avg_weight']:.3f}")
    ax4.legend(fontsize=7, labelcolor=TXT, facecolor=BG, edgecolor="none")
    style(ax4, "Edge Weight Distribution")

    # ── Panel 5: Top 10 hub proteins ──────────────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    labels = [n for n, _ in stats["top_proteins"]]
    values = [d for _, d in stats["top_proteins"]]
    colors = plt.cm.plasma(np.linspace(0.3, 0.9, len(labels)))
    bars   = ax5.barh(labels[::-1], values[::-1],
                      color=colors, edgecolor="none")
    for bar, val in zip(bars, values[::-1]):
        ax5.text(bar.get_width() + 0.5,
                 bar.get_y() + bar.get_height() / 2,
                 str(val), va="center", color=TXT, fontsize=7)
    ax5.set_xlabel("Total Degree")
    style(ax5, "Top 10 Hub Proteins")

    # ── Panel 6: Summary text ─────────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_facecolor(BG)
    ax6.axis("off")
    summary = (
        f"NETWORK SUMMARY\n"
        f"{'─'*30}\n"
        f"Proteins (nodes)   : {stats['num_nodes']:,}\n"
        f"Interactions       : {stats['num_edges']:,}\n"
        f"Graph density      : {stats['density']:.2e}\n"
        f"Avg degree         : {stats['avg_degree']:.2f}\n"
        f"Median degree      : {stats['median_degree']:.0f}\n"
        f"Max degree         : {stats['max_degree']}\n"
        f"Avg edge weight    : {stats['avg_weight']:.4f}\n"
        f"Max edge weight    : {stats['max_weight']:.4f}\n"
        f"Min edge weight    : {stats['min_weight']:.4f}\n"
        f"Connected comps    : {stats['num_wcc']:,}\n"
        f"Largest comp size  : {stats['largest_wcc_size']:,}\n"
        f"Directed graph     : Yes\n"
    )
    ax6.text(0.05, 0.97, summary, transform=ax6.transAxes,
             fontsize=9.5, va="top", fontfamily="monospace",
             color=CYAN,
             bbox=dict(facecolor="#0f1117", edgecolor=CYAN,
                       boxstyle="round,pad=0.6", lw=1))
    ax6.set_title("Statistics Summary", color=TXT,
                  fontsize=11, fontweight="bold", pad=8)

    fig.suptitle(
        "PathLinker Human PPI Interactome — Network Statistics",
        fontsize=14, fontweight="bold", color=TXT, y=1.01
    )

    out = os.path.join(FIG_DIR, "graph_statistics.png")
    plt.savefig(out, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Statistics figure saved → {out}")


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 4 — Plot ego subgraph (separate figure, sampled for readability)
# ══════════════════════════════════════════════════════════════════════════════
# def _plot_ego_network(G: nx.DiGraph, stats: dict, max_neighbors: int = 60):
#     """
#     Saves a clean ego subgraph of the top hub protein to:
#         outputs/figures/graph_network.png

#     To keep the figure readable, only up to `max_neighbors` neighbors
#     are shown (highest-confidence edges kept, rest trimmed).

#     Parameters
#     ----------
#     G             : nx.DiGraph
#     stats         : dict — output of compute_statistics()
#     max_neighbors : int  — max number of neighbors to display (default 60)
#     """
#     print(f"[plot_network] Saving ego subgraph "
#           f"(max {max_neighbors} neighbors)...")

#     seed = stats["top_proteins"][0][0]   # top hub protein

#     # ── Sample neighbors: keep only top N by edge confidence ──────────────────
#     # Collect all successors + predecessors with their confidence
#     neighbor_edges = []
#     for nbr in G.successors(seed):
#         neighbor_edges.append((seed, nbr, G[seed][nbr]["weight"]))
#     for nbr in G.predecessors(seed):
#         neighbor_edges.append((nbr, seed, G[nbr][seed]["weight"]))

#     # Sort by confidence descending, keep top max_neighbors unique neighbors
#     neighbor_edges.sort(key=lambda x: x[2], reverse=True)

#     seen_neighbors = set()
#     kept_edges     = []
#     for u, v, w in neighbor_edges:
#         nbr = v if u == seed else u
#         if nbr not in seen_neighbors:
#             seen_neighbors.add(nbr)
#             kept_edges.append((u, v))
#         if len(seen_neighbors) >= max_neighbors:
#             break

#     # Build the display subgraph
#     kept_nodes = {seed} | seen_neighbors
#     sub        = G.subgraph(kept_nodes).copy()
#     # Keep only edges that connect seed to its sampled neighbors
#     sub        = nx.edge_subgraph(
#                      sub,
#                      [(u, v) for u, v in kept_edges if sub.has_edge(u, v)]
#                  ).copy()

#     # ── Layout: radial / circular around seed ─────────────────────────────────
#     # Use spring layout with seed node pinned at center
#     pos         = nx.spring_layout(sub, seed=42, k=2.5)
#     pos[seed]   = np.array([0.0, 0.0])   # force seed to center

#     # ── Figure ────────────────────────────────────────────────────────────────
#     fig, ax = plt.subplots(figsize=(14, 12))
#     fig.patch.set_facecolor("#0f1117")
#     ax.set_facecolor("#1a1d27")

#     # Node colors
#     node_colors = ["#e74c3c" if n == seed else "#00d4ff"
#                    for n in sub.nodes()]
#     node_sizes  = [1400 if n == seed else 280
#                    for n in sub.nodes()]

#     # Edge colors and widths by confidence
#     edge_weights = [sub[u][v]["weight"] for u, v in sub.edges()]
#     edge_colors  = [plt.cm.YlOrRd(w) for w in edge_weights]
#     edge_widths  = [0.8 + 2.5 * w     for w in edge_weights]

#     # Draw nodes
#     nx.draw_networkx_nodes(sub, pos, node_color=node_colors,
#                            node_size=node_sizes, alpha=0.92, ax=ax)

#     # Draw edges individually (needed for per-edge width)
#     for (u, v), col, wid in zip(sub.edges(), edge_colors, edge_widths):
#         nx.draw_networkx_edges(sub, pos,
#                                edgelist=[(u, v)],
#                                edge_color=[col],
#                                width=wid,
#                                arrows=True,
#                                arrowsize=10,
#                                connectionstyle="arc3,rad=0.1",
#                                ax=ax)

#     # Labels: show all (they are few enough now)
#     nx.draw_networkx_labels(sub, pos, font_size=6,
#                             font_color="white", ax=ax)

#     # Colorbar
#     sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd,
#                                norm=plt.Normalize(vmin=0, vmax=1))
#     sm.set_array([])
#     cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
#     cbar.set_label("Edge Confidence Score", color="#e0e0e0", fontsize=9)
#     cbar.ax.yaxis.set_tick_params(color="#e0e0e0", labelsize=7)
#     plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#e0e0e0")

#     # Legend
#     import matplotlib.patches as mpatches
#     legend_patches = [
#         mpatches.Patch(color="#e74c3c", label=f"Hub: {seed}"),
#         mpatches.Patch(color="#00d4ff", label="Direct neighbor"),
#         mpatches.Patch(color="#888888", label="Edge thickness ∝ confidence"),
#     ]
#     ax.legend(handles=legend_patches, loc="upper left",
#               facecolor="#1a1d27", edgecolor="#2e3250",
#               labelcolor="#e0e0e0", fontsize=9)

#     ax.set_title(
#         f"Ego Network of Top Hub Protein: {seed}\n"
#         f"Showing top {len(seen_neighbors)} neighbors "
#         f"(by confidence)  |  Edge color & thickness = confidence score",
#         color="#e0e0e0", fontsize=12, fontweight="bold", pad=12
#     )
#     ax.axis("off")
#     plt.tight_layout()

#     out = os.path.join(FIG_DIR, "graph_network.png")
#     plt.savefig(out, dpi=150, bbox_inches="tight",
#                 facecolor=fig.get_facecolor())
#     plt.close()
#     print(f"  Network figure saved    → {out}\n")
def _plot_ego_network(G: nx.DiGraph, stats: dict, max_neighbors: int = 60):
    print(f"[plot_network] Saving ego subgraph (max {max_neighbors} neighbors)...")

    seed = stats["top_proteins"][0][0]

    # ── Sample neighbors: keep top N by edge confidence ───────────────────────
    neighbor_edges = []
    for nbr in G.successors(seed):
        neighbor_edges.append((seed, nbr, G[seed][nbr]["weight"]))
    for nbr in G.predecessors(seed):
        neighbor_edges.append((nbr, seed, G[nbr][seed]["weight"]))

    neighbor_edges.sort(key=lambda x: x[2], reverse=True)

    seen_neighbors = set()
    for u, v, w in neighbor_edges:
        nbr = v if u == seed else u
        if nbr not in seen_neighbors:
            seen_neighbors.add(nbr)
        if len(seen_neighbors) >= max_neighbors:
            break

    # ── ✅ KEY FIX: keep the FULL subgraph — all edges among all kept nodes ────
    kept_nodes = {seed} | seen_neighbors
    sub = G.subgraph(kept_nodes).copy()   # ALL edges within this node set

    print(f"  Nodes: {sub.number_of_nodes()}  |  "
          f"Edges (including inter-neighbor): {sub.number_of_edges()}")

    # ── Layout ────────────────────────────────────────────────────────────────
    # Use kamada_kawai — much better than spring for showing inter-connections
    try:
        pos = nx.kamada_kawai_layout(sub, weight="weight")
    except Exception:
        pos = nx.spring_layout(sub, seed=42, k=1.8)
    pos[seed] = np.array([0.0, 0.0])

    # ── Separate edge types for different styling ─────────────────────────────
    hub_edges      = [(u, v) for u, v in sub.edges() if seed in (u, v)]
    internal_edges = [(u, v) for u, v in sub.edges() if seed not in (u, v)]

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 14))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    # Node styling
    node_colors = ["#e74c3c" if n == seed else "#00d4ff" for n in sub.nodes()]
    node_sizes  = [1600     if n == seed else 320        for n in sub.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(sub, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.92, ax=ax)

    # ── Draw hub edges (seed ↔ neighbor) — thick, colored by confidence ───────
    hub_weights = [sub[u][v]["weight"] for u, v in hub_edges]
    hub_colors  = [plt.cm.YlOrRd(w)   for w in hub_weights]
    hub_widths  = [0.8 + 2.5 * w      for w in hub_weights]
    for (u, v), col, wid in zip(hub_edges, hub_colors, hub_widths):
        nx.draw_networkx_edges(sub, pos,
                               edgelist=[(u, v)],
                               edge_color=[col],
                               width=wid,
                               arrows=True,
                               arrowsize=10,
                               connectionstyle="arc3,rad=0.08",
                               ax=ax)

    # ── Draw inter-neighbor edges — thinner, distinct color ───────────────────
    if internal_edges:
        int_weights = [sub[u][v]["weight"] for u, v in internal_edges]
        int_colors  = [plt.cm.cool(w)      for w in int_weights]
        int_widths  = [0.4 + 1.2 * w       for w in int_weights]
        for (u, v), col, wid in zip(internal_edges, int_colors, int_widths):
            nx.draw_networkx_edges(sub, pos,
                                   edgelist=[(u, v)],
                                   edge_color=[col],
                                   width=wid,
                                   arrows=True,
                                   arrowsize=7,
                                   connectionstyle="arc3,rad=0.15",
                                   ax=ax)

    # Labels
    nx.draw_networkx_labels(sub, pos, font_size=6,
                            font_color="white", ax=ax)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.YlOrRd,
                               norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Hub Edge Confidence Score", color="#e0e0e0", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="#e0e0e0", labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#e0e0e0")

    # Legend
    import matplotlib.patches as mpatches
    legend_patches = [
        mpatches.Patch(color="#e74c3c", label=f"Hub: {seed}"),
        mpatches.Patch(color="#00d4ff", label="Direct neighbor"),
        mpatches.Patch(color="#ff6600", label="Hub ↔ neighbor edge (YlOrRd)"),
        mpatches.Patch(color="#00ccff", label="Neighbor ↔ neighbor edge (cool)"),
    ]
    ax.legend(handles=legend_patches, loc="upper left",
              facecolor="#1a1d27", edgecolor="#2e3250",
              labelcolor="#e0e0e0", fontsize=9)

    ax.set_title(
        f"Ego Network of Top Hub Protein: {seed}\n"
        f"Nodes: {sub.number_of_nodes()}  |  "
        f"Hub edges: {len(hub_edges)}  |  "
        f"Inter-neighbor edges: {len(internal_edges)}",
        color="#e0e0e0", fontsize=12, fontweight="bold", pad=12
    )
    ax.axis("off")
    plt.tight_layout()

    out = os.path.join(FIG_DIR, "graph_network.png")
    plt.savefig(out, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Network figure saved → {out}\n")

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 5 — Public wrapper called by main.py
# ══════════════════════════════════════════════════════════════════════════════
def plot_statistics(G: nx.DiGraph, stats: dict, max_neighbors: int = 60):
    """
    Generates and saves both figures:
        1. graph_statistics.png  — statistics dashboard
        2. graph_network.png     — ego subgraph

    Parameters
    ----------
    G             : nx.DiGraph
    stats         : dict — output of compute_statistics()
    max_neighbors : int  — neighbors shown in the ego plot (default 60)
    """
    _plot_stats_dashboard(stats)
    _plot_ego_network(G, stats, max_neighbors=max_neighbors)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    G     = load_graph(DATA_PATH)
    stats = compute_statistics(G)
    plot_statistics(G, stats)
    print("✅ load_graph.py complete!\n")