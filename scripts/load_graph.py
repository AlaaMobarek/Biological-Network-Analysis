"""
load_graph.py
─────────────
Member 1 — Graph Construction & Visualization

Responsibilities:
    - Load the full PathLinker human PPI interactome into a directed NetworkX graph
    - Compute and print full network statistics
    - Visualize a sample ego subgraph + full statistics dashboard
    - Save the figure to outputs/figures/graph_overview.png

Usage (standalone):
    python scripts/load_graph.py

Usage (from main.py):
    from scripts.load_graph import load_graph
    G = load_graph("data/interactome.txt")
"""

import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


# ── Constants ──────────────────────────────────────────────────────────────────
DATA_PATH = "data/interactome.txt"
FIG_DIR   = "outputs/figures"
os.makedirs(FIG_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Load the interactome into a NetworkX DiGraph
# ══════════════════════════════════════════════════════════════════════════════
def load_graph(path: str) -> nx.DiGraph:
    """
    Reads the PathLinker interactome file and builds a directed weighted graph.

    Parameters
    ----------
    path : str
        Path to the interactome .txt file (tab-separated, 4 columns).

    Returns
    -------
    G : nx.DiGraph
        Directed graph where:
            - nodes  = proteins (UniProt IDs)
            - edges  = interactions with 'weight' and 'method' attributes
    """
    print(f"\n[load_graph] Reading file: {path}")

    G = nx.DiGraph()  # Directed graph: tail → head

    with open(path, "r") as f:
        for line in f:
            # Skip header / comment lines
            if line.startswith("#") or line.strip() == "":
                continue

            parts = line.strip().split("\t")

            # Need at least 3 columns: tail, head, weight
            if len(parts) < 3:
                continue

            tail   = parts[0]                               # Source protein
            head   = parts[1]                               # Target protein
            weight = float(parts[2])                        # Confidence [0,1]
            method = parts[3] if len(parts) > 3 else "unknown"  # Method(s)

            # Add directed edge with attributes stored on the edge
            G.add_edge(tail, head, weight=weight, method=method)

    print(f"[load_graph] Graph loaded successfully.\n")
    return G


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Compute and print full network statistics
# ══════════════════════════════════════════════════════════════════════════════
def compute_statistics(G: nx.DiGraph) -> dict:
    """
    Computes comprehensive statistics for the full interactome graph.

    Parameters
    ----------
    G : nx.DiGraph
        The loaded directed graph.

    Returns
    -------
    stats : dict
        Dictionary containing all computed statistics.
    """
    print("[compute_statistics] Computing network statistics...")

    # ── Basic counts ──────────────────────────────────────────────────────────
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # ── Degree statistics ─────────────────────────────────────────────────────
    degrees     = [d for _, d in G.degree()]       # total degree (in + out)
    in_degrees  = [d for _, d in G.in_degree()]    # in-degree only
    out_degrees = [d for _, d in G.out_degree()]   # out-degree only

    avg_degree    = np.mean(degrees)
    max_degree    = max(degrees)
    min_degree    = min(degrees)
    median_degree = np.median(degrees)

    # ── Edge weight statistics ────────────────────────────────────────────────
    weights    = [data["weight"] for _, _, data in G.edges(data=True)]
    avg_weight = np.mean(weights)
    max_weight = max(weights)
    min_weight = min(weights)

    # ── Connectivity ──────────────────────────────────────────────────────────
    weakly_connected = nx.is_weakly_connected(G)
    num_wcc          = nx.number_weakly_connected_components(G)
    largest_wcc_size = len(max(nx.weakly_connected_components(G), key=len))

    # ── Density ───────────────────────────────────────────────────────────────
    density = nx.density(G)

    # ── Top 10 hub proteins by total degree ───────────────────────────────────
    top_proteins = sorted(G.degree(), key=lambda x: x[1], reverse=True)[:10]

    # ── Bundle everything into a dictionary ───────────────────────────────────
    stats = {
        "num_nodes"        : num_nodes,
        "num_edges"        : num_edges,
        "avg_degree"       : avg_degree,
        "max_degree"       : max_degree,
        "min_degree"       : min_degree,
        "median_degree"    : median_degree,
        "avg_weight"       : avg_weight,
        "max_weight"       : max_weight,
        "min_weight"       : min_weight,
        "weakly_connected" : weakly_connected,
        "num_wcc"          : num_wcc,
        "largest_wcc_size" : largest_wcc_size,
        "density"          : density,
        "top_proteins"     : top_proteins,
        "degrees"          : degrees,
        "in_degrees"       : in_degrees,
        "out_degrees"      : out_degrees,
        "weights"          : weights,
    }

    # ── Pretty-print to console ───────────────────────────────────────────────
    print(f"\n{'═'*52}")
    print(f"   FULL NETWORK STATISTICS")
    print(f"{'═'*52}")
    print(f"  Nodes (proteins)            : {num_nodes:,}")
    print(f"  Edges (interactions)        : {num_edges:,}")
    print(f"  Graph density               : {density:.2e}")
    print(f"  Is directed                 : Yes")
    print(f"  Weakly connected            : {weakly_connected}")
    print(f"  Weakly connected components : {num_wcc:,}")
    print(f"  Largest component size      : {largest_wcc_size:,}")
    print(f"{'─'*52}")
    print(f"  Degree Statistics:")
    print(f"    Average degree            : {avg_degree:.2f}")
    print(f"    Median  degree            : {median_degree:.0f}")
    print(f"    Max     degree            : {max_degree}")
    print(f"    Min     degree            : {min_degree}")
    print(f"{'─'*52}")
    print(f"  Edge Weight Statistics:")
    print(f"    Average weight            : {avg_weight:.4f}")
    print(f"    Max     weight            : {max_weight:.4f}")
    print(f"    Min     weight            : {min_weight:.4f}")
    print(f"{'─'*52}")
    print(f"  Top 5 Hub Proteins (by degree):")
    for i, (node, deg) in enumerate(top_proteins[:5], 1):
        print(f"    {i}. {node}  →  degree {deg}")
    print(f"{'═'*52}\n")

    return stats


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Plot full network statistics dashboard + ego subgraph
# ══════════════════════════════════════════════════════════════════════════════
def plot_statistics(G: nx.DiGraph, stats: dict):
    """
    Creates a 7-panel statistics dashboard for the full interactome:
        Panel 1 : Total degree distribution (log scale)
        Panel 2 : In-degree distribution
        Panel 3 : Out-degree distribution
        Panel 4 : Edge weight distribution
        Panel 5 : Top 10 hub proteins (bar chart)
        Panel 6 : Network summary text box
        Panel 7 : Ego subgraph of the top hub protein (full bottom row)

    Saves figure to: outputs/figures/graph_overview.png

    Parameters
    ----------
    G     : nx.DiGraph  — the loaded graph
    stats : dict        — output of compute_statistics()
    """
    print("[plot_statistics] Generating statistics dashboard...")

    # ── Figure & dark theme setup ─────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor("#0f1117")

    # 3-row, 3-column grid; bottom row spans all 3 columns for ego graph
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.35)

    # Color palette
    ACCENT   = "#00d4ff"   # cyan
    ACCENT2  = "#ff6b6b"   # coral/red  (used for seed node)
    ACCENT3  = "#a8ff78"   # lime green
    BG_PANEL = "#1a1d27"   # panel background
    TEXT_COL = "#e0e0e0"   # axis text

    def style_ax(ax, title):
        """Apply consistent dark-theme styling to an axes panel."""
        ax.set_facecolor(BG_PANEL)
        ax.set_title(title, color=TEXT_COL, fontsize=11,
                     fontweight="bold", pad=8)
        ax.tick_params(colors=TEXT_COL, labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor("#2e3250")

    # ── Panel 1 : Total degree distribution ───────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(stats["degrees"], bins=60, color=ACCENT,
             edgecolor="none", alpha=0.85, log=True)
    ax1.set_xlabel("Degree", color=TEXT_COL, fontsize=8)
    ax1.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
    ax1.axvline(stats["avg_degree"], color=ACCENT2, linewidth=1.5,
                linestyle="--", label=f"Mean = {stats['avg_degree']:.1f}")
    ax1.legend(fontsize=7, labelcolor=TEXT_COL,
               facecolor=BG_PANEL, edgecolor="none")
    style_ax(ax1, "Total Degree Distribution")

    # ── Panel 2 : In-degree distribution ──────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(stats["in_degrees"], bins=60, color=ACCENT3,
             edgecolor="none", alpha=0.85, log=True)
    ax2.set_xlabel("In-Degree", color=TEXT_COL, fontsize=8)
    ax2.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
    style_ax(ax2, "In-Degree Distribution")

    # ── Panel 3 : Out-degree distribution ─────────────────────────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(stats["out_degrees"], bins=60, color=ACCENT2,
             edgecolor="none", alpha=0.85, log=True)
    ax3.set_xlabel("Out-Degree", color=TEXT_COL, fontsize=8)
    ax3.set_ylabel("Count (log scale)", color=TEXT_COL, fontsize=8)
    style_ax(ax3, "Out-Degree Distribution")

    # ── Panel 4 : Edge weight distribution ────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.hist(stats["weights"], bins=50, color="#f7971e",
             edgecolor="none", alpha=0.85)
    ax4.set_xlabel("Edge Weight (confidence score)", color=TEXT_COL, fontsize=8)
    ax4.set_ylabel("Count", color=TEXT_COL, fontsize=8)
    ax4.axvline(stats["avg_weight"], color=ACCENT, linewidth=1.5,
                linestyle="--", label=f"Mean = {stats['avg_weight']:.3f}")
    ax4.legend(fontsize=7, labelcolor=TEXT_COL,
               facecolor=BG_PANEL, edgecolor="none")
    style_ax(ax4, "Edge Weight Distribution")

    # ── Panel 5 : Top 10 hub proteins horizontal bar chart ────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    top10  = stats["top_proteins"]                      # already sorted top-10
    labels = [n for n, _ in top10]
    values = [d for _, d in top10]
    colors = plt.cm.plasma(np.linspace(0.3, 0.9, len(labels)))
    bars   = ax5.barh(labels[::-1], values[::-1],
                      color=colors, edgecolor="none")
    ax5.set_xlabel("Total Degree", color=TEXT_COL, fontsize=8)
    for bar, val in zip(bars, values[::-1]):
        ax5.text(bar.get_width() + 0.3,
                 bar.get_y() + bar.get_height() / 2,
                 str(val), va="center", color=TEXT_COL, fontsize=7)
    style_ax(ax5, "Top 10 Hub Proteins")

    # ── Panel 6 : Summary text box ────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_facecolor(BG_PANEL)
    ax6.axis("off")
    summary_text = (
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
    ax6.text(0.05, 0.97, summary_text,
             transform=ax6.transAxes,
             fontsize=9, verticalalignment="top",
             fontfamily="monospace", color=ACCENT,
             bbox=dict(facecolor="#0f1117", edgecolor=ACCENT,
                       boxstyle="round,pad=0.6", linewidth=1))
    ax6.set_title("Statistics Summary", color=TEXT_COL,
                  fontsize=11, fontweight="bold", pad=8)

    # ── Panel 7 : Ego subgraph — full bottom row ───────────────────────────────
    ax7 = fig.add_subplot(gs[2, :])   # spans all 3 columns
    ax7.set_facecolor(BG_PANEL)

    # Use the highest-degree node as the ego seed
    seed = stats["top_proteins"][0][0]
    ego  = nx.ego_graph(G, seed, radius=1)

    pos = nx.spring_layout(ego, seed=42, k=0.65)

    # Color: seed = coral, neighbors = cyan
    node_colors = [ACCENT2 if n == seed else ACCENT for n in ego.nodes()]
    node_sizes  = [900    if n == seed else 200    for n in ego.nodes()]

    # Edge color mapped to weight
    edge_weights = [ego[u][v]["weight"] for u, v in ego.edges()]
    edge_colors  = [plt.cm.YlOrRd(w) for w in edge_weights]

    nx.draw_networkx_nodes(ego, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9, ax=ax7)
    nx.draw_networkx_edges(ego, pos, edge_color=edge_colors,
                           width=1.2, arrows=True, arrowsize=10,
                           connectionstyle="arc3,rad=0.08", ax=ax7)
    nx.draw_networkx_labels(ego, pos, font_size=5.5,
                            font_color="white", ax=ax7)

    # Colorbar for edge weights
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.YlOrRd,
        norm=plt.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax7, fraction=0.015, pad=0.01)
    cbar.set_label("Edge Weight (confidence)", color=TEXT_COL, fontsize=8)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COL, labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COL)

    ax7.set_title(
        f"Ego Subgraph of Top Hub Protein: {seed}  "
        f"({ego.number_of_nodes()} proteins, {ego.number_of_edges()} interactions)",
        color=TEXT_COL, fontsize=11, fontweight="bold", pad=8
    )
    ax7.axis("off")

    # ── Main figure title ─────────────────────────────────────────────────────
    fig.suptitle(
        "PathLinker Human PPI Interactome — Full Network Analysis Dashboard",
        fontsize=15, fontweight="bold", color=TEXT_COL, y=0.99
    )

    # ── Save figure ───────────────────────────────────────────────────────────
    out_path = os.path.join(FIG_DIR, "graph_overview.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[plot_statistics] Figure saved → {out_path}\n")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN — runs when executed directly: python scripts/load_graph.py
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    G     = load_graph(DATA_PATH)       # Step 1: load graph
    stats = compute_statistics(G)       # Step 2: compute statistics
    plot_statistics(G, stats)           # Step 3: plot dashboard + ego subgraph
    print("✅ load_graph.py complete!\n")