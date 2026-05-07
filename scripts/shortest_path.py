"""
shortest_path.py
────────────────
Given two proteins (source & target), finds ALL simple (acyclic) shortest
paths between them in the directed PPI graph, reports results to a text file,
and draws the sub-network using NetworkX + Matplotlib.

Output files:
    outputs/txt/shortest_paths.txt
    outputs/figures/shortest_path_subnetwork.png

Usage:
    python scripts/shortest_path.py --source P04637 --target P00533
    python scripts/shortest_path.py --source P04637 --target P00533 \\
        --interactome data/interactome.txt
"""

import os
import sys
import argparse

import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from load_graph import load_graph


# ── Path helpers ──────────────────────────────────────────────────────────────

def path_total_weight(G, path):
    return sum(G[u][v]["weight"] for u, v in zip(path, path[1:]))

def path_total_confidence(G, path):
    return sum(G[u][v]["confidence"] for u, v in zip(path, path[1:]))


def find_all_shortest_simple_paths(G, source, target):
    """
    Return every simple (acyclic) path that achieves the minimum total weight.
    Uses nx.shortest_simple_paths (Yen's-style generator) which yields paths
    in non-decreasing weight order.
    """
    try:
        gen = nx.shortest_simple_paths(G, source, target, weight="weight")
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return []

    best_w  = None
    results = []
    for path in gen:
        w = path_total_weight(G, path)
        if best_w is None:
            best_w = w
        if w > best_w + 1e-9:
            break
        results.append(path)
    return results


# ── Text report ───────────────────────────────────────────────────────────────

def write_report(G, paths, source, target, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("  SHORTEST PATH(S) REPORT\n")
        fh.write("=" * 70 + "\n")
        fh.write(f"  Source protein : {source}\n")
        fh.write(f"  Target protein : {target}\n")
        fh.write(f"  Paths found    : {len(paths)}\n")
        fh.write("=" * 70 + "\n\n")

        if not paths:
            fh.write("No path found between the two proteins.\n")
            return

        best_w = path_total_weight(G, paths[0])
        fh.write(f"Minimum total weight  (sum of 1/confidence) : {best_w:.6f}\n\n")

        for idx, path in enumerate(paths, 1):
            tw = path_total_weight(G, path)
            tc = path_total_confidence(G, path)
            fh.write(f"-- Path {idx}  ({len(path)-1} edge(s)) {'-'*40}\n")
            fh.write(f"   Nodes              : {' -> '.join(path)}\n")
            fh.write(f"   Total weight       : {tw:.6f}\n")
            fh.write(f"   Total confidence   : {tc:.6f}\n")
            fh.write("   Edge breakdown:\n")
            for u, v in zip(path, path[1:]):
                d = G[u][v]
                fh.write(
                    f"     {u} -> {v}"
                    f"  |  confidence={d['confidence']:.4f}"
                    f"  weight(1/conf)={d['weight']:.4f}"
                    f"  method={d.get('method','unknown')}\n"
                )
            fh.write("\n")

    print(f"[INFO] Report saved  → {out_path}")


# ── Visualisation ─────────────────────────────────────────────────────────────

def draw_subnetwork(G, paths, source, target, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Collect all nodes / edges used by the paths
    all_nodes = set()
    path_edges = {}
    for pidx, path in enumerate(paths):
        all_nodes.update(path)
        for u, v in zip(path, path[1:]):
            path_edges.setdefault((u, v), []).append(pidx + 1)

    # Build sub-graph
    SG = nx.DiGraph()
    for n in all_nodes:
        SG.add_node(n)
    for (u, v), pidxs in path_edges.items():
        d = G[u][v]
        SG.add_edge(u, v,
                    confidence=d["confidence"],
                    weight=d["weight"],
                    paths=pidxs)

    # Layout
    try:
        pos = nx.kamada_kawai_layout(SG, weight="weight")
    except Exception:
        pos = nx.spring_layout(SG, seed=42)

    # Node colours
    node_colours = []
    for n in SG.nodes():
        if n == source:   node_colours.append("#2ecc71")
        elif n == target: node_colours.append("#e74c3c")
        else:             node_colours.append("#3498db")

    # Edge colours by confidence
    confs       = [SG[u][v]["confidence"] for u, v in SG.edges()]
    norm        = mcolors.Normalize(vmin=min(confs), vmax=max(confs))
    cmap        = cm.plasma
    edge_colours = [cmap(norm(c)) for c in confs]

    edge_labels = {(u, v): f"{SG[u][v]['confidence']:.3f}" for u, v in SG.edges()}

    # Draw
    fig, ax = plt.subplots(figsize=(14, 9))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    nx.draw_networkx_nodes(SG, pos, ax=ax,
                           node_color=node_colours, node_size=1400, alpha=0.92)
    nx.draw_networkx_labels(SG, pos, ax=ax,
                            font_size=8, font_color="white", font_weight="bold")
    nx.draw_networkx_edges(SG, pos, ax=ax,
                           edge_color=edge_colours, arrows=True,
                           arrowsize=22, arrowstyle="-|>", width=2.2,
                           connectionstyle="arc3,rad=0.08", alpha=0.88,
                           min_source_margin=28, min_target_margin=28)
    nx.draw_networkx_edge_labels(SG, pos, edge_labels=edge_labels, ax=ax,
                                 font_size=7, font_color="#f0f0f0",
                                 bbox=dict(boxstyle="round,pad=0.2",
                                           fc="#222", ec="none", alpha=0.65))

    # Colorbar for edge confidence
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.01)
    cbar.set_label("Edge confidence", color="white", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white", fontsize=8)

    legend_handles = [
        mpatches.Patch(color="#2ecc71", label=f"Source: {source}"),
        mpatches.Patch(color="#e74c3c", label=f"Target: {target}"),
        mpatches.Patch(color="#3498db", label="Intermediate"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              framealpha=0.4, labelcolor="white",
              facecolor="#1a1a2e", edgecolor="none", fontsize=9)

    ax.set_title(
        f"Shortest-path Sub-network  |  {source} -> {target}\n"
        f"{len(paths)} path(s)  -  edge labels = confidence",
        color="white", fontsize=12, pad=12)
    ax.axis("off")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[INFO] Figure saved  → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Find all acyclic shortest paths between two PPI proteins.")
    parser.add_argument("--source",      required=True, help="UniProt ID of source protein")
    parser.add_argument("--target",      required=True, help="UniProt ID of target protein")
    parser.add_argument("--interactome", default="data/interactome.txt")
    parser.add_argument("--out-txt",     default="outputs/txt/shortest_paths.txt")
    parser.add_argument("--out-fig",     default="outputs/figures/shortest_path_subnetwork.png")
    args = parser.parse_args()

    G = load_graph(args.interactome)

    for prot, role in [(args.source, "source"), (args.target, "target")]:
        if prot not in G:
            sys.exit(f"[ERROR] {role} protein '{prot}' not found in graph.")

    print(f"[INFO] Searching shortest paths: {args.source} → {args.target}")
    paths = find_all_shortest_simple_paths(G, args.source, args.target)

    if not paths:
        print(f"[WARN] No path found between {args.source} and {args.target}.")
    else:
        print(f"[INFO] Found {len(paths)} shortest path(s).")

    write_report(G, paths, args.source, args.target, args.out_txt)

    if paths:
        draw_subnetwork(G, paths, args.source, args.target, args.out_fig)

    print("\n[DONE] shortest_path.py finished.")


if __name__ == "__main__":
    main()