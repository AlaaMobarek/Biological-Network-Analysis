"""
adjacency_matrix.py
────────────────────
Converts the directed PPI graph into an UNWEIGHTED directed graph and
saves it as an adjacency matrix (CSV file).

Because the full interactome has thousands of nodes, the resulting matrix
can be very large. This script therefore supports:
    • Full export  (may be gigabytes – use carefully)
    • Subset export: provide a list of proteins to extract their sub-matrix
    • Optional sparse COO format  (row, col, value)  as an alternative

Output:
    outputs/matrices/adjacency_matrix.csv

Usage:
    # Export adjacency matrix for a subset of proteins
    python scripts/adjacency_matrix.py --proteins P04637 P00533 Q9Y243 O14763

    # Export the full interactome adjacency matrix (large!)
    python scripts/adjacency_matrix.py --full

    # Export sparse COO format for the full graph
    python scripts/adjacency_matrix.py --full --sparse
"""

import os
import sys
import argparse
import csv

import networkx as nx
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from load_graph import load_graph


# ── Graph conversion ──────────────────────────────────────────────────────────

def to_unweighted(G: nx.DiGraph) -> nx.DiGraph:
    """
    Return a copy of G with all edge weights removed (unweighted digraph).
    Edges become binary: 1 = interaction exists, 0 = no interaction.
    """
    UG = nx.DiGraph()
    UG.add_nodes_from(G.nodes())
    for u, v in G.edges():
        UG.add_edge(u, v)
    return UG


# ── Dense adjacency matrix ────────────────────────────────────────────────────

def save_adjacency_matrix_dense(G: nx.DiGraph, nodes: list, out_path: str) -> None:
    """
    Save the dense adjacency matrix of the sub-graph induced by `nodes`.
    Rows/columns are ordered as in `nodes`.
    Cell value: 1 if edge exists (tail→head), 0 otherwise.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Build sub-graph
    SG = G.subgraph(nodes).copy()

    n = len(nodes)
    idx = {node: i for i, node in enumerate(nodes)}

    # Build matrix
    matrix = np.zeros((n, n), dtype=int)
    for u, v in SG.edges():
        if u in idx and v in idx:
            matrix[idx[u]][idx[v]] = 1

    # Write CSV  (first row = header with node names, first col = row node name)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + nodes)              # header row
        for i, node in enumerate(nodes):
            writer.writerow([node] + matrix[i].tolist())

    print(f"[INFO] Dense adjacency matrix ({n}x{n}) saved -> {out_path}")
    print(f"       Non-zero entries: {int(matrix.sum())}")


# ── Sparse COO adjacency ──────────────────────────────────────────────────────

def save_adjacency_matrix_sparse(G: nx.DiGraph, out_path: str) -> None:
    """
    Save a sparse edge-list (COO format) adjacency representation:
        source, target, value
    This is memory-efficient for large graphs.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["source", "target", "value"])
        for u, v in G.edges():
            writer.writerow([u, v, 1])

    edge_count = G.number_of_edges()
    print(f"[INFO] Sparse COO adjacency matrix ({edge_count} entries) saved -> {out_path}")


# ── Graph overview figure ─────────────────────────────────────────────────────

def draw_graph_overview(G: nx.DiGraph, nodes: list, out_path: str) -> None:
    """Draw the sub-graph for the given nodes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    SG = G.subgraph(nodes).copy()

    try:
        pos = nx.kamada_kawai_layout(SG)
    except Exception:
        pos = nx.spring_layout(SG, seed=42)

    degrees     = dict(SG.degree())
    max_deg     = max(degrees.values()) if degrees else 1
    node_sizes  = [300 + 1200 * (degrees.get(n, 0) / max_deg) for n in SG.nodes()]
    norm        = mcolors.Normalize(vmin=0, vmax=max_deg)
    node_colours= [cm.cool(norm(degrees.get(n, 0))) for n in SG.nodes()]

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    nx.draw_networkx_nodes(SG, pos, ax=ax,
                           node_color=node_colours,
                           node_size=node_sizes, alpha=0.9)
    nx.draw_networkx_labels(SG, pos, ax=ax,
                            font_size=8, font_color="white", font_weight="bold")
    nx.draw_networkx_edges(SG, pos, ax=ax,
                           edge_color="#888", arrows=True,
                           arrowsize=15, arrowstyle="-|>", width=1.2,
                           connectionstyle="arc3,rad=0.06", alpha=0.6,
                           min_source_margin=20, min_target_margin=20)

    sm = cm.ScalarMappable(cmap=cm.cool, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.01)
    cbar.set_label("Node degree", color="white", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white", fontsize=8)

    ax.set_title(
        f"PPI Sub-graph Overview  |  {len(nodes)} proteins  ·  "
        f"{SG.number_of_edges()} interactions",
        color="white", fontsize=12, pad=12)
    ax.axis("off")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[INFO] Overview figure saved → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Save unweighted adjacency matrix of the PPI graph.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--proteins",   nargs="+",
                      help="UniProt IDs for a sub-matrix (recommended for large graphs)")
    mode.add_argument("--proteins-file", metavar="FILE",
                      help="File with one UniProt ID per line")
    mode.add_argument("--full",       action="store_true",
                      help="Export the full interactome (can be very large)")
    parser.add_argument("--sparse",   action="store_true",
                        help="Use sparse COO format instead of dense matrix")
    parser.add_argument("--interactome", default="data/interactome.txt")
    parser.add_argument("--out-csv",  default="outputs/matrices/adjacency_matrix.csv")
    parser.add_argument("--out-fig",  default="outputs/figures/graph_overview.png")
    args = parser.parse_args()

    G  = load_graph(args.interactome)
    UG = to_unweighted(G)
    print(f"[INFO] Unweighted graph: {UG.number_of_nodes():,} nodes, "
          f"{UG.number_of_edges():,} edges")

    # Decide which nodes to use
    if args.full:
        nodes = list(UG.nodes())
        print(f"[WARN] Full export: {len(nodes):,} nodes. "
              f"Dense matrix will be {len(nodes)**2 * 4 / 1e6:.1f} MB. "
              f"Consider --sparse for large graphs.")
        if args.sparse:
            save_adjacency_matrix_sparse(UG, args.out_csv)
        else:
            if len(nodes) > 5000:
                print("[WARN] Graph is very large (>5000 nodes). "
                      "Switching to sparse format automatically.")
                sparse_path = args.out_csv.replace(".csv", "_sparse.csv")
                save_adjacency_matrix_sparse(UG, sparse_path)
            else:
                save_adjacency_matrix_dense(UG, nodes, args.out_csv)

    else:
        if args.proteins:
            proteins = args.proteins
        else:
            if not os.path.isfile(args.proteins_file):
                sys.exit(f"[ERROR] File not found: {args.proteins_file}")
            with open(args.proteins_file) as fh:
                proteins = [l.strip() for l in fh
                            if l.strip() and not l.startswith("#")]

        # Filter to proteins present in the graph
        nodes = [p for p in proteins if p in UG]
        missing = [p for p in proteins if p not in UG]
        if missing:
            print(f"[WARN] Not in graph (skipped): {missing}")
        if not nodes:
            sys.exit("[ERROR] None of the provided proteins found in graph.")

        print(f"[INFO] Building sub-matrix for {len(nodes)} proteins: {nodes}")

        if args.sparse:
            SG = UG.subgraph(nodes).copy()
            save_adjacency_matrix_sparse(SG, args.out_csv)
        else:
            save_adjacency_matrix_dense(UG, nodes, args.out_csv)

        draw_graph_overview(UG, nodes, args.out_fig)

    print("[DONE] adjacency_matrix.py finished.")


if __name__ == "__main__":
    main()