import os
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse

# ─────────────────────────────────────────────────────────────────────────────
# Output directories
# ─────────────────────────────────────────────────────────────────────────────
MATRIX_DIR = "outputs/matrices"
FIG_DIR = "outputs/figures"

os.makedirs(MATRIX_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


# ═════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Convert to unweighted graph
# ═════════════════════════════════════════════════════════════════════════════
def convert_to_unweighted(G: nx.DiGraph) -> nx.DiGraph:
    """
    Convert weighted directed graph to unweighted directed graph.

    Parameters
    ----------
    G : nx.DiGraph
        Original weighted interactome graph

    Returns
    -------
    nx.DiGraph
        Unweighted directed graph
    """

    print("\n[1] Converting weighted graph to unweighted graph...")

    G_unweighted = nx.DiGraph()

    for u, v in G.edges():
        G_unweighted.add_edge(u, v)

    print(f"    Nodes : {G_unweighted.number_of_nodes():,}")
    print(f"    Edges : {G_unweighted.number_of_edges():,}")

    return G_unweighted


# ═════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Save full adjacency matrix
# ═════════════════════════════════════════════════════════════════════════════
def save_full_adjacency_matrix(
    G_unweighted: nx.DiGraph,
    save_csv: bool = True
):
    """
    Generate and save FULL adjacency matrix.

    Saves:
        1. Dense CSV matrix
        2. Sparse NPZ matrix
        3. Node order mapping

    Parameters
    ----------
    G_unweighted : nx.DiGraph
        Unweighted directed graph

    save_csv : bool
        Whether to save huge CSV matrix
    """

    print("\n[2] Building FULL adjacency matrix...")

    # Get all nodes in consistent order
    all_nodes = list(G_unweighted.nodes())

    # Build sparse adjacency matrix
    adj_sparse = nx.adjacency_matrix(
        G_unweighted,
        nodelist=all_nodes
    )

    print("    Sparse matrix created")

    # ── Save sparse matrix ────────────────────────────────────────────────
    sparse_path = os.path.join(
        MATRIX_DIR,
        "full_adjacency_matrix_sparse.npz"
    )

    sparse.save_npz(sparse_path, adj_sparse)

    print(f"    Sparse matrix saved -> {sparse_path}")

    # ── Save node order mapping ───────────────────────────────────────────
    nodes_path = os.path.join(
        MATRIX_DIR,
        "matrix_node_order.txt"
    )

    with open(nodes_path, "w") as f:
        for i, node in enumerate(all_nodes):
            f.write(f"{i}\t{node}\n")

    print(f"    Node mapping saved  -> {nodes_path}")

    # ── Optionally save FULL CSV matrix ──────────────────────────────────
    if save_csv:

        print("\n    Converting sparse matrix to dense array...")
        print("    This may take time for large interactomes...")

        adj_dense = adj_sparse.toarray().astype(np.uint8)

        df = pd.DataFrame(
            adj_dense,
            index=all_nodes,
            columns=all_nodes
        )

        csv_path = os.path.join(
            MATRIX_DIR,
            "full_adjacency_matrix.csv"
        )

        df.to_csv(csv_path)

        print(f"    Full CSV saved      -> {csv_path}")

        print(f"    Matrix shape        : "
              f"{df.shape[0]} × {df.shape[1]}")

        print(f"    Non-zero entries    : "
              f"{int(adj_dense.sum()):,}")

    return adj_sparse, all_nodes


# ═════════════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Save unweighted edge list
# ═════════════════════════════════════════════════════════════════════════════
def save_edge_list(G_unweighted: nx.DiGraph):
    """
    Save unweighted graph as edge list.

    Format:
        source<TAB>target
    """

    print("\n[3] Saving unweighted edge list...")

    edge_list_path = os.path.join(
        MATRIX_DIR,
        "adjacency_edgelist.txt"
    )

    with open(edge_list_path, "w") as f:

        f.write("# Unweighted adjacency edge list\n")
        f.write("# source\ttarget\n")

        for u, v in G_unweighted.edges():
            f.write(f"{u}\t{v}\n")

    print(f"    Edge list saved -> {edge_list_path}")


# ═════════════════════════════════════════════════════════════════════════════
# FUNCTION 4 — Plot sample adjacency heatmap
# ═════════════════════════════════════════════════════════════════════════════
def plot_sample_heatmap(
    G_unweighted: nx.DiGraph,
    sample_size: int = 50
):
    """
    Plot heatmap for top hub proteins only.

    Full 18k×18k matrix is impossible to visualize clearly,
    so we visualize only a smaller subgraph sample.
    """

    print(f"\n[4] Building sample heatmap "
          f"(top {sample_size} hub proteins)...")

    # Top proteins by degree
    top_nodes = [
        n for n, _ in
        sorted(
            G_unweighted.degree(),
            key=lambda x: x[1],
            reverse=True
        )[:sample_size]
    ]

    # Create subgraph
    sub = G_unweighted.subgraph(top_nodes).copy()

    # Create adjacency matrix
    matrix = nx.to_numpy_array(
        sub,
        nodelist=top_nodes,
        dtype=np.uint8
    )

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 12))

    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    im = ax.imshow(
    matrix,
    cmap="Blues",
    aspect="auto",
    interpolation="nearest"
) #white = no interaction (0)
#blue = interaction exists (1)

    # Labels
    fontsize = max(4, int(80 / sample_size))

    ax.set_xticks(range(len(top_nodes)))
    ax.set_yticks(range(len(top_nodes)))

    ax.set_xticklabels(
        top_nodes,
        rotation=90,
        fontsize=fontsize,
        color="#e0e0e0"
    )

    ax.set_yticklabels(
        top_nodes,
        fontsize=fontsize,
        color="#e0e0e0"
    )

    # Colorbar
    cbar = plt.colorbar(
        im,
        ax=ax,
        fraction=0.03,
        pad=0.02
    )

    cbar.set_label(
        "Interaction Exists (1) / Absent (0)",
        color="#e0e0e0",
        fontsize=9
    )

    plt.setp(
        cbar.ax.yaxis.get_ticklabels(),
        color="#e0e0e0"
    )

    # Grid
    ax.set_xticks(
        np.arange(-0.5, len(top_nodes), 1),
        minor=True
    )

    ax.set_yticks(
        np.arange(-0.5, len(top_nodes), 1),
        minor=True
    )

    ax.grid(
        which="minor",
        color="#2e3250",
        linewidth=0.4
    )

    # Titles
    ax.set_title(
        f"Adjacency Matrix Heatmap — "
        f"Top {sample_size} Hub Proteins",
        color="#e0e0e0",
        fontsize=13,
        fontweight="bold",
        pad=12
    )

    ax.set_xlabel(
        "Target Protein",
        color="#e0e0e0"
    )

    ax.set_ylabel(
        "Source Protein",
        color="#e0e0e0"
    )

    plt.tight_layout()

    # Save figure
    out_path = os.path.join(
        FIG_DIR,
        "adjacency_matrix_heatmap.png"
    )

    plt.savefig(
        out_path,
        dpi=150,
        bbox_inches="tight",
        facecolor=fig.get_facecolor()
    )

    plt.close()

    print(f"    Heatmap saved -> {out_path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":

    from load_graph import load_graph

    print("\nLoading interactome graph...")

    G = load_graph()

    print("\nInteractome loaded successfully!")
    print(f"Nodes : {G.number_of_nodes():,}")
    print(f"Edges : {G.number_of_edges():,}")

    # ─────────────────────────────────────────────────────────────────────
    # Step 1 — Convert to unweighted graph
    # ─────────────────────────────────────────────────────────────────────
    G_unweighted = convert_to_unweighted(G)

    # ─────────────────────────────────────────────────────────────────────
    # Step 2 — Save adjacency matrix
    # ─────────────────────────────────────────────────────────────────────
    adj_sparse, nodes = save_full_adjacency_matrix(
        G_unweighted,
        save_csv=False
    )

    # NOTE:
    # save_csv=False is recommended for huge interactomes.
    # Change to True ONLY if your machine has enough RAM/storage.

    # ─────────────────────────────────────────────────────────────────────
    # Step 3 — Save edge list
    # ─────────────────────────────────────────────────────────────────────
    save_edge_list(G_unweighted)

    # ─────────────────────────────────────────────────────────────────────
    # Step 4 — Plot sample heatmap
    # ─────────────────────────────────────────────────────────────────────
    plot_sample_heatmap(
        G_unweighted,
        sample_size=50
    )

    print("\n✅ adjacency_matrix.py completed successfully!\n")
    # Load sparse matrix
    matrix = sparse.load_npz(
        "outputs/matrices/full_adjacency_matrix_sparse.npz"
    )

    print(matrix)

    print("Shape:", matrix.shape)

    print("Non-zero entries:", matrix.nnz)