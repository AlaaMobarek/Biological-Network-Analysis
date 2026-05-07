"""
main.py
───────
Project Entry Point — Protein-Protein Interaction Network Analysis

Runs the full analysis pipeline in order:
    1. Load the interactome graph
    2. Compute & visualize full network statistics
    3. Find shortest paths between two proteins
    4. Analyze neighbors of one protein
    5. Analyze degrees of a set of proteins
    6. Map UniProt IDs to gene names
    7. Save the adjacency matrix

Usage:
    python main.py
"""

import sys
import os

# Add scripts/ folder to path so imports work
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "Code", "scripts"))

from load_graph        import load_graph, compute_statistics, plot_statistics
from shortest_path     import find_shortest_paths
from neighbors         import get_neighbors
from degree_analysis   import analyze_degrees
from mapping           import map_uniprot_to_gene
from adjacency_matrix  import save_adjacency_matrix


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — Edit these values to customize the analysis
# ══════════════════════════════════════════════════════════════════════════════

# Path to the interactome file
DATA_PATH = "data/interactome.txt"

# Shortest path: choose any two proteins in the dataset
SOURCE_PROTEIN = "P04637"    # TP53
TARGET_PROTEIN = "P00533"    # EGFR

# Neighbor analysis: choose one protein
QUERY_PROTEIN  = "P04637"    # TP53

# Degree analysis: choose a set of proteins
PROTEIN_SET = [
    "P04637",   # TP53
    "P00533",   # EGFR
    "P06400",   # RB1
    "Q09472",   # EP300
    "P38936",   # CDKN1A
    "P12931",   # SRC
    "P27986",   # PIK3R1
    "P62993",   # GRB2
    "O14543",   # SOCS3
    "P43403",   # ZAP70
]

# UniProt → gene mapping: can be any set of IDs
MAPPING_IDS = PROTEIN_SET   # reuse the same set, or define a different list

# Adjacency matrix: number of top hub proteins to include in the sample matrix
MATRIX_SAMPLE_SIZE = 50


# ══════════════════════════════════════════════════════════════════════════════
# PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("\n" + "═" * 60)
    print("   PPI NETWORK ANALYSIS — FULL PIPELINE")
    print("═" * 60)

    # ── Step 1: Load graph ────────────────────────────────────────────────────
    print("\n[1/6] Loading interactome graph...")
    G = load_graph(DATA_PATH)

    # ── Step 2: Network statistics & visualization ────────────────────────────
    print("[2/6] Computing and plotting network statistics...")
    stats = compute_statistics(G)
    plot_statistics(G, stats)

    # ── Step 3: Shortest path analysis ────────────────────────────────────────
    print("[3/6] Finding shortest paths...")
    find_shortest_paths(G, source=SOURCE_PROTEIN, target=TARGET_PROTEIN)

    # ── Step 4: Neighbor analysis ─────────────────────────────────────────────
    print("[4/6] Analyzing neighbors...")
    get_neighbors(G, protein=QUERY_PROTEIN)

    # ── Step 5: Degree analysis ───────────────────────────────────────────────
    print("[5/6] Analyzing protein degrees...")
    analyze_degrees(G, proteins=PROTEIN_SET)

    # ── Step 6: UniProt → Gene mapping ────────────────────────────────────────
    print("[6/6] Mapping UniProt IDs to gene names...")
    map_uniprot_to_gene(MAPPING_IDS)

    # ── Step 7: Adjacency matrix ──────────────────────────────────────────────
    print("[7/7] Saving adjacency matrix...")
    save_adjacency_matrix(G, sample_size=MATRIX_SAMPLE_SIZE)

    # ── Done ──────────────────────────────────────────────────────────────────
    print("\n" + "═" * 60)
    print("   ✅ ALL STEPS COMPLETE")
    print("═" * 60)
    print("\n  Output files:")
    print("    Code/outputs/figures/graph_overview.png")
    print("    Code/outputs/figures/shortest_path_subnetwork.png")
    print("    Code/outputs/figures/neighbors_subgraph.png")
    print("    Code/outputs/figures/degree_histogram.png")
    print("    Code/outputs/figures/adjacency_matrix_heatmap.png")
    print("    Code/outputs/txt/shortest_paths.txt")
    print("    Code/outputs/txt/neighbors.txt")
    print("    Code/outputs/txt/degree_ranking.txt")
    print("    Code/outputs/txt/uniprot_gene_map.txt")
    print("    Code/outputs/matrices/adjacency_matrix.csv")
    print("    Code/outputs/matrices/adjacency_edgelist.txt")
    print()


if __name__ == "__main__":
    main()