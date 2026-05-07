"""
main.py
───────
Master runner — executes all 5 PPI network analyses in sequence.

Edit the CONFIGURATION section below to set your proteins, then run:
    python main.py
"""

import os
import sys

# ── Add scripts/ to path so we can import from it ────────────────────────────
SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")
sys.path.insert(0, SCRIPTS_DIR)

# ══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION  ← Edit this section before running
# ══════════════════════════════════════════════════════════════════════════════

INTERACTOME = "data/interactome.txt"

# Requirement 1 – Shortest path between two proteins
SOURCE_PROTEIN = "P04637"   # TP53
TARGET_PROTEIN = "P00533"   # EGFR

# Requirement 2 – Neighbors of one protein
QUERY_PROTEIN  = "P04637"   # TP53

# Requirements 3, 4, 5 – Set of proteins
PROTEIN_SET = [
    "P04637",   # TP53
    "P00533",   # EGFR
    "Q9Y243",   # AKT3
    "O14763",   # TNFRSF10B
    "P45983",   # MAPK8  (JNK1)
]

# Output paths
OUT_SHORTEST_TXT = "outputs/txt/shortest_paths.txt"
OUT_SHORTEST_FIG = "outputs/figures/shortest_path_subnetwork.png"
OUT_NEIGHBORS    = "outputs/txt/neighbors.txt"
OUT_DEGREE_TXT   = "outputs/txt/degree_ranking.txt"
OUT_DEGREE_FIG   = "outputs/figures/degree_histogram.png"
OUT_MAPPING      = "outputs/txt/uniprot_gene_map.txt"
OUT_MATRIX_CSV   = "outputs/matrices/adjacency_matrix.csv"
OUT_MATRIX_FIG   = "outputs/figures/graph_overview.png"

# ══════════════════════════════════════════════════════════════════════════════


def make_output_dirs():
    for d in ["outputs/txt", "outputs/figures", "outputs/matrices"]:
        os.makedirs(d, exist_ok=True)


def banner(title):
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


# ──────────────────────────────────────────────────────────────────────────────
# Step 0 — Load the graph (shared by all steps)
# ──────────────────────────────────────────────────────────────────────────────
def step0_load_graph():
    banner("STEP 0 — Load PPI Graph")
    from load_graph import load_graph
    G = load_graph(INTERACTOME)
    return G


# ──────────────────────────────────────────────────────────────────────────────
# Step 1 — Shortest paths
# ──────────────────────────────────────────────────────────────────────────────
def step1_shortest_paths(G):
    banner("STEP 1 — Shortest Path(s)")
    print(f"  Source : {SOURCE_PROTEIN}")
    print(f"  Target : {TARGET_PROTEIN}")

    from shortest_path import (
        find_all_shortest_simple_paths,
        write_report,
        draw_subnetwork,
    )

    for prot, role in [(SOURCE_PROTEIN, "source"), (TARGET_PROTEIN, "target")]:
        if prot not in G:
            print(f"  [SKIP] {role} protein '{prot}' not found in graph.")
            return

    paths = find_all_shortest_simple_paths(G, SOURCE_PROTEIN, TARGET_PROTEIN)

    if not paths:
        print(f"  [WARN] No path found between {SOURCE_PROTEIN} and {TARGET_PROTEIN}.")
    else:
        print(f"  Found {len(paths)} shortest path(s).")

    write_report(G, paths, SOURCE_PROTEIN, TARGET_PROTEIN, OUT_SHORTEST_TXT)

    if paths:
        draw_subnetwork(G, paths, SOURCE_PROTEIN, TARGET_PROTEIN, OUT_SHORTEST_FIG)


# ──────────────────────────────────────────────────────────────────────────────
# Step 2 — Neighbors
# ──────────────────────────────────────────────────────────────────────────────
def step2_neighbors(G):
    banner("STEP 2 — Neighbors")
    print(f"  Query protein : {QUERY_PROTEIN}")

    from neighbors import get_neighbors, write_neighbors_report

    if QUERY_PROTEIN not in G:
        print(f"  [SKIP] Protein '{QUERY_PROTEIN}' not found in graph.")
        return

    nbrs = get_neighbors(G, QUERY_PROTEIN)
    write_neighbors_report(G, QUERY_PROTEIN, nbrs, OUT_NEIGHBORS)


# ──────────────────────────────────────────────────────────────────────────────
# Step 3 — Degree analysis
# ──────────────────────────────────────────────────────────────────────────────
def step3_degree_analysis(G):
    banner("STEP 3 — Degree Analysis")
    print(f"  Proteins : {PROTEIN_SET}")

    from degree_analysis import get_degrees, write_degree_ranking, draw_histogram

    degrees = get_degrees(G, PROTEIN_SET)
    if not degrees:
        print("  [SKIP] None of the proteins found in graph.")
        return

    write_degree_ranking(degrees, OUT_DEGREE_TXT)
    draw_histogram(degrees, OUT_DEGREE_FIG)


# ──────────────────────────────────────────────────────────────────────────────
# Step 4 — UniProt → Gene name mapping
# ──────────────────────────────────────────────────────────────────────────────
def step4_mapping():
    banner("STEP 4 — UniProt → Gene Name Mapping")
    print(f"  Proteins : {PROTEIN_SET}")
    print("  (Requires internet access — queries UniProt REST API)")

    from mapping import fetch_gene_names, print_mapping, write_mapping_report

    mapping = fetch_gene_names(PROTEIN_SET)
    print_mapping(mapping)
    write_mapping_report(mapping, OUT_MAPPING)


# ──────────────────────────────────────────────────────────────────────────────
# Step 5 — Unweighted adjacency matrix
# ──────────────────────────────────────────────────────────────────────────────
def step5_adjacency_matrix(G):
    banner("STEP 5 — Unweighted Adjacency Matrix")
    print(f"  Proteins : {PROTEIN_SET}")

    from adjacency_matrix import (
        to_unweighted,
        save_adjacency_matrix_dense,
        draw_graph_overview,
    )

    UG    = to_unweighted(G)
    nodes = [p for p in PROTEIN_SET if p in UG]
    missing = [p for p in PROTEIN_SET if p not in UG]

    if missing:
        print(f"  [WARN] Not found (skipped): {missing}")
    if not nodes:
        print("  [SKIP] No proteins found in graph.")
        return

    save_adjacency_matrix_dense(UG, nodes, OUT_MATRIX_CSV)
    draw_graph_overview(UG, nodes, OUT_MATRIX_FIG)


# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────
def print_summary():
    banner("ALL DONE — Output Files")
    outputs = [
        ("Shortest paths text",   OUT_SHORTEST_TXT),
        ("Shortest path figure",  OUT_SHORTEST_FIG),
        ("Neighbors text",        OUT_NEIGHBORS),
        ("Degree ranking text",   OUT_DEGREE_TXT),
        ("Degree histogram",      OUT_DEGREE_FIG),
        ("UniProt gene map",      OUT_MAPPING),
        ("Adjacency matrix CSV",  OUT_MATRIX_CSV),
        ("Graph overview figure", OUT_MATRIX_FIG),
    ]
    for label, path in outputs:
        status = "✓" if os.path.isfile(path) else "✗ missing"
        print(f"  {status}  {label:<28}  {path}")
    print()


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    make_output_dirs()

    G = step0_load_graph()

    step1_shortest_paths(G)
    step2_neighbors(G)
    step3_degree_analysis(G)
    step4_mapping()            # needs internet
    step5_adjacency_matrix(G)

    print_summary()