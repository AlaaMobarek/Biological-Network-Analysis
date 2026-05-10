import os

# Import custom analysis and utility functions
from adjacency_matrix import (
    convert_to_unweighted,
    save_edge_list,
    save_full_adjacency_matrix,
    plot_sample_heatmap,
)
from degree_analysis import analyze_degrees
from graph_visualization import run_pathway_visualization, draw_top_hubs
from load_graph import (
    DEFAULT_DATA_FILE,
    DEFAULT_GRAPHML_PATH,
    DEFAULT_GPICKLE_PATH,
    load_graph,
    save_full_graph,
)
from mapping import map_uniprot_to_gene
from neighbors import get_neighbors
from shortest_path import run_shortest_paths_pipeline

# Initialize project directory structure
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_ROOT = os.path.join(PROJECT_ROOT, "outputs")
os.makedirs(os.path.join(OUTPUT_ROOT, "graphs"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_ROOT, "txt"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_ROOT, "figures"), exist_ok=True)

# Define target proteins and pathway for the pipeline
neighbor = "P10721"
SOURCE_PROTEIN = "P04637"  # TP53
TARGET_PROTEIN = "P00533"  # EGFR
PATHWAY_ID = "hsa04151"    # PI3K-AKT

PROTEINS_TO_ANALYZE = [
    "P04637", "P00533", "P06400", "Q09472", "P38936",
    "P12931", "P27986", "P62993", "O14543", "P43403",
]

def main():
    """Executes the full bioinformatics pipeline sequentially."""
    print(f"[main] Building graph from: {DEFAULT_DATA_FILE}")
    G = load_graph(DEFAULT_DATA_FILE, DEFAULT_GPICKLE_PATH, rebuild=True, save_on_build=False)
    if G is None:
        print("[main] Failed to load the interactome graph.")
        return

    # 1. Persistence: Save the graph for future use
    save_full_graph(G, DEFAULT_GRAPHML_PATH, DEFAULT_GPICKLE_PATH)

    # 2. Visualization: Create KEGG-based pathway map
    print("\n[main] Running pathway visualization...")
    run_pathway_visualization(
        G_full=G,
        pathway_id=PATHWAY_ID,
        source_id=SOURCE_PROTEIN,
        target_id=TARGET_PROTEIN,
        conf_threshold=0.75,
    )

    # 3. Network Metrics: Hub identification
    print("\n[main] Drawing top hub proteins...")
    draw_top_hubs(G, top_n=20)

    # 4. Path Analysis: Source-Target interaction paths
    print("\n[main] Running shortest-path analysis...")
    shortest_txt = os.path.join(OUTPUT_ROOT, "txt", "shortest_paths_result.txt")
    shortest_img = os.path.join(OUTPUT_ROOT, "figures", "shortest_paths_subnetwork.png")
    run_shortest_paths_pipeline(
        G=G,
        source=SOURCE_PROTEIN,
        target=TARGET_PROTEIN,
        output_txt=shortest_txt,
        output_img=shortest_img,
    )

    # 5. Local Neighborhood: 1-hop connections
    print("\n[main] Running neighbors analysis...")
    get_neighbors(G, protein=neighbor)

    # 6. Global Stats: Degree distribution and ranking
    print("\n[main] Running degree analysis...")
    analyze_degrees(G, PROTEINS_TO_ANALYZE)

    # 7. Annotation: Map IDs to human-readable names via API
    print("\n[main] Running mapping report...")
    map_uniprot_to_gene(PROTEINS_TO_ANALYZE)

    # 8. Matrix Export: Adjacency representations
    print("\n[main] Running adjacency-matrix export...")
    G_unweighted = convert_to_unweighted(G)
    save_full_adjacency_matrix(G_unweighted, save_csv=False)
    save_edge_list(G_unweighted)
    plot_sample_heatmap(G_unweighted, sample_size=50)

    print("\n[main] Pipeline complete.")

if __name__ == "__main__":
    main()