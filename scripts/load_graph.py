import os
import math
import pickle
import networkx as nx

# Project paths
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DEFAULT_DATA_FILE = os.path.join(
    PROJECT_ROOT,
    "data",
    "PathLinker_2018_human-ppi-weighted-cap0_75.txt",
)

DEFAULT_GRAPHML_PATH = os.path.join(
    PROJECT_ROOT,
    "outputs",
    "graphs",
    "full_interactome.graphml"
)

DEFAULT_GPICKLE_PATH = os.path.join(
    PROJECT_ROOT,
    "outputs",
    "graphs",
    "full_interactome.gpickle"
)


def _resolve_path(path):
    """Convert relative paths to absolute paths."""
    if path is None:
        return None

    if os.path.isabs(path):
        return path

    return os.path.join(PROJECT_ROOT, path)


def build_full_graph(filepath,
                     gene_names_dict=None,
                     pathway_proteins_set=None):
    """
    Build the full directed interactome graph from the dataset.

    Edge attributes:
        weight         -> confidence score
        confidence     -> same as weight
        cost           -> 1 - confidence
        neg_log_weight -> -log(confidence)

    Node attributes:
        gene_name
        in_pathway
    """

    G = nx.DiGraph()

    gene_names = gene_names_dict or {}
    pathway_proteins = pathway_proteins_set or set()

    filepath = _resolve_path(filepath)

    if not os.path.exists(filepath):
        print(f"Error: File {filepath} not found.")
        return None

    with open(filepath) as fh:

        for line in fh:

            line = line.strip()

            # Skip empty/comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) < 3:
                continue

            src, dst = parts[0], parts[1]

            try:
                conf = float(parts[2])

            except ValueError:
                continue

            # Add directed interaction
            G.add_edge(
                src,
                dst,
                weight=conf,
                confidence=conf,
                cost=round(1.0 - conf, 6),
                neg_log_weight=round(
                    -math.log(max(conf, 1e-9)),
                    6
                )
            )

    # Add node metadata
    for node in G.nodes():

        G.nodes[node]["gene_name"] = (
            gene_names.get(node, node)
        )

        G.nodes[node]["in_pathway"] = (
            node in pathway_proteins
        )

    return G


def save_full_graph(G,
                    out_graphml,
                    out_gpickle):
    """
    Save graph as:
        - GraphML (portable)
        - GPickle (fast reload)
    """

    out_graphml = _resolve_path(out_graphml)
    out_gpickle = _resolve_path(out_gpickle)

    os.makedirs(
        os.path.dirname(out_graphml),
        exist_ok=True
    )

    nx.write_graphml(G, out_graphml)

    with open(out_gpickle, "wb") as f:
        pickle.dump(
            G,
            f,
            protocol=pickle.HIGHEST_PROTOCOL
        )

    print(
        f"Successfully saved:\n"
        f" - {out_graphml}\n"
        f" - {out_gpickle}"
    )


def load_graph(txt_filepath=None,
               gpickle_path=None,
               gene_names_dict=None,
               pathway_proteins_set=None,
               rebuild=False,
               save_on_build=True):
    """
    Load cached graph if available,
    otherwise rebuild from dataset.
    """

    txt_filepath = _resolve_path(
        txt_filepath or DEFAULT_DATA_FILE
    )

    gpickle_path = _resolve_path(
        gpickle_path or DEFAULT_GPICKLE_PATH
    )

    # Load cached graph
    if (
        not rebuild
        and gpickle_path
        and os.path.exists(gpickle_path)
    ):

        with open(gpickle_path, "rb") as f:
            G = pickle.load(f)

        print(
            f"Loaded graph from pickle -> {gpickle_path}"
        )

        return G

    # Build graph from dataset
    G = build_full_graph(
        txt_filepath,
        gene_names_dict=gene_names_dict,
        pathway_proteins_set=pathway_proteins_set,
    )

    if G is None:
        return None

    # Save graph cache
    if gpickle_path and save_on_build:

        save_full_graph(
            G,
            DEFAULT_GRAPHML_PATH,
            gpickle_path
        )

    return G