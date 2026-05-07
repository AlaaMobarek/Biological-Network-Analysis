"""
load_graph.py
─────────────
Constructs a directed, weighted PPI graph from the PathLinker interactome.

Interactome columns (tab-separated):
    0 – tail / source protein  (UniProt ID)
    1 – head / destination protein (UniProt ID)
    2 – interaction confidence  [0 – 1]
    3 – detection method

Edge attributes stored:
    confidence  – raw value from file
    weight      – 1 / confidence  (Dijkstra minimises this)
    method      – detection method string

Usage:
    from load_graph import load_graph
    G = load_graph("data/interactome.txt")
"""

import os
import sys
import networkx as nx


def load_graph(filepath: str) -> nx.DiGraph:
    """Parse interactome and return a directed NetworkX graph."""
    if not os.path.isfile(filepath):
        sys.exit(f"[ERROR] Interactome file not found: {filepath}")

    G = nx.DiGraph()
    skipped = 0

    with open(filepath, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                print(f"[WARN] Line {lineno} has fewer than 3 columns – skipped.")
                skipped += 1
                continue

            tail   = parts[0].strip()
            head   = parts[1].strip()
            method = parts[3].strip() if len(parts) >= 4 else "unknown"

            try:
                confidence = float(parts[2])
            except ValueError:
                print(f"[WARN] Line {lineno}: cannot parse confidence '{parts[2]}' – skipped.")
                skipped += 1
                continue

            weight = 1.0 / confidence if confidence > 0 else float("inf")

            G.add_edge(tail, head,
                       confidence=confidence,
                       weight=weight,
                       method=method)

    print(f"[INFO] Graph loaded from '{filepath}'")
    print(f"       Nodes  : {G.number_of_nodes():,}")
    print(f"       Edges  : {G.number_of_edges():,}")
    if skipped:
        print(f"       Skipped: {skipped} line(s)")

    return G


# ── Standalone sanity-check ───────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Load PPI interactome and print statistics.")
    parser.add_argument("--interactome", default="data/interactome.txt")
    args = parser.parse_args()

    G = load_graph(args.interactome)

    degrees    = dict(G.degree())
    avg_deg    = sum(degrees.values()) / len(degrees) if degrees else 0
    max_node   = max(degrees, key=degrees.get) if degrees else "N/A"
    confs      = [d["confidence"] for _, _, d in G.edges(data=True)]

    print("\n── Graph Statistics ─────────────────────────────────────────────")
    print(f"  Is directed          : {G.is_directed()}")
    print(f"  Nodes                : {G.number_of_nodes():,}")
    print(f"  Edges                : {G.number_of_edges():,}")
    print(f"  Average degree       : {avg_deg:.2f}")
    print(f"  Highest-degree node  : {max_node}  (degree {degrees.get(max_node, 0):,})")
    if confs:
        print(f"  Min confidence       : {min(confs):.4f}")
        print(f"  Max confidence       : {max(confs):.4f}")
        print(f"  Mean confidence      : {sum(confs)/len(confs):.4f}")
    print("\n[OK] load_graph.py finished.")