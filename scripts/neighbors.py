"""
neighbors.py
────────────
Given one protein (UniProt ID), lists all directly connected proteins
(predecessors + successors in the directed graph) and reports the degree.

Output:
    outputs/txt/neighbors.txt

Usage:
    python scripts/neighbors.py --protein P04637
    python scripts/neighbors.py --protein P04637 --interactome data/interactome.txt
"""

import os
import sys
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from load_graph import load_graph


def get_neighbors(G, protein):
    """
    Return a dict with:
        'in'  – proteins that point TO this protein  (predecessors)
        'out' – proteins this protein points TO      (successors)
    Each entry maps neighbor_id → edge data dict.
    """
    result = {"in": {}, "out": {}}

    for pred in G.predecessors(protein):
        result["in"][pred] = G[pred][protein]

    for succ in G.successors(protein):
        result["out"][succ] = G[protein][succ]

    return result


def write_neighbors_report(G, protein, neighbors, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    in_deg    = G.in_degree(protein)
    out_deg   = G.out_degree(protein)
    total_deg = G.degree(protein)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("  NEIGHBOR REPORT\n")
        fh.write("=" * 70 + "\n")
        fh.write(f"  Protein            : {protein}\n")
        fh.write(f"  Total degree       : {total_deg}\n")
        fh.write(f"  In-degree          : {in_deg}  (proteins pointing TO this node)\n")
        fh.write(f"  Out-degree         : {out_deg} (proteins this node points TO)\n")
        fh.write("=" * 70 + "\n\n")

        # ── Outgoing neighbors ────────────────────────────────────────────────
        fh.write(f"OUTGOING (successors, {out_deg} proteins)\n")
        fh.write("-" * 70 + "\n")
        fh.write(f"  {'Neighbor':<20}  {'Confidence':>12}  {'Weight(1/conf)':>15}  Method\n")
        fh.write(f"  {'-'*20}  {'-'*12}  {'-'*15}  {'-'*20}\n")
        for nbr, data in sorted(neighbors["out"].items(),
                                 key=lambda x: x[1]["confidence"], reverse=True):
            fh.write(
                f"  {nbr:<20}  {data['confidence']:>12.4f}"
                f"  {data['weight']:>15.4f}  {data.get('method','unknown')}\n"
            )

        fh.write("\n")

        # ── Incoming neighbors ────────────────────────────────────────────────
        fh.write(f"INCOMING (predecessors, {in_deg} proteins)\n")
        fh.write("-" * 70 + "\n")
        fh.write(f"  {'Neighbor':<20}  {'Confidence':>12}  {'Weight(1/conf)':>15}  Method\n")
        fh.write(f"  {'-'*20}  {'-'*12}  {'-'*15}  {'-'*20}\n")
        for nbr, data in sorted(neighbors["in"].items(),
                                  key=lambda x: x[1]["confidence"], reverse=True):
            fh.write(
                f"  {nbr:<20}  {data['confidence']:>12.4f}"
                f"  {data['weight']:>15.4f}  {data.get('method','unknown')}\n"
            )

    print(f"[INFO] Report saved  → {out_path}")
    print(f"[INFO] Protein {protein}:  degree={total_deg}  "
          f"(in={in_deg}, out={out_deg})")


def main():
    parser = argparse.ArgumentParser(
        description="List all directly connected proteins (neighbors) for a given protein.")
    parser.add_argument("--protein",     required=True, help="UniProt ID of the query protein")
    parser.add_argument("--interactome", default="data/interactome.txt")
    parser.add_argument("--out-txt",     default="outputs/txt/neighbors.txt")
    args = parser.parse_args()

    G = load_graph(args.interactome)

    if args.protein not in G:
        sys.exit(f"[ERROR] Protein '{args.protein}' not found in graph.")

    neighbors = get_neighbors(G, args.protein)
    write_neighbors_report(G, args.protein, neighbors, args.out_txt)

    print("[DONE] neighbors.py finished.")


if __name__ == "__main__":
    main()