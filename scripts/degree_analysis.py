"""
degree_analysis.py
──────────────────
Given a set of proteins, draws a degree histogram and ranks them from
highest to lowest connectivity in a text file.

Output:
    outputs/txt/degree_ranking.txt
    outputs/figures/degree_histogram.png

Usage:
    # Pass proteins as space-separated arguments
    python scripts/degree_analysis.py --proteins P04637 P00533 Q9Y243 O14763

    # Or point to a file with one UniProt ID per line
    python scripts/degree_analysis.py --proteins-file my_proteins.txt
"""

import os
import sys
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from load_graph import load_graph


# ── Degree helpers ────────────────────────────────────────────────────────────

def get_degrees(G, proteins):
    """Return dict: protein → (total_degree, in_degree, out_degree)."""
    result = {}
    for p in proteins:
        if p in G:
            result[p] = (G.degree(p), G.in_degree(p), G.out_degree(p))
        else:
            print(f"[WARN] Protein '{p}' not found in graph – skipped.")
    return result


# ── Text report ───────────────────────────────────────────────────────────────

def write_degree_ranking(degrees, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    ranked = sorted(degrees.items(), key=lambda x: x[1][0], reverse=True)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("=" * 60 + "\n")
        fh.write("  DEGREE RANKING  (highest -> lowest connectivity)\n")
        fh.write("=" * 60 + "\n")
        fh.write(f"  {'Rank':<6}  {'Protein':<20}  "
                 f"{'Total':>7}  {'In':>6}  {'Out':>6}\n")
        fh.write(f"  {'-'*6}  {'-'*20}  {'-'*7}  {'-'*6}  {'-'*6}\n")
        for rank, (prot, (tot, ind, outd)) in enumerate(ranked, 1):
            fh.write(f"  {rank:<6}  {prot:<20}  {tot:>7}  {ind:>6}  {outd:>6}\n")
        fh.write("\n")
        fh.write(f"  Total proteins analysed : {len(ranked)}\n")

    print(f"[INFO] Ranking saved → {out_path}")


# ── Histogram ─────────────────────────────────────────────────────────────────

def draw_histogram(degrees, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    proteins = list(degrees.keys())
    totals   = [degrees[p][0] for p in proteins]
    ins      = [degrees[p][1] for p in proteins]
    outs     = [degrees[p][2] for p in proteins]

    x    = np.arange(len(proteins))
    w    = 0.28

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.patch.set_facecolor("#0f1117")

    # ── Bar chart (per-protein breakdown) ────────────────────────────────────
    ax = axes[0]
    ax.set_facecolor("#0f1117")
    bars_tot = ax.bar(x - w,   totals, w, label="Total",  color="#3498db", alpha=0.9)
    bars_in  = ax.bar(x,       ins,    w, label="In",     color="#2ecc71", alpha=0.9)
    bars_out = ax.bar(x + w,   outs,   w, label="Out",    color="#e74c3c", alpha=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(proteins, rotation=45, ha="right",
                        color="white", fontsize=8)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_xlabel("Protein (UniProt ID)", color="white", fontsize=10)
    ax.set_ylabel("Degree", color="white", fontsize=10)
    ax.set_title("Per-Protein Degree Breakdown", color="white", fontsize=12, pad=10)
    ax.tick_params(colors="white")
    ax.spines[:].set_color("#444")
    ax.legend(facecolor="#1a1a2e", edgecolor="none",
              labelcolor="white", fontsize=9)

    # Value labels on bars
    for bar in [*bars_tot, *bars_in, *bars_out]:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.15,
                    str(int(h)), ha="center", va="bottom",
                    color="white", fontsize=7)

    # ── Sorted total-degree bar (ranking view) ────────────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor("#0f1117")

    sorted_items = sorted(zip(proteins, totals), key=lambda x: x[1], reverse=True)
    sorted_prots = [i[0] for i in sorted_items]
    sorted_vals  = [i[1] for i in sorted_items]

    cmap   = plt.cm.cool
    colors = [cmap(v / max(sorted_vals)) for v in sorted_vals]
    bars2  = ax2.bar(range(len(sorted_prots)), sorted_vals, color=colors, alpha=0.9)

    ax2.set_xticks(range(len(sorted_prots)))
    ax2.set_xticklabels(sorted_prots, rotation=45, ha="right",
                         color="white", fontsize=8)
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax2.set_xlabel("Protein (ranked)", color="white", fontsize=10)
    ax2.set_ylabel("Total Degree", color="white", fontsize=10)
    ax2.set_title("Degree Ranking (High → Low)", color="white", fontsize=12, pad=10)
    ax2.tick_params(colors="white")
    ax2.spines[:].set_color("#444")

    for bar, val in zip(bars2, sorted_vals):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.15,
                 str(int(val)), ha="center", va="bottom",
                 color="white", fontsize=8)

    fig.suptitle("Protein Degree Analysis", color="white", fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[INFO] Histogram saved → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Degree histogram and ranking for a set of proteins.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--proteins",      nargs="+",
                       help="Space-separated UniProt IDs")
    group.add_argument("--proteins-file", metavar="FILE",
                       help="Path to file with one UniProt ID per line")
    parser.add_argument("--interactome",  default="data/interactome.txt")
    parser.add_argument("--out-txt",      default="outputs/txt/degree_ranking.txt")
    parser.add_argument("--out-fig",      default="outputs/figures/degree_histogram.png")
    args = parser.parse_args()

    # Collect protein list
    if args.proteins:
        proteins = args.proteins
    else:
        if not os.path.isfile(args.proteins_file):
            sys.exit(f"[ERROR] Proteins file not found: {args.proteins_file}")
        with open(args.proteins_file) as fh:
            proteins = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    print(f"[INFO] Proteins to analyse: {proteins}")

    G = load_graph(args.interactome)
    degrees = get_degrees(G, proteins)

    if not degrees:
        sys.exit("[ERROR] None of the provided proteins were found in the graph.")

    write_degree_ranking(degrees, args.out_txt)
    draw_histogram(degrees, args.out_fig)

    print("[DONE] degree_analysis.py finished.")


if __name__ == "__main__":
    main()