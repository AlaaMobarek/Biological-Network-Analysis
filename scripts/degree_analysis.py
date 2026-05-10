import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────────
TXT_DIR = "outputs/txt"
FIG_DIR = "outputs/figures"
os.makedirs(TXT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Analyze degrees of a set of proteins
# ══════════════════════════════════════════════════════════════════════════════
def analyze_degrees(G: nx.DiGraph, proteins: list):
    """
    Computes the degree (total, in, out) of each protein in the given set,
    ranks them from most to least connected, saves ranking to a text file,
    and draws a degree histogram.

    Parameters
    ----------
    G        : nx.DiGraph — the loaded interactome graph
    proteins : list       — list of UniProt IDs to analyze
    """
    print(f"\n[degree_analysis] Analyzing {len(proteins)} proteins...")

    # ── Validate and collect degrees ──────────────────────────────────────────
    results   = []
    not_found = []

    for protein in proteins:
        if protein not in G:
            not_found.append(protein)
            continue

        total_deg = G.degree(protein)
        in_deg    = G.in_degree(protein)
        out_deg   = G.out_degree(protein)

        results.append({
            "protein"   : protein,
            "total"     : total_deg,
            "in_degree" : in_deg,
            "out_degree": out_deg,
        })

    if not_found:
        print(f"  ⚠ Not found in graph: {', '.join(not_found)}")

    if not results:
        print("  ✗ No valid proteins found.")
        return

    # ── Rank by total degree (highest first) ──────────────────────────────────
    results.sort(key=lambda x: x["total"], reverse=True)

    # Print ranking to console
    print(f"  {'Rank':<5} {'Protein':<14} {'Total':>7} "
          f"{'In':>7} {'Out':>7}")
    print(f"  {'-'*46}")
    for i, r in enumerate(results, 1):
        print(f"  {i:<5} {r['protein']:<14} {r['total']:>7} "
              f"{r['in_degree']:>7} {r['out_degree']:>7}")

    # ── Save ranking to text file ─────────────────────────────────────────────
    out_txt = os.path.join(TXT_DIR, "degree_ranking.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"Degree Ranking Analysis\n")
        f.write(f"{'='*55}\n")
        f.write(f"Proteins analyzed : {len(results)}\n")
        f.write(f"{'='*55}\n\n")
        f.write(f"{'Rank':<5} {'Protein':<14} {'Total':>7} "
                f"{'In':>7} {'Out':>7}\n")
        f.write(f"{'-'*46}\n")
        for i, r in enumerate(results, 1):
            f.write(f"{i:<5} {r['protein']:<14} {r['total']:>7} "
                    f"{r['in_degree']:>7} {r['out_degree']:>7}\n")

    print(f"\n  Ranking saved -> {out_txt}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    _draw_degree_histogram(results)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Draw degree histogram (3-panel: total, in, out)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_degree_histogram(results: list):
    """
    Draws a 3-panel figure:
        Panel 1 : Grouped bar chart (total / in / out degree per protein)
        Panel 2 : Stacked bar chart (in vs out)
        Panel 3 : Ranked total-degree bar chart with color gradient
    """
    proteins  = [r["protein"]    for r in results]
    totals    = [r["total"]      for r in results]
    in_degs   = [r["in_degree"]  for r in results]
    out_degs  = [r["out_degree"] for r in results]
    x         = np.arange(len(proteins))

    # Dark theme
    fig = plt.figure(figsize=(18, 12))
    fig.patch.set_facecolor("#0f1117")
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    BG_PANEL = "#1a1d27"
    TEXT_COL = "#e0e0e0"

    def style_ax(ax, title):
        ax.set_facecolor(BG_PANEL)
        ax.set_title(title, color=TEXT_COL, fontsize=11,
                     fontweight="bold", pad=8)
        ax.tick_params(colors=TEXT_COL, labelsize=8)
        ax.xaxis.label.set_color(TEXT_COL)
        ax.yaxis.label.set_color(TEXT_COL)
        for spine in ax.spines.values():
            spine.set_edgecolor("#2e3250")

    # ── Panel 1 : Grouped bar (total / in / out) ──────────────────────────────
    ax1  = fig.add_subplot(gs[0, 0])
    w    = 0.28
    ax1.bar(x - w, totals,  width=w, color="#00d4ff", label="Total",  alpha=0.9)
    ax1.bar(x,     in_degs, width=w, color="#a8ff78", label="In",     alpha=0.9)
    ax1.bar(x + w, out_degs,width=w, color="#ff6b6b", label="Out",    alpha=0.9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(proteins, rotation=45, ha="right", fontsize=7)
    ax1.set_ylabel("Degree")
    ax1.legend(facecolor=BG_PANEL, edgecolor="none",
               labelcolor=TEXT_COL, fontsize=8)
    style_ax(ax1, "Grouped Degree (Total / In / Out)")

    # ── Panel 2 : Stacked bar (in vs out) ─────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.bar(x, in_degs,  color="#a8ff78", label="In-degree",  alpha=0.9)
    ax2.bar(x, out_degs, bottom=in_degs,
            color="#ff6b6b", label="Out-degree", alpha=0.9)
    ax2.set_xticks(x)
    ax2.set_xticklabels(proteins, rotation=45, ha="right", fontsize=7)
    ax2.set_ylabel("Degree")
    ax2.legend(facecolor=BG_PANEL, edgecolor="none",
               labelcolor=TEXT_COL, fontsize=8)
    style_ax(ax2, "Stacked In vs Out Degree")

    # ── Panel 3 : Ranked total degree (color gradient) ────────────────────────
    ax3 = fig.add_subplot(gs[1, :])   # full bottom row
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(proteins)))
    bars   = ax3.bar(x, totals, color=colors, edgecolor="none", alpha=0.9)
    ax3.set_xticks(x)
    ax3.set_xticklabels(proteins, rotation=45, ha="right", fontsize=8)
    ax3.set_ylabel("Total Degree")
    # Value labels on top of bars
    for bar, val in zip(bars, totals):
        ax3.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 0.3, str(val),
                 ha="center", va="bottom", color=TEXT_COL, fontsize=7)
    style_ax(ax3, "Ranked Total Degree (High → Low)")

    fig.suptitle(
        "Protein Degree Analysis",
        fontsize=15, fontweight="bold", color=TEXT_COL, y=0.99
    )

    out_path = os.path.join(FIG_DIR, "degree_histogram.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Figure saved  -> {out_path}\n")


