"""
neighbors.py
────────────
Layout:
    - Query protein at CENTER (large red circle)
    - 6 clusters arranged around it like a clock face
    - Nodes are RANDOMLY SCATTERED inside each cluster disk (not in a ring)
    - Cluster background = soft colored disk
    - Arrow color matches the cluster color
    - Node color matches the cluster

Usage (standalone):   python scripts/neighbors.py
Usage (from main.py): from scripts.neighbors import get_neighbors
                      get_neighbors(G, protein="P04637")
"""

import os
import math
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────────
TXT_DIR   = "outputs/txt"
FIG_DIR   = "outputs/figures"
os.makedirs(TXT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

HIGH_CONF = 0.75
MED_CONF  = 0.50

# ── Cluster config: angle (°), base orbit radius, colors, label ───────────────
CLUSTERS = {
    "out_high": dict(angle=120, radius=5.0,
                     node_fc="#1a9e5c", node_ec="#0e6b3a",
                     edge_col="#27ae60",
                     label="Outgoing HIGH  (conf ≥ 0.75)"),
    "out_med":  dict(angle=60,  radius=5.0,
                     node_fc="#17a589", node_ec="#0e7260",
                     edge_col="#1abc9c",
                     label="Outgoing MED   (0.50–0.75)"),
    "out_low":  dict(angle=0,   radius=5.0,
                     node_fc="#d4ac0d", node_ec="#9a7d0a",
                     edge_col="#f1c40f",
                     label="Outgoing LOW   (conf < 0.50)"),
    "in_high":  dict(angle=300, radius=5.0,
                     node_fc="#1a5276", node_ec="#0e3457",
                     edge_col="#2980b9",
                     label="Incoming HIGH  (conf ≥ 0.75)"),
    "in_med":   dict(angle=240, radius=5.0,
                     node_fc="#6c3483", node_ec="#4a235a",
                     edge_col="#9b59b6",
                     label="Incoming MED   (0.50–0.75)"),
    "in_low":   dict(angle=180, radius=5.0,
                     node_fc="#a93226", node_ec="#7b241c",
                     edge_col="#e74c3c",
                     label="Incoming LOW   (conf < 0.50)"),
}


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Main entry point
# ══════════════════════════════════════════════════════════════════════════════
def get_neighbors(G: nx.DiGraph, protein: str):
    print(f"\n[neighbors] Analyzing neighbors of: {protein}")

    if protein not in G:
        print(f"  ⚠  Protein '{protein}' not found in graph.")
        return

    successors = []
    for nbr in G.successors(protein):
        w = G[protein][nbr]["weight"]
        m = G[protein][nbr].get("method", "unknown")
        successors.append((nbr, w, m))

    predecessors = []
    for nbr in G.predecessors(protein):
        w = G[nbr][protein]["weight"]
        m = G[nbr][protein].get("method", "unknown")
        predecessors.append((nbr, w, m))

    successors   = sorted(successors,   key=lambda x: x[1], reverse=True)
    predecessors = sorted(predecessors, key=lambda x: x[1], reverse=True)

    in_deg  = G.in_degree(protein)
    out_deg = G.out_degree(protein)
    print(f"  Total degree       : {in_deg+out_deg}  (in={in_deg}, out={out_deg})")
    print(f"  Outgoing neighbors : {len(successors)}")
    print(f"  Incoming neighbors : {len(predecessors)}")

    _save_results(protein, successors, predecessors, in_deg, out_deg)
    _draw_cluster_graph(G, protein, successors, predecessors)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Save text file
# ══════════════════════════════════════════════════════════════════════════════
def _save_results(protein, successors, predecessors, in_deg, out_deg):
    out_txt = os.path.join(TXT_DIR, "neighbors.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("Neighbor Analysis\n")
        f.write(f"{'='*60}\n")
        f.write(f"Query protein  : {protein}\n")
        f.write(f"Total degree   : {in_deg+out_deg}\n")
        f.write(f"  In-degree    : {in_deg}\n")
        f.write(f"  Out-degree   : {out_deg}\n")
        f.write(f"{'='*60}\n\n")

        f.write(f"OUTGOING ({len(successors)})  -- {protein} -> neighbor:\n")
        f.write(f"{'─'*60}\n")
        for nbr, w, m in successors:
            tier = "HIGH" if w >= HIGH_CONF else "MED" if w >= MED_CONF else "LOW"
            f.write(f"  {nbr:<12}  weight: {w:.6f}  [{tier}]  method: {m}\n")

        f.write(f"\nINCOMING ({len(predecessors)})  -- neighbor -> {protein}:\n")
        f.write(f"{'─'*60}\n")
        for nbr, w, m in predecessors:
            tier = "HIGH" if w >= HIGH_CONF else "MED" if w >= MED_CONF else "LOW"
            f.write(f"  {nbr:<12}  weight: {w:.6f}  [{tier}]  method: {m}\n")

    print(f"  Results saved -> {out_txt}")


# ══════════════════════════════════════════════════════════════════════════════
# HELPER — Assign neighbors into buckets
# ══════════════════════════════════════════════════════════════════════════════
def _assign_buckets(successors, predecessors):
    buckets = {k: [] for k in CLUSTERS}
    for nbr, w, _ in successors:
        if   w >= HIGH_CONF: buckets["out_high"].append((nbr, w))
        elif w >= MED_CONF:  buckets["out_med"].append((nbr, w))
        else:                buckets["out_low"].append((nbr, w))
    for nbr, w, _ in predecessors:
        if   w >= HIGH_CONF: buckets["in_high"].append((nbr, w))
        elif w >= MED_CONF:  buckets["in_med"].append((nbr, w))
        else:                buckets["in_low"].append((nbr, w))
    return buckets


# ══════════════════════════════════════════════════════════════════════════════
# HELPER — Scatter nodes randomly inside a disk, no overlaps
# ══════════════════════════════════════════════════════════════════════════════
def _scatter_in_disk(
        nodes,
        cx,
        cy,
        disk_r,
        node_r,
        rng,
        min_margin=0.45,
        max_tries=400):
    """
    Scatter nodes randomly INSIDE a safe inner region
    of the disk while preventing overlaps.

    Parameters
    ----------
    min_margin : float
        Minimum distance from outer boundary.
    """

    pos = {}
    placed = []

    # Effective usable radius
    effective_r = max(disk_r - min_margin - node_r, node_r)

    for node in nodes:

        placed_ok = False

        for _ in range(max_tries):

            # Uniform random point inside INNER disk
            r = effective_r * math.sqrt(rng.random())
            theta = rng.random() * 2 * math.pi

            x = cx + r * math.cos(theta)
            y = cy + r * math.sin(theta)

            # Prevent overlaps
            too_close = any(
                math.hypot(x - px, y - py) <= 2.0 * node_r
                for px, py in placed
            )

            if not too_close:
                placed_ok = True
                break

        # fallback
        if not placed_ok:
            angle = rng.random() * 2 * math.pi
            r = effective_r * 0.8

            x = cx + r * math.cos(angle)
            y = cy + r * math.sin(angle)

        pos[node] = (x, y)
        placed.append((x, y))

    return pos

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Draw the cluster network
# ══════════════════════════════════════════════════════════════════════════════
def _draw_cluster_graph(G: nx.DiGraph, protein: str,
                         successors: list, predecessors: list):

    buckets   = _assign_buckets(successors, predecessors)
    rng       = np.random.default_rng(seed=42)   # reproducible scatter

    NODE_R    = 0.18    # radius of each protein circle
    QUERY_R   = 0.50    # radius of the center (query) circle

    # ── Compute disk radius for each cluster proportional to node count ────────
    def _disk_radius(n):
        # Enough area to fit n circles of radius NODE_R with some breathing room
        if n == 0:   return 0.5
        if n == 1:   return NODE_R * 2.5
        return max(1.2, NODE_R * 2.5 * math.sqrt(n))

    # ── Build all positions ────────────────────────────────────────────────────
    pos       = {protein: (0.0, 0.0)}
    node_meta = {protein: ("#e74c3c", "#922b21", "query")}

    cluster_centers = {}   # ckey -> (cx, cy)
    cluster_disk_r  = {}   # ckey -> disk_r

    for ckey, cfg in CLUSTERS.items():
        nodes_in = [n for n, _ in buckets[ckey]]
        n        = len(nodes_in)
        d_r      = _disk_radius(n)

        # Orbit: push cluster center further out so disk doesn't overlap center
        orbit = cfg["radius"] + d_r * 0.9
        ang   = math.radians(cfg["angle"])
        cx    = orbit * math.cos(ang)
        cy    = orbit * math.sin(ang)

        cluster_centers[ckey] = (cx, cy)
        cluster_disk_r[ckey]  = d_r

        if n == 0:
            continue

        # Scatter nodes randomly inside the disk
        scatter = _scatter_in_disk(nodes_in, cx, cy, d_r, NODE_R, rng)
        pos.update(scatter)

        for node in nodes_in:
            node_meta[node] = (cfg["node_fc"], cfg["node_ec"], ckey)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(24, 24))
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")
    ax.set_aspect("equal")

    # ── Draw cluster background disks ─────────────────────────────────────────
    for ckey, cfg in CLUSTERS.items():
        n = len(buckets[ckey])
        if n == 0:
            continue
        cx, cy = cluster_centers[ckey]
        d_r    = cluster_disk_r[ckey]

        # Soft filled disk
        bg = plt.Circle((cx, cy), d_r,
                         facecolor=cfg["node_fc"],
                         edgecolor=cfg["node_ec"],
                         alpha=0.10,
                         linewidth=1.5,
                         linestyle="--",
                         zorder=0)
        ax.add_patch(bg)

        # Cluster label just outside the disk
        ang = math.radians(cfg["angle"])
        lx  = cx + (d_r + 0.5) * math.cos(ang)
        ly  = cy + (d_r + 0.5) * math.sin(ang)
        ax.text(lx, ly,
                f"{cfg['label']}\n({n} proteins)",
                ha="center", va="center",
                fontsize=7.5, color=cfg["node_ec"],
                fontweight="bold",
                bbox=dict(facecolor="white",
                          edgecolor=cfg["node_ec"],
                          boxstyle="round,pad=0.3",
                          alpha=0.90),
                zorder=6)

    # ── Draw edges ────────────────────────────────────────────────────────────
    succ_set = {n for n, _, _ in successors}
    pred_set = {n for n, _, _ in predecessors}

    all_edges = (
        [(protein, n, w) for n, w, _ in successors   if n in pos] +
        [(n, protein, w) for n, w, _ in predecessors if n in pos]
    )

    for u, v, w in all_edges:
        xu, yu = pos[u]
        xv, yv = pos[v]

        # Color from the neighbor's cluster
        nbr        = v if u == protein else u
        _, _, ckey = node_meta.get(nbr, ("#888", "#555", None))
        edge_col   = CLUSTERS[ckey]["edge_col"] if ckey in CLUSTERS else "#888888"

        alpha = 0.45 + 0.50 * w
        lw    = 1.0  + 2.2  * w

        # Node radii
        start_r = QUERY_R if u == protein else NODE_R
        end_r   = QUERY_R if v == protein else NODE_R

        # Direction vector
        dx = xv - xu
        dy = yv - yu
        dist = math.hypot(dx, dy)


        if dist == 0:
            continue

        ux = dx / dist
        uy = dy / dist

        # Push arrow outside node borders
        start_x = xu + ux * start_r
        start_y = yu + uy * start_r

        end_x = xv - ux * end_r
        end_y = yv - uy * end_r

        # Draw directed edge
        ax.annotate(
            "",
            xy=(end_x, end_y),
            xytext=(start_x, start_y),

            arrowprops=dict(
                arrowstyle="-|>",      # better arrow head
                color=edge_col,
                lw=lw,
                alpha=alpha,

                shrinkA=0,
                shrinkB=0,

                mutation_scale=16,     # BIGGER arrows
                connectionstyle="arc3,rad=0.08",
            ),

            zorder=2
        )

    # ── Draw nodes ────────────────────────────────────────────────────────────
    for node, (x, y) in pos.items():
        fc, ec, ckey = node_meta.get(node, ("#cccccc", "#888888", None))

        if node == protein:
            r, lw_c, fs, tc = QUERY_R, 3.0, 9, "white"
            fc, ec           = "#e74c3c", "#922b21"
            zo               = 10
        elif node in (succ_set & pred_set):
            r, lw_c, fs, tc = NODE_R * 1.15, 1.8, 5.5, "white"
            fc, ec           = "#8e44ad", "#6c3483"
            zo               = 5
        else:
            r, lw_c, fs, tc = NODE_R, 1.0, 5.0, "white"
            zo               = 4

        circ = plt.Circle((x, y), r,
                           facecolor=fc, edgecolor=ec,
                           linewidth=lw_c, zorder=zo)
        ax.add_patch(circ)

        ax.text(x, y, node,
                ha="center", va="center",
                fontsize=fs, fontweight="bold",
                color=tc, zorder=zo + 1,
                clip_on=True)

    # ── Legend ────────────────────────────────────────────────────────────────
    in_deg  = len(predecessors)
    out_deg = len(successors)
    both_n  = len(succ_set & pred_set)

    handles = [
        mpatches.Patch(facecolor="#e74c3c", edgecolor="#922b21",
                       label=f"Query: {protein}  (degree={in_deg+out_deg})"),
        mpatches.Patch(facecolor="#8e44ad", edgecolor="#6c3483",
                       label=f"Bidirectional  ({both_n} proteins)"),
    ] + [
        mpatches.Patch(facecolor=cfg["node_fc"], edgecolor=cfg["node_ec"],
                       label=f"{cfg['label']}  ({len(buckets[ckey])})")
        for ckey, cfg in CLUSTERS.items()
        if len(buckets[ckey]) > 0
    ]

    ax.legend(handles=handles,
              loc="lower center",
              bbox_to_anchor=(0.5, -0.04),
              fontsize=8.5, framealpha=0.95, ncol=4,
              title=(f"Neighbor Network — {protein}  |  "
                     f"in={in_deg}  out={out_deg}  total={in_deg+out_deg}"),
              title_fontsize=9.5,
              edgecolor="#cccccc")

    ax.set_title(
        f"Neighbor Network of  {protein}\n"
        f"Total degree: {in_deg+out_deg}  (in={in_deg}, out={out_deg})  "
        f"·  Clusters by direction & confidence tier",
        fontsize=14, fontweight="bold", pad=16, color="#2c3e50")

    # Fit axes tightly around all placed nodes
    all_x = [p[0] for p in pos.values()]
    all_y = [p[1] for p in pos.values()]
    pad   = 1.5
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - pad, max(all_y) + pad)
    ax.axis("off")

    plt.tight_layout(rect=[0, 0.04, 1, 1])

    out_path = os.path.join(FIG_DIR, "neighbors_subgraph.png")
    plt.savefig(out_path, dpi=160, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Figure saved  -> {out_path}\n")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    from load_graph import load_graph

    DATA_PATH = "data/PathLinker_2018_human-ppi-weighted-cap0_75.txt"
    G = load_graph(DATA_PATH)

    get_neighbors(G, protein="P10721")

    print("✅ neighbors.py complete!\n")