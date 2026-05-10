
"""
neighbors.py

"""

import os
import math
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Output directories ─────────────────────────────────────────────────────────
TXT_DIR = "outputs/txt"
FIG_DIR = "outputs/figures"
os.makedirs(TXT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ── Confidence thresholds ──────────────────────────────────────────────────────
HIGH_CONF = 0.75
MED_CONF  = 0.50

# ── Node / disk sizing constants ───────────────────────────────────────────────
NODE_R       = 0.18    # radius of every neighbor circle
QUERY_R      = 0.50    # radius of the center query circle
MAX_DISK_R   = 3.50    # hard cap on cluster background disk radius
MIN_DISK_R   = 0.55    # minimum disk radius (for clusters with 1–2 nodes)
ORBIT_GAP    = 0.80    # guaranteed clear gap between disk edge and center area
DISK_SPACING = 0.50    # extra gap between adjacent cluster disks

# ── Cluster config: angle (°), colors, label ──────────────────────────────────
# 'radius' is now a BASE only — actual orbit is computed dynamically below
CLUSTERS = {
    "out_high": dict(angle=120,
                     node_fc="#1a9e5c", node_ec="#0e6b3a",
                     edge_col="#27ae60",
                     label="Outgoing HIGH  (conf >= 0.75)"),
    "out_med":  dict(angle=60,
                     node_fc="#17a589", node_ec="#0e7260",
                     edge_col="#1abc9c",
                     label="Outgoing MED   (0.50-0.75)"),
    "out_low":  dict(angle=0,
                     node_fc="#d4ac0d", node_ec="#9a7d0a",
                     edge_col="#f1c40f",
                     label="Outgoing LOW   (conf < 0.50)"),
    "in_high":  dict(angle=300,
                     node_fc="#1a5276", node_ec="#0e3457",
                     edge_col="#2980b9",
                     label="Incoming HIGH  (conf >= 0.75)"),
    "in_med":   dict(angle=240,
                     node_fc="#6c3483", node_ec="#4a235a",
                     edge_col="#9b59b6",
                     label="Incoming MED   (0.50-0.75)"),
    "in_low":   dict(angle=180,
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
        print(f"  Warning: Protein '{protein}' not found in graph.")
        return

    # Collect successors (protein -> neighbor)
    successors = []
    for nbr in G.successors(protein):
        w = G[protein][nbr]["weight"]
        m = G[protein][nbr].get("method", "unknown")
        successors.append((nbr, w, m))

    # Collect predecessors (neighbor -> protein)
    predecessors = []
    for nbr in G.predecessors(protein):
        w = G[nbr][protein]["weight"]
        m = G[nbr][protein].get("method", "unknown")
        predecessors.append((nbr, w, m))

    successors   = sorted(successors,   key=lambda x: x[1], reverse=True)
    predecessors = sorted(predecessors, key=lambda x: x[1], reverse=True)

    in_deg  = G.in_degree(protein)
    out_deg = G.out_degree(protein)
    print(f"  Total degree       : {in_deg + out_deg}  "
          f"(in={in_deg}, out={out_deg})")
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
        f.write(f"Total degree   : {in_deg + out_deg}\n")
        f.write(f"  In-degree    : {in_deg}\n")
        f.write(f"  Out-degree   : {out_deg}\n")
        f.write(f"{'='*60}\n\n")

        f.write(f"OUTGOING ({len(successors)})  -- {protein} -> neighbor:\n")
        f.write(f"{'─'*60}\n")
        for nbr, w, m in successors:
            tier = ("HIGH" if w >= HIGH_CONF else
                    "MED"  if w >= MED_CONF  else "LOW")
            f.write(f"  {nbr:<12}  weight: {w:.6f}  [{tier}]  "
                    f"method: {m}\n")

        f.write(f"\nINCOMING ({len(predecessors)})  -- neighbor -> {protein}:\n")
        f.write(f"{'─'*60}\n")
        for nbr, w, m in predecessors:
            tier = ("HIGH" if w >= HIGH_CONF else
                    "MED"  if w >= MED_CONF  else "LOW")
            f.write(f"  {nbr:<12}  weight: {w:.6f}  [{tier}]  "
                    f"method: {m}\n")

    print(f"  Results saved -> {out_txt}")


# ══════════════════════════════════════════════════════════════════════════════
# HELPER — Assign neighbors into confidence-tier buckets
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
# HELPER — Compute CAPPED disk radius for a cluster of n nodes
# ══════════════════════════════════════════════════════════════════════════════
def _disk_radius(n: int) -> float:
    """
    Ideal radius grows with sqrt(n) so the disk area is proportional
    to the number of nodes. Hard-capped at MAX_DISK_R so large clusters
    don't dominate the figure.

    For very small clusters the disk is set to MIN_DISK_R so there is
    always a visible background circle.
    """
    if n == 0:
        return 0.0
    ideal = NODE_R * 2.8 * math.sqrt(n)
    return max(MIN_DISK_R, min(ideal, MAX_DISK_R))


# ══════════════════════════════════════════════════════════════════════════════
# HELPER — Scatter nodes randomly inside a disk without overlaps
# ══════════════════════════════════════════════════════════════════════════════
def _scatter_in_disk(nodes, cx, cy, disk_r, node_r, rng, max_tries=500):
    """
    Place each node at a uniformly random position inside the disk of
    radius `disk_r`, keeping nodes at least 2*node_r apart (rejection
    sampling). Falls back gracefully if no non-overlapping spot is found
    after max_tries attempts.
    """
    pos    = {}
    placed = []

    # Leave a margin so nodes don't sit on the disk boundary
    usable_r = max(disk_r - node_r - 0.10, node_r)

    for node in nodes:
        accepted = False
        for _ in range(max_tries):
            r     = usable_r * math.sqrt(rng.random())
            theta = rng.random() * 2 * math.pi
            x     = cx + r * math.cos(theta)
            y     = cy + r * math.sin(theta)

            too_close = any(
                math.hypot(x - px, y - py) < 2.0 * node_r
                for px, py in placed
            )
            if not too_close:
                accepted = True
                break

        if not accepted:
            # Fallback: place anywhere inside a slightly larger radius
            theta = rng.random() * 2 * math.pi
            x     = cx + usable_r * 0.9 * math.cos(theta)
            y     = cy + usable_r * 0.9 * math.sin(theta)

        pos[node] = (x, y)
        placed.append((x, y))

    return pos


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Draw the cluster network
# ══════════════════════════════════════════════════════════════════════════════
def _draw_cluster_graph(G: nx.DiGraph, protein: str,
                         successors: list, predecessors: list):

    buckets = _assign_buckets(successors, predecessors)
    rng     = np.random.default_rng(seed=42)   # reproducible layout

    # ── Step 1: compute each cluster's disk radius ────────────────────────────
    disk_r = {}
    for ckey, bucket in buckets.items():
        disk_r[ckey] = _disk_radius(len(bucket))

    # ── Step 2: compute orbit distance for each cluster ───────────────────────
    # Orbit = distance from origin to cluster CENTER.
    # We want:   orbit > QUERY_R + ORBIT_GAP + disk_r   (no overlap with center)
    # We also want each cluster disk NOT to overlap its two angular neighbors.
    #
    # For 6 equally-spaced clusters (60 deg apart), the chord between adjacent
    # centers at orbit distance R is:  chord = 2 * R * sin(30 deg) = R
    # So the minimum orbit so adjacent disks don't touch is:
    #   R >= disk_r[A] + disk_r[B] + DISK_SPACING
    #
    # We compute a single orbit radius that satisfies ALL constraints.

    ckeys  = list(CLUSTERS.keys())
    n_clus = len(ckeys)

    # Constraint 1: clear of center
    orbit_min_center = [
        QUERY_R + ORBIT_GAP + disk_r[ck]
        for ck in ckeys
    ]

    # Constraint 2: clear of angular neighbors (clusters are 60 deg apart)
    # chord = orbit * 2 * sin(pi/6) = orbit * 1.0  =>  orbit >= dr_A + dr_B + gap
    orbit_min_neighbors = []
    for i, ck in enumerate(ckeys):
        left_ck  = ckeys[(i - 1) % n_clus]
        right_ck = ckeys[(i + 1) % n_clus]
        needed   = (disk_r[ck] + disk_r[left_ck]  + DISK_SPACING,
                    disk_r[ck] + disk_r[right_ck] + DISK_SPACING)
        orbit_min_neighbors.append(max(needed))

    # Use one shared orbit radius = max of all constraints (keeps layout tidy)
    orbit = max(max(orbit_min_center), max(orbit_min_neighbors))

    # ── Step 3: cluster centers ───────────────────────────────────────────────
    cluster_centers = {}
    for ckey, cfg in CLUSTERS.items():
        ang = math.radians(cfg["angle"])
        cluster_centers[ckey] = (orbit * math.cos(ang),
                                  orbit * math.sin(ang))

    # ── Step 4: scatter nodes inside each cluster disk ────────────────────────
    pos       = {protein: (0.0, 0.0)}
    node_meta = {protein: ("#e74c3c", "#922b21", "query")}

    for ckey, cfg in CLUSTERS.items():
        nodes_in = [n for n, _ in buckets[ckey]]
        if not nodes_in:
            continue
        cx, cy = cluster_centers[ckey]
        scatter = _scatter_in_disk(nodes_in, cx, cy,
                                   disk_r[ckey], NODE_R, rng)
        pos.update(scatter)
        for node in nodes_in:
            node_meta[node] = (cfg["node_fc"], cfg["node_ec"], ckey)

    # ── Figure setup ──────────────────────────────────────────────────────────
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
        dr     = disk_r[ckey]

        bg = plt.Circle((cx, cy), dr,
                         facecolor=cfg["node_fc"],
                         edgecolor=cfg["node_ec"],
                         alpha=0.12,
                         linewidth=1.8,
                         linestyle="--",
                         zorder=0)
        ax.add_patch(bg)

        # Label just outside the disk, pushed further along the same angle
        ang = math.radians(cfg["angle"])
        lx  = cx + (dr + 0.65) * math.cos(ang)
        ly  = cy + (dr + 0.65) * math.sin(ang)
        ax.text(lx, ly,
                f"{cfg['label']}\n({n} proteins)",
                ha="center", va="center",
                fontsize=8, color=cfg["node_ec"],
                fontweight="bold",
                bbox=dict(facecolor="white",
                          edgecolor=cfg["node_ec"],
                          boxstyle="round,pad=0.35",
                          alpha=0.92),
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

        nbr        = v if u == protein else u
        _, _, ckey = node_meta.get(nbr, ("#888", "#555", None))
        edge_col   = (CLUSTERS[ckey]["edge_col"]
                      if ckey in CLUSTERS else "#888888")

        alpha = 0.40 + 0.55 * w
        lw    = 0.8  + 2.0  * w

        dx   = xv - xu
        dy   = yv - yu
        dist = math.hypot(dx, dy)
        if dist == 0:
            continue

        ux = dx / dist
        uy = dy / dist

        start_r = QUERY_R if u == protein else NODE_R
        end_r   = QUERY_R if v == protein else NODE_R

        ax.annotate(
            "",
            xy    =(xu + ux * (dist - end_r),
                    yu + uy * (dist - end_r)),
            xytext=(xu + ux * start_r,
                    yu + uy * start_r),
            arrowprops=dict(
                arrowstyle="-|>",
                color=edge_col,
                lw=lw,
                alpha=alpha,
                shrinkA=0,
                shrinkB=0,
                mutation_scale=14,
                connectionstyle="arc3,rad=0.06",
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

        ax.add_patch(plt.Circle((x, y), r,
                                facecolor=fc, edgecolor=ec,
                                linewidth=lw_c, zorder=zo))
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
                       label=f"Query: {protein}  "
                             f"(degree={in_deg + out_deg})"),
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
              bbox_to_anchor=(0.5, -0.03),
              fontsize=9, framealpha=0.95, ncol=4,
              title=(f"Neighbor Network -- {protein}  |  "
                     f"in={in_deg}  out={out_deg}  "
                     f"total={in_deg + out_deg}"),
              title_fontsize=10,
              edgecolor="#cccccc")

    ax.set_title(
        f"Neighbor Network of  {protein}\n"
        f"Total degree: {in_deg + out_deg}  "
        f"(in={in_deg}, out={out_deg})  --  "
        f"Clusters by direction & confidence tier",
        fontsize=14, fontweight="bold", pad=16, color="#2c3e50")

    all_x = [p[0] for p in pos.values()]
    all_y = [p[1] for p in pos.values()]
    pad   = MAX_DISK_R + 1.5
    ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax.set_ylim(min(all_y) - pad, max(all_y) + pad)
    ax.axis("off")

    plt.tight_layout(rect=[0, 0.04, 1, 1])

    out_path = os.path.join(FIG_DIR, "neighbors_subgraph.png")
    plt.savefig(out_path, dpi=160, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Figure saved  -> {out_path}\n")

