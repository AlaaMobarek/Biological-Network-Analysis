"""
pathway_viz.py  —  KEGG Pathway Sub-Network Visualizer
=======================================================
Fixes vs. old version
  1.  gene mapping  → runs ONLY on pathway proteins (not all 17 k nodes)
  2.  functional groups / colors  → fetched from KEGG dynamically
  3.  layer layout  → computed automatically via BFS from source
  4.  zero manual UniProt lists  in this file
  5.  nodes → filled circles  (label centred inside)
"""

import os
import time
import requests
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from collections import defaultdict
import mygene


# Output dirs
OUT_FIGS = os.path.join("outputs", "figures")
OUT_TXTS = os.path.join("outputs", "txt")
os.makedirs(OUT_FIGS, exist_ok=True)
os.makedirs(OUT_TXTS, exist_ok=True)

# Defaults
DEFAULT_PATHWAY_ID    = "hsa04151"   # PI3K-AKT
DEFAULT_SOURCE_ID     = "P04637"     # TP53
DEFAULT_TARGET_ID     = "P31749"     # AKT1
DEFAULT_CONF_THRESH   = 0.75
MAX_EDGES             = 200

# ── A small but visible qualitative palette (12 distinct colours) ──────────────
_PALETTE = [
    "#a8d8a8", "#f9c784", "#a0c4ff", "#ff9aa2", "#2196F3",
    "#b5ead7", "#c9b1ff", "#ffd6e7", "#ffe0ac", "#ff6b6b",
    "#d4a5a5", "#d0d0d0",
]

# ── Protein grouping & colors (for draw_top_hubs) ──────────────────────────────
PROTEIN_GROUP = {}
GROUP_COLORS = {
    "Kinase": "#FF6B6B",
    "TF": "#4ECDC4",
    "Adapter": "#95E1D3",
    "Other": "#CCCCCC",
}


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 0  —  Bulk UniProt ID → Gene Name mapping
# ══════════════════════════════════════════════════════════════════════════════

def build_gene_dict(uniprot_ids: list) -> dict:
    mg      = mygene.MyGeneInfo()
    results = mg.querymany(uniprot_ids, scopes="uniprot",
                           fields="symbol", species="human", verbose=False)
    mapping = {uid: uid for uid in uniprot_ids}
    for r in results:
        if "symbol" in r:
            mapping[r["query"]] = r["symbol"]
    found = sum(1 for k, v in mapping.items() if v != k)
    print(f"[Mapping] Done -- {found}/{len(uniprot_ids)} IDs mapped OK")
    return mapping


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1  —  Fetch pathway proteins from KEGG
# ══════════════════════════════════════════════════════════════════════════════

def fetch_kegg_uniprot_ids(pathway_id: str) -> set:
    print(f"\n[KEGG] Fetching proteins for  {pathway_id}  …")
    url_a = f"https://rest.kegg.jp/link/hsa/{pathway_id}"
    try:
        r_a = requests.get(url_a, timeout=20)
        r_a.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"[KEGG] Step-A ERROR — {e}")
        return set()

    pathway_hsa_genes = set()
    for line in r_a.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2:
            pathway_hsa_genes.add(parts[1].strip())

    if not pathway_hsa_genes:
        print(f"[KEGG] No genes returned for pathway {pathway_id}")
        return set()
    print(f"[KEGG] Step-A: {len(pathway_hsa_genes)} hsa gene IDs in pathway")

    uniprot_ids = set()
    gene_list   = sorted(pathway_hsa_genes)
    BATCH       = 100
    for start in range(0, len(gene_list), BATCH):
        batch = gene_list[start: start + BATCH]
        url_b = f"https://rest.kegg.jp/conv/uniprot/{'+'.join(batch)}"
        try:
            r_b = requests.get(url_b, timeout=20)
            r_b.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"[KEGG] Step-B batch ERROR — {e}")
            continue
        for line in r_b.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                uniprot_ids.add(parts[1].strip().split(":")[-1])
        time.sleep(0.1)

    print(f"[KEGG] Step-B: {len(uniprot_ids)} UniProt IDs resolved")
    return uniprot_ids


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1b  —  Fetch KEGG modules/groups dynamically
# ══════════════════════════════════════════════════════════════════════════════

def fetch_kegg_modules(pathway_id: str, valid_proteins: set) -> dict:
    try:
        r = requests.get(f"https://rest.kegg.jp/link/module/{pathway_id}",
                         timeout=15)
        r.raise_for_status()
    except Exception:
        return {}

    gene_to_module = {}
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2:
            gene_to_module[parts[0].strip()] = parts[1].strip().split(":")[-1]

    if not gene_to_module:
        return {}

    try:
        r2 = requests.get(f"https://rest.kegg.jp/link/hsa/{pathway_id}",
                          timeout=15)
        r2.raise_for_status()
        pathway_genes = [ln.split("\t")[1].strip()
                         for ln in r2.text.strip().split("\n")
                         if len(ln.split("\t")) == 2]
    except Exception:
        return {}

    kegg_to_uniprot = {}
    BATCH = 100
    for start in range(0, len(pathway_genes), BATCH):
        batch = pathway_genes[start: start + BATCH]
        try:
            r3 = requests.get(
                f"https://rest.kegg.jp/conv/uniprot/{'+'.join(batch)}",
                timeout=20)
            r3.raise_for_status()
            for ln in r3.text.strip().split("\n"):
                p = ln.split("\t")
                if len(p) == 2:
                    kegg_to_uniprot[p[0].strip()] = p[1].strip().split(":")[-1]
        except Exception:
            pass
        time.sleep(0.1)

    module_ids  = list(set(gene_to_module.values()))[:30]
    module_name = {}
    for mid in module_ids:
        try:
            r4 = requests.get(f"https://rest.kegg.jp/get/{mid}", timeout=10)
            for ln in r4.text.split("\n"):
                if ln.startswith("NAME"):
                    module_name[mid] = ln.split(None, 1)[1].strip()[:25]
                    break
        except Exception:
            pass
        time.sleep(0.05)

    protein_module = {}
    for kg, uid in kegg_to_uniprot.items():
        if uid in valid_proteins and kg in gene_to_module:
            mid = gene_to_module[kg]
            protein_module[uid] = module_name.get(mid, mid)

    print(f"[Modules] {len(protein_module)} proteins -> "
          f"{len(set(protein_module.values()))} modules")
    return protein_module


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2  —  Filter to PathLinker dataset
# ══════════════════════════════════════════════════════════════════════════════

def filter_to_dataset(kegg_proteins: set, G_full: nx.DiGraph) -> set:
    valid   = kegg_proteins & set(G_full.nodes())
    missing = kegg_proteins - valid
    print(f"\n[Filter] KEGG proteins : {len(kegg_proteins)}")
    print(f"[Filter] In dataset    : {len(valid)}")
    print(f"[Filter] Missing       : {len(missing)}")
    return valid


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3  —  Build confidence-filtered subgraph
# ══════════════════════════════════════════════════════════════════════════════

def build_pathway_subgraph(G_full: nx.DiGraph,
                            pathway_proteins: set,
                            conf_threshold: float) -> nx.DiGraph:
    G_sub = nx.DiGraph()
    for u, v, data in G_full.edges(data=True):
        if u in pathway_proteins and v in pathway_proteins:
            w = data.get("weight", 0)
            if w >= conf_threshold:
                G_sub.add_edge(u, v, weight=w, cost=round(1.0 - w, 6))
    print(f"\n[Subgraph] Threshold {conf_threshold}  ->  "
          f"{G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")
    return G_sub


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4  —  Helpers
# ══════════════════════════════════════════════════════════════════════════════

def prune_to_top_edges(G: nx.DiGraph,
                        max_edges: int,
                        must_keep: set = None) -> nx.DiGraph:
    must_keep = must_keep or set()
    all_e  = sorted(G.edges(data=True),
                    key=lambda e: e[2].get("weight", 0), reverse=True)
    pinned = [(u, v, d) for u, v, d in all_e
              if u in must_keep or v in must_keep]
    rest   = [(u, v, d) for u, v, d in all_e
              if u not in must_keep and v not in must_keep]
    budget = max(max_edges, len(pinned))
    chosen = pinned + rest[: budget - len(pinned)]

    G2 = nx.DiGraph()
    G2.add_nodes_from(G.nodes())
    G2.add_edges_from((u, v, d) for u, v, d in chosen)
    isolated = [n for n in list(G2.nodes())
                if G2.degree(n) == 0 and n not in must_keep]
    G2.remove_nodes_from(isolated)
    print(f"[Prune]  {G.number_of_edges()} -> {G2.number_of_edges()} edges  "
          f"({len(isolated)} isolated removed)")
    return G2


def _bfs_layers(G: nx.DiGraph, source_id: str, nodes: set) -> dict:
    layers = {}
    G_und  = G.to_undirected()
    source = source_id if source_id in G_und.nodes() else next(iter(G_und.nodes()))
    queue  = [source]
    layers[source] = 0
    while queue:
        current = queue.pop(0)
        for nb in G_und.neighbors(current):
            if nb not in layers:
                layers[nb] = layers[current] + 1
                queue.append(nb)
    for n in nodes:
        if n not in layers:
            layers[n] = 3
    return layers


def auto_hierarchical_layout(G: nx.DiGraph,
                              source_id: str,
                              gene_dict: dict) -> dict:
    layers     = _bfs_layers(G, source_id, set(G.nodes()))
    max_layer  = max(layers.values()) if layers else 1
    by_layer   = defaultdict(list)
    for node, lyr in layers.items():
        by_layer[lyr].append(node)
    pos = {}
    for lyr, nodes in by_layer.items():
        nodes = sorted(nodes, key=lambda n: gene_dict.get(n, n))
        n     = len(nodes)
        y     = 1.0 - lyr / max(max_layer, 1)
        for i, node in enumerate(nodes):
            pos[node] = ((i + 0.5) / n, y)
    return pos


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 5  —  Draw   (circles only — no rectangles)
# ══════════════════════════════════════════════════════════════════════════════

def draw_pathway(G: nx.DiGraph,
                 gene_dict: dict,
                 protein_module: dict,
                 source_id: str,
                 target_id: str,
                 pathway_id: str,
                 pathway_name: str,
                 conf_threshold: float,
                 out_path: str):

    if G.number_of_nodes() == 0:
        print("[Draw] Subgraph empty — nothing to draw.")
        return

    # ── module colour map ──────────────────────────────────────────────────
    all_modules = sorted(set(protein_module.get(n, "Other") for n in G.nodes()))
    mod_color   = {m: _PALETTE[i % len(_PALETTE)]
                   for i, m in enumerate(all_modules)}

    # ── layout (unchanged — BFS hierarchical) ─────────────────────────────
    pos = auto_hierarchical_layout(G, source_id, gene_dict)

    fig, ax = plt.subplots(figsize=(26, 18))
    ax.set_facecolor("#f8f9fa")
    fig.patch.set_facecolor("#f8f9fa")

    # ── faint horizontal bands per depth level ─────────────────────────────
    y_values = sorted(set(y for _, y in pos.values()), reverse=True)
    band_bg  = ["#eaf4fb", "#eafaf1", "#fef9e7", "#e8f8f5",
                 "#f4ecf7", "#fdfefe", "#fdedec"]
    for i, yc in enumerate(y_values):
        dy = (1.0 / max(len(y_values) - 1, 1)) * 0.46
        ax.axhspan(yc - dy, yc + dy,
                   color=band_bg[i % len(band_bg)], alpha=0.35, zorder=0)

    # ── node colours & sizes ──────────────────────────────────────────────
    node_list = list(G.nodes())
    node_colors, node_sizes, font_colors = [], [], []
    for n in node_list:
        if n == source_id:
            node_colors.append("#e74c3c")
            node_sizes.append(2800)
            font_colors.append("white")
        elif n == target_id:
            node_colors.append("#27ae60")
            node_sizes.append(2800)
            font_colors.append("white")
        else:
            node_colors.append(mod_color.get(protein_module.get(n, "Other"), "#d0d0d0"))
            node_sizes.append(1400)
            font_colors.append("#111111")

    # ── edges ─────────────────────────────────────────────────────────────
    edge_weights = [G[u][v].get("weight", 0) for u, v in G.edges()]
    if edge_weights:
        wmin, wmax = min(edge_weights), max(edge_weights)
        wr         = wmax - wmin if wmax > wmin else 1.0
        for (u, v), w in zip(G.edges(), edge_weights):
            nw  = (w - wmin) / wr
            col = (0.15 + 0.65 * nw, 0.25 - 0.15 * nw, 0.85 - 0.65 * nw, 0.80)
            lw  = 1.0 + 2.0 * nw
            ax.annotate(
                "", xy=pos[v], xytext=pos[u],
                arrowprops=dict(
                    arrowstyle="-|>", color=col, lw=lw,
                    mutation_scale=14,
                    connectionstyle="arc3,rad=0.06",
                    shrinkA=22, shrinkB=22,   # ← يبعد السهم عن مركز الدايرة
                ),
                zorder=2,
            )

    # ── draw circles (cluster nodes) ──────────────────────────────────────
    reg_nodes  = [n for n in node_list if n not in {source_id, target_id}]
    reg_colors = [mod_color.get(protein_module.get(n, "Other"), "#d0d0d0")
                  for n in reg_nodes]

    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        nodelist=reg_nodes,
        node_color=reg_colors,
        node_size=1400,
        node_shape="o",
        linewidths=1.2,
        edgecolors="#555555",
    )

    # ── source circle (bigger, red) ───────────────────────────────────────
    if source_id in G.nodes():
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=[source_id],
            node_color="#e74c3c",
            node_size=2800,
            node_shape="o",
            linewidths=2.0,
            edgecolors="#922b21",
        )

    # ── target circle (bigger, green) ─────────────────────────────────────
    if target_id in G.nodes():
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            nodelist=[target_id],
            node_color="#27ae60",
            node_size=2800,
            node_shape="o",
            linewidths=2.0,
            edgecolors="#1a6b3c",
        )

    # ── labels inside circles ─────────────────────────────────────────────
    labels = {n: gene_dict.get(n, n) for n in G.nodes()}

    # cluster nodes — small dark label
    nx.draw_networkx_labels(
        G, pos,
        {n: labels[n] for n in reg_nodes},
        ax=ax,
        font_size=6.5,
        font_weight="bold",
        font_color="#111111",
        font_family="DejaVu Sans Mono",
    )

    # source & target — white label, slightly bigger
    spec = {}
    if source_id in G.nodes(): spec[source_id] = labels[source_id]
    if target_id in G.nodes(): spec[target_id] = labels[target_id]
    if spec:
        nx.draw_networkx_labels(
            G, pos, spec, ax=ax,
            font_size=9,
            font_weight="bold",
            font_color="white",
            font_family="DejaVu Sans Mono",
        )

    # ── colorbar ──────────────────────────────────────────────────────────
    cmap = LinearSegmentedColormap.from_list(
        "conf", [(0.15, 0.25, 0.85), (0.80, 0.10, 0.20)], N=256)
    sm   = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=conf_threshold, vmax=1.0))
    sm.set_array([])
    cb_ax = fig.add_axes([0.92, 0.15, 0.012, 0.35])
    cb = fig.colorbar(sm, cax=cb_ax)
    cb.set_label("Interaction\nConfidence", fontsize=8, labelpad=6)
    cb.ax.tick_params(labelsize=7)

    # ── legend ────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(facecolor=mod_color[m], edgecolor="#555",
                       linewidth=0.6, label=m)
        for m in all_modules
    ] + [
        mpatches.Patch(facecolor="#e74c3c", edgecolor="#922b21", linewidth=1.2,
                       label=f"{gene_dict.get(source_id, source_id)} – Source"),
        mpatches.Patch(facecolor="#27ae60", edgecolor="#1a6b3c", linewidth=1.2,
                       label=f"{gene_dict.get(target_id, target_id)} – Target"),
    ]
    ax.legend(handles=legend_patches, loc="lower center",
              bbox_to_anchor=(0.45, -0.10),
              fontsize=6.5, framealpha=0.9, ncol=5,
              title="KEGG Module / Group", title_fontsize=7, edgecolor="#aaa")

    # ── title ─────────────────────────────────────────────────────────────
    ax.set_title(
        f"PPI Sub-Network: {pathway_name}  ({pathway_id})\n"
        f"{G.number_of_nodes()} proteins  ·  {G.number_of_edges()} interactions"
        f"  ·  Confidence ≥ {conf_threshold}  ·  BFS-hierarchical layout",
        fontsize=13, fontweight="bold", pad=16)

    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.07, 1.07)
    ax.axis("off")
    plt.tight_layout(rect=[0, 0.07, 0.91, 1])
    plt.savefig(out_path, dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[Draw] Saved -> {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 6  —  Save text report
# ══════════════════════════════════════════════════════════════════════════════

def save_pathway_report(valid_proteins: set, gene_dict: dict,
                         protein_module: dict,
                         source_id: str, target_id: str,
                         pathway_id: str, pathway_name: str,
                         out_path: str):
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"Pathway  : {pathway_name} ({pathway_id})\n")
        fh.write(f"Proteins : {len(valid_proteins)}\n\n")
        fh.write(f"{'UniProt ID':<14} {'Gene Name':<14} {'Module':<26} {'Note'}\n")
        fh.write("─" * 70 + "\n")
        for pid in sorted(valid_proteins):
            tag = (" ← Source" if pid == source_id else
                   " ← Target" if pid == target_id else "")
            fh.write(f"{pid:<14} {gene_dict.get(pid, pid):<14} "
                     f"{protein_module.get(pid, 'Other'):<26}{tag}\n")
    print(f"[Report] Saved -> {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
#  PUBLIC ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def run_pathway_visualization(G_full: nx.DiGraph,
                               pathway_id: str         = DEFAULT_PATHWAY_ID,
                               source_id: str          = DEFAULT_SOURCE_ID,
                               target_id: str          = DEFAULT_TARGET_ID,
                               conf_threshold: float   = DEFAULT_CONF_THRESH,
                               pathway_name: str       = None,
                               gene_dict: dict         = None):

    if pathway_name is None:
        try:
            r = requests.get(f"https://rest.kegg.jp/get/{pathway_id}", timeout=10)
            for line in r.text.split("\n"):
                if line.startswith("NAME"):
                    pathway_name = line.split(None, 1)[1].strip()
                    break
        except Exception:
            pass
    pathway_name = pathway_name or pathway_id

    print(f"\n{'='*60}")
    print(f"  Pathway : {pathway_name} ({pathway_id})")
    print(f"  Source  : {source_id}")
    print(f"  Target  : {target_id}")
    print(f"  Conf >= : {conf_threshold}")
    print(f"{'='*60}")

    kegg_proteins = fetch_kegg_uniprot_ids(pathway_id)
    if not kegg_proteins:
        print("[ERROR] No proteins from KEGG.")
        return

    valid_proteins = filter_to_dataset(kegg_proteins, G_full)
    if not valid_proteins:
        print("[ERROR] No pathway proteins found in the dataset.")
        return

    if gene_dict is None:
        gene_dict = build_gene_dict(sorted(valid_proteins))
    else:
        gene_dict = {uid: gene_dict.get(uid, uid) for uid in valid_proteins}

    protein_module = fetch_kegg_modules(pathway_id, valid_proteins)
    if not protein_module:
        print("[Modules] KEGG module info unavailable -- all nodes -> 'Other'")

    G_sub = build_pathway_subgraph(G_full, valid_proteins, conf_threshold)
    if G_sub.number_of_nodes() == 0:
        print(f"[ERROR] Empty subgraph — try lowering conf_threshold "
              f"(current: {conf_threshold})")
        return

    if G_sub.number_of_edges() > MAX_EDGES:
        G_draw = prune_to_top_edges(G_sub, MAX_EDGES,
                                     must_keep={source_id, target_id})
    else:
        G_draw = G_sub

    out_txt = os.path.join(OUT_TXTS, f"{pathway_id}_proteins.txt")
    save_pathway_report(valid_proteins, gene_dict, protein_module,
                         source_id, target_id,
                         pathway_id, pathway_name, out_txt)

    out_fig = os.path.join(OUT_FIGS, f"{pathway_id}_network.png")
    draw_pathway(G_draw, gene_dict, protein_module,
                 source_id, target_id,
                 pathway_id, pathway_name,
                 conf_threshold, out_fig)

    print(f"\n[DONE]  Figure  -> {out_fig}")
    print(f"        Report  -> {out_txt}")


# ══════════════════════════════════════════════════════════════════════════════
#  TOP HUBS VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════
def draw_top_hubs(
        G: nx.DiGraph,
        top_n: int = 20,
        out_path: str = "outputs/figures/top_hubs.png"):
    """
    Draw a graph of the Top-N hub proteins
    (highest total degree in the interactome).

    Features
    --------
    - Light theme
    - Node size proportional to total degree
    - Node color based on functional group
    - Edge color reflects interaction confidence
    - Degree shown below each node
    - Graph statistics shown in title
    """

    print(f"\n[hubs] Selecting Top-{top_n} hub proteins...")

    # ─────────────────────────────────────────────────────────────
    # 1) Compute total degree for all proteins
    # ─────────────────────────────────────────────────────────────
    degree_dict = dict(G.degree())

    # Sort descending by degree
    top_hubs = sorted(
        degree_dict.items(),
        key=lambda x: x[1],
        reverse=True
    )[:top_n]

    hub_nodes = [node for node, deg in top_hubs]

    print(f"[hubs] Top hubs selected:")

    for i, (node, deg) in enumerate(top_hubs, 1):

        gene = G.nodes[node].get("gene_name", node)

        print(
            f"   {i:>2}. "
            f"{gene:<10} ({node})   "
            f"degree={deg}"
        )

    # ─────────────────────────────────────────────────────────────
    # 2) Create subgraph
    # ─────────────────────────────────────────────────────────────
    H = G.subgraph(hub_nodes).copy()

    print(
        f"[hubs] Subgraph: "
        f"{H.number_of_nodes()} nodes, "
        f"{H.number_of_edges()} edges"
    )

    # ─────────────────────────────────────────────────────────────
    # 3) Layout
    # ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 12))

    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f8f9fa")

    pos = nx.spring_layout(
        H,
        seed=42,
        k=1.2
    )

    # ─────────────────────────────────────────────────────────────
    # 4) Node sizes proportional to degree
    # ─────────────────────────────────────────────────────────────
    max_deg = max(degree_dict[n] for n in H.nodes())

    node_sizes = [
        500 + 3500 * (degree_dict[n] / max_deg)
        for n in H.nodes()
    ]

    # ─────────────────────────────────────────────────────────────
    # 5) Node colors by functional group
    # ─────────────────────────────────────────────────────────────
    node_colors = []

    for n in H.nodes():

        group = PROTEIN_GROUP.get(n, "Other")

        node_colors.append(
            GROUP_COLORS.get(group, "#cccccc")
        )

    # ─────────────────────────────────────────────────────────────
    # 6) Edge confidence weights
    # ─────────────────────────────────────────────────────────────
    edge_weights = [
        H[u][v].get("weight", 0.5)
        for u, v in H.edges()
    ]

    edge_widths = [
        1.0 + 4.0 * w
        for w in edge_weights
    ]

    edge_colors = [
        plt.cm.coolwarm(w)
        for w in edge_weights
    ]

    # ─────────────────────────────────────────────────────────────
    # 7) Draw edges
    # ─────────────────────────────────────────────────────────────
    for (u, v), col, wid in zip(
            H.edges(),
            edge_colors,
            edge_widths):

        nx.draw_networkx_edges(
            H,
            pos,
            edgelist=[(u, v)],
            edge_color=[col],
            width=wid,
            alpha=0.75,
            arrows=True,
            arrowsize=18,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.08",
            ax=ax
        )

    # ─────────────────────────────────────────────────────────────
    # 8) Draw nodes
    # ─────────────────────────────────────────────────────────────
    nx.draw_networkx_nodes(
        H,
        pos,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="black",
        linewidths=1.2,
        alpha=0.95,
        ax=ax
    )

    # ─────────────────────────────────────────────────────────────
    # 9) Labels (gene names)
    # ─────────────────────────────────────────────────────────────
    labels = {
        n: G.nodes[n].get("gene_name", n)
        for n in H.nodes()
    }

    nx.draw_networkx_labels(
        H,
        pos,
        labels,
        font_size=9,
        font_weight="bold",
        font_color="black",
        ax=ax
    )

    # ─────────────────────────────────────────────────────────────
    # 10) Degree annotation BELOW each node
    # ─────────────────────────────────────────────────────────────
    for node, (x, y) in pos.items():

        ax.text(
            x,
            y - 0.07,
            f"deg={degree_dict[node]}",
            ha="center",
            va="top",
            fontsize=6,
            color="#444444",
            bbox=dict(
                facecolor="white",
                edgecolor="#cccccc",
                alpha=0.85,
                boxstyle="round,pad=0.2"
            ),
            zorder=5
        )

    # ─────────────────────────────────────────────────────────────
    # 11) Functional group legend
    # ─────────────────────────────────────────────────────────────
    shown_groups = sorted({
        PROTEIN_GROUP.get(n, "Other")
        for n in H.nodes()
    })

    legend_patches = [
        mpatches.Patch(
            facecolor=GROUP_COLORS.get(g, "#cccccc"),
            edgecolor="#555555",
            linewidth=0.8,
            label=g.replace("_", " ")
        )
        for g in shown_groups
    ]

    ax.legend(
        handles=legend_patches,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        fontsize=8,
        ncol=4,
        framealpha=0.95,
        facecolor="white",
        edgecolor="#999999",
        title="Functional Group",
        title_fontsize=9
    )

    # ─────────────────────────────────────────────────────────────
    # 12) Colorbar for edge confidence
    # ─────────────────────────────────────────────────────────────
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.coolwarm,
        norm=plt.Normalize(vmin=0, vmax=1)
    )

    sm.set_array([])

    cbar = plt.colorbar(
        sm,
        ax=ax,
        fraction=0.025,
        pad=0.02
    )

    cbar.set_label(
        "Interaction Confidence",
        fontsize=9
    )

    cbar.ax.tick_params(labelsize=7)

    # ─────────────────────────────────────────────────────────────
    # 13) Graph statistics in title
    # ─────────────────────────────────────────────────────────────
    density = nx.density(H)

    plt.title(
        f"Top-{top_n} Hub Proteins in Interactome\n"
        f"Nodes: {H.number_of_nodes()}   |   "
        f"Interactions: {H.number_of_edges()}   |   "
        f"Density: {density:.3f}\n"
        f"Node Size = Degree   |   "
        f"Edge Color = Confidence",
        fontsize=15,
        fontweight="bold",
        color="black",
        pad=16
    )

    # ─────────────────────────────────────────────────────────────
    # 14) Final formatting
    # ─────────────────────────────────────────────────────────────
    plt.axis("off")

    plt.tight_layout(
        rect=[0, 0.05, 1, 1]
    )

    # ─────────────────────────────────────────────────────────────
    # 15) Save figure
    # ─────────────────────────────────────────────────────────────
    plt.savefig(
        out_path,
        dpi=200,
        bbox_inches="tight",
        facecolor=fig.get_facecolor()
    )

    plt.close()

    print(f"[hubs] Figure saved -> {out_path}")


    
# def draw_top_hubs(G: nx.DiGraph,
#                   top_n: int = 20,
#                   out_path: str = None):
#     """
#     Draw a graph of the Top-N hub proteins
#     (highest total degree in the interactome).

#     Parameters
#     ----------
#     G : nx.DiGraph
#         Full interactome graph

#     top_n : int
#         Number of hub proteins to visualize

#     out_path : str
#         Output figure path (default: outputs/figures/top_hubs.png)
#     """
#     if out_path is None:
#         out_path = os.path.join(OUT_FIGS, "top_hubs.png")

#     print(f"\n[hubs] Selecting Top-{top_n} hub proteins...")

#     # ── Compute total degree for every node ────────────────────────────────
#     degree_dict = dict(G.degree())

#     # ── Sort descending by degree ──────────────────────────────────────────
#     top_hubs = sorted(
#         degree_dict.items(),
#         key=lambda x: x[1],
#         reverse=True
#     )[:top_n]

#     hub_nodes = [node for node, deg in top_hubs]

#     print(f"[hubs] Top hubs selected:")
#     for i, (node, deg) in enumerate(top_hubs, 1):
#         gene = G.nodes[node].get("gene_name", node)
#         print(f"   {i:>2}. {gene:<10} ({node})  degree={deg}")

#     # ── Create subgraph ────────────────────────────────────────────────────
#     H = G.subgraph(hub_nodes).copy()

#     print(f"[hubs] Subgraph: "
#           f"{H.number_of_nodes()} nodes, "
#           f"{H.number_of_edges()} edges")

#     # ── Layout ─────────────────────────────────────────────────────────────
#     fig, ax = plt.subplots(figsize=(14, 10))
#     pos = nx.spring_layout(H, seed=42, k=1.2)

#     # ── Node sizes proportional to degree ─────────────────────────────────
#     node_sizes = [
#         degree_dict[n] * 2
#         for n in H.nodes()
#     ]

#     # ── Node colors ────────────────────────────────────────────────────────
#     node_colors = []

#     for n in H.nodes():
#         group = PROTEIN_GROUP.get(n, "Other")
#         node_colors.append(
#             GROUP_COLORS.get(group, "#cccccc")
#         )

#     # ── Edge weights for visualization ────────────────────────────────────
#     edge_weights = [
#         H[u][v].get("weight", 0.5)
#         for u, v in H.edges()
#     ]

#     # ── Draw edges ─────────────────────────────────────────────────────────
#     nx.draw_networkx_edges(
#         H,
#         pos,
#         ax=ax,
#         edge_color=edge_weights,
#         edge_cmap=plt.cm.coolwarm,
#         width=1.5,
#         alpha=0.7,
#         arrows=True,
#         arrowsize=14
#     )

#     # ── Draw nodes ─────────────────────────────────────────────────────────
#     nx.draw_networkx_nodes(
#         H,
#         pos,
#         ax=ax,
#         node_size=node_sizes,
#         node_color=node_colors,
#         edgecolors="black",
#         linewidths=1.2,
#         alpha=0.95
#     )

#     # ── Labels ─────────────────────────────────────────────────────────────
#     labels = {
#         n: G.nodes[n].get("gene_name", n)
#         for n in H.nodes()
#     }

#     nx.draw_networkx_labels(
#         H,
#         pos,
#         labels,
#         ax=ax,
#         font_size=9,
#         font_weight="bold"
#     )

#     # ── Title ──────────────────────────────────────────────────────────────
#     ax.set_title(
#         f"Top-{top_n} Hub Proteins in Interactome",
#         fontsize=16,
#         fontweight="bold"
#     )

#     ax.axis("off")
#     plt.tight_layout()

#     # ── Save ───────────────────────────────────────────────────────────────
#     os.makedirs(os.path.dirname(out_path), exist_ok=True)
#     plt.savefig(out_path, dpi=200, bbox_inches="tight")
#     plt.close(fig)

#     print(f"[hubs] Figure saved -> {out_path}")

