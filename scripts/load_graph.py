# """
# load_graph.py
# =============
# Two things in one file:

# 1. KEGG-style PI3K-AKT sub-network visualization
#    - Colored rectangular boxes grouped by functional category
#    - Directed edges with visible arrowheads (offset to box edges)
#    - Hierarchical layer layout, confidence colorbar, legend below graph
#    - Top-N edge pruning with isolated-node removal

# 2. Full PathLinker interactome DiGraph builder  (build_full_graph)
#    - Loads every edge in the file, no protein filter, no confidence cutoff
#    - Stores weight / cost / neg_log_weight on every edge
#    - Annotates nodes with gene_name and in_pathway flag
#    - Saves as GraphML + GPickle for fast reload
#    - Runs a TP53 → AKT1 shortest-path smoke test

# Fixes vs previous version:
#   - Arrowheads now guaranteed on ALL edges (zero-length vector guard added)
#   - Isolated nodes removed after pruning
#   - Source/target edges always preserved
#   - Legend moved below axes (no overlap with RTK row)
#   - Zorder stack: bands(0) → edges(2) → boxes(3) → labels(4)
# """

# import os
# import math
# import pickle
# import networkx as nx
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from matplotlib.patches import FancyBboxPatch
# import numpy as np
# from collections import defaultdict

# # ── paths ─────────────────────────────────────────────────────────────────────
# INTERACTOME     = os.path.join("data", "PathLinker_2018_human-ppi-weighted-cap0_75.txt")
# OUT_FIG         = os.path.join("outputs", "figures", "graph_overview.png")
# OUT_PATHWAY_TXT = os.path.join("outputs", "txt", "pathway_proteins.txt")
# OUT_GRAPHML     = os.path.join("outputs", "graphs", "full_interactome.graphml")
# OUT_GPICKLE     = os.path.join("outputs", "graphs", "full_interactome.gpickle")
# os.makedirs(os.path.dirname(OUT_FIG),         exist_ok=True)
# os.makedirs(os.path.dirname(OUT_PATHWAY_TXT), exist_ok=True)
# os.makedirs(os.path.dirname(OUT_GRAPHML),     exist_ok=True)

# KEGG_PATHWAY_ID = "hsa04151"
# PATHWAY_NAME    = "PI3K-AKT Signaling Pathway"
# SOURCE_ID       = "P04637"   # TP53
# TARGET_ID       = "P31749"   # AKT1
# CONF_THRESHOLD  = 0.75

# # Edge pruning: keep top-N edges; isolated nodes are then removed automatically.
# # Set PRUNE_EDGES = False to draw all edges (denser but noisier).
# PRUNE_EDGES = True
# MAX_EDGES   = 200

# # ── UniProt → Gene name map ────────────────────────────────────────────────────
# GENE_NAMES = {
#     "P31749": "AKT1",    "P31751": "AKT2",    "Q9Y243": "AKT3",
#     "P42336": "PIK3CA",  "P42338": "PIK3CB",  "O00329": "PIK3CD",  "P48736": "PIK3CG",
#     "P27986": "PIK3R1",  "O00459": "PIK3R2",  "Q92569": "PIK3R3",  "Q8WYR1": "PIK3R5",
#     "P60484": "PTEN",    "Q6ZVD8": "PHLPP1",  "O60346": "PHLPP2",
#     "O15327": "INPP4B",  "Q92835": "INPP5D",
#     "P42345": "MTOR",    "Q8TB45": "RICTOR",  "Q6R327": "MLST8",
#     "Q8N122": "DEPTOR",  "Q9BPZ7": "RPTOR",
#     "Q92574": "TSC1",    "P49815": "TSC2",    "Q9Y3Q8": "RHEB",
#     "O15530": "PDPK1",
#     "P00533": "EGFR",    "P04626": "ERBB2",   "P21860": "ERBB3",   "Q15303": "ERBB4",
#     "P35968": "KDR",     "P17948": "FLT1",    "P36888": "FLT3",
#     "P09619": "PDGFRB",  "P09110": "PDGFRA",  "P10721": "KIT",
#     "P08069": "IGF1R",   "P06213": "INSR",
#     "P04085": "PDGFA",   "P01127": "PDGFB",   "P15692": "VEGFA",
#     "P01116": "KRAS",    "P01112": "HRAS",    "P01111": "NRAS",
#     "P15056": "BRAF",    "P04049": "RAF1",
#     "Q02750": "MAP2K1",  "P36507": "MAP2K2",
#     "P27361": "MAPK3",   "P28482": "MAPK1",
#     "P20936": "RASA1",   "Q99490": "NF1",
#     "P29353": "SHC1",    "P62993": "GRB2",    "Q07889": "SOS1",    "Q9Y6I9": "SOS2",
#     "P35568": "IRS1",    "P35570": "IRS2",
#     "P23443": "RPS6KB1", "P23444": "RPS6KB2", "Q13541": "EIF4EBP1",
#     "P49840": "GSK3A",   "P49841": "GSK3B",
#     "Q12778": "FOXO1",   "O43524": "FOXO3",   "Q9UKT9": "FOXO4",
#     "P04637": "TP53",    "Q00987": "MDM2",    "O15151": "MDM4",
#     "P24385": "CCND1",   "P11802": "CDK4",    "P24941": "CDK2",
#     "P38936": "CDKN1A",  "P46527": "CDKN1B",
#     "Q07812": "BAX",     "Q07817": "BCL2L1",  "P10415": "BCL2",
#     "P19838": "NFKB1",   "Q04206": "RELA",
#     "O14920": "IKBKB",   "Q15653": "NFKBIA",
#     "P17612": "PRKACA",  "P22694": "PRKACB",  "P16220": "CREB1",
#     "Q13131": "PRKAA1",  "P54646": "PRKAA2",
#     "P51812": "RPS6KA1", "Q15208": "RPS6KA2",
#     "Q16665": "HIF1A",   "P27540": "ARNT",
#     "P51955": "NEK6",    "Q8WYQ9": "THEM4",
# }

# # ── Functional group colors ────────────────────────────────────────────────────
# GROUP_COLORS = {
#     "RTK":        "#a8d8a8",
#     "RAS":        "#f9c784",
#     "PI3K":       "#a0c4ff",
#     "PTEN":       "#ff9aa2",
#     "AKT":        "#2196F3",
#     "mTOR":       "#b5ead7",
#     "TSC":        "#c9b1ff",
#     "FOXO":       "#ffd6e7",
#     "Cell_cycle": "#ffe0ac",
#     "Apoptosis":  "#ff6b6b",
#     "NFkB":       "#d4a5a5",
#     "TP53":       "#e74c3c",
#     "Adaptor":    "#e8e8e8",
#     "Other":      "#d0d0d0",
# }

# PROTEIN_GROUP = {
#     **{p: "RTK" for p in [
#         "P00533","P04626","P21860","Q15303","P35968","P17948",
#         "P36888","P09619","P09110","P10721","P08069","P06213",
#         "P04085","P01127","P15692"]},
#     **{p: "RAS" for p in [
#         "P01116","P01112","P01111","P15056","P04049",
#         "Q02750","P36507","P27361","P28482","P20936","Q99490"]},
#     **{p: "PI3K" for p in [
#         "P42336","P42338","O00329","P48736",
#         "P27986","O00459","Q92569","Q8WYR1",
#         "O15530","P35568","P35570"]},
#     **{p: "PTEN" for p in ["P60484","Q6ZVD8","O60346","O15327","Q92835"]},
#     **{p: "AKT"  for p in ["P31749","P31751","Q9Y243"]},
#     **{p: "mTOR" for p in ["P42345","Q8TB45","Q6R327","Q8N122","Q9BPZ7"]},
#     **{p: "TSC"  for p in ["Q92574","P49815","Q9Y3Q8"]},
#     **{p: "Adaptor" for p in ["P29353","P62993","Q07889","Q9Y6I9"]},
#     **{p: "FOXO" for p in [
#         "Q12778","O43524","Q9UKT9","P49840","P49841",
#         "P23443","P23444","Q13541",
#         "P17612","P22694","P16220",
#         "Q13131","P54646","P51812","Q15208",
#         "Q16665","P27540","P51955","Q8WYQ9"]},
#     **{p: "Cell_cycle" for p in ["P24385","P11802","P24941","P38936","P46527"]},
#     **{p: "Apoptosis"  for p in ["Q07812","Q07817","P10415"]},
#     **{p: "NFkB"       for p in ["P19838","Q04206","O14920","Q15653"]},
#     **{p: "TP53"       for p in ["P04637","Q00987","O15151"]},
# }

# # ── Layer assignment ───────────────────────────────────────────────────────────
# LAYER_MAP = {
#     **{p: 0 for p in [
#         "P00533","P04626","P21860","Q15303","P35968","P17948",
#         "P36888","P09619","P09110","P10721","P08069","P06213",
#         "P04085","P01127","P15692"]},
#     **{p: 1 for p in [
#         "P29353","P62993","Q07889","Q9Y6I9",
#         "P01116","P01112","P01111","P20936","Q99490",
#         "P35568","P35570"]},
#     **{p: 2 for p in [
#         "P42336","P42338","O00329","P48736",
#         "P27986","O00459","Q92569","Q8WYR1",
#         "P60484","O15327","Q92835","Q6ZVD8","O60346"]},
#     **{p: 3 for p in ["P31749","P31751","Q9Y243","O15530","Q13131","P54646"]},
#     **{p: 4 for p in [
#         "P42345","Q8TB45","Q6R327","Q8N122","Q9BPZ7",
#         "Q92574","P49815","Q9Y3Q8",
#         "P15056","P04049","Q02750","P36507",
#         "O14920","Q15653"]},
#     **{p: 5 for p in [
#         "P23443","P23444","Q13541",
#         "Q12778","O43524","Q9UKT9",
#         "P49840","P49841","P27361","P28482",
#         "P24385","P11802","P24941",
#         "P17612","P22694","P16220",
#         "P51812","Q15208","Q16665","P27540",
#         "P51955","Q8WYQ9"]},
#     **{p: 6 for p in [
#         "P04637","Q00987","O15151",
#         "P38936","P46527",
#         "Q07812","Q07817","P10415",
#         "P19838","Q04206"]},
# }


# # ── Helpers ────────────────────────────────────────────────────────────────────

# def get_pi3k_akt_uniprot_ids():
#     return set(GENE_NAMES.keys())


# def load_pathway_subgraph(filepath, pathway_proteins, conf_threshold=0.75):
#     G_pathway = nx.DiGraph()
#     G_full    = nx.DiGraph()
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 3:
#                 continue
#             src, dst, conf = parts[0], parts[1], float(parts[2])
#             G_full.add_node(src)
#             G_full.add_node(dst)
#             if (src in pathway_proteins and dst in pathway_proteins
#                     and conf >= conf_threshold):
#                 G_pathway.add_edge(src, dst, weight=conf, cost=1.0 - conf)
#     found = pathway_proteins & set(G_full.nodes())
#     print(f"[viz] Pathway proteins in interactome : {len(found)}")
#     print(f"[viz] After conf >= {conf_threshold} : "
#           f"{G_pathway.number_of_nodes()} nodes, "
#           f"{G_pathway.number_of_edges()} edges")
#     return G_full, G_pathway, found


# def save_pathway_proteins(found_proteins, source, target, out_path):
#     with open(out_path, "w") as fh:
#         fh.write(f"Pathway : {PATHWAY_NAME} ({KEGG_PATHWAY_ID})\n")
#         fh.write(f"Proteins found in interactome: {len(found_proteins)}\n\n")
#         for pid in sorted(found_proteins):
#             tag = (" <- TP53 (source)" if pid == source else
#                    " <- AKT1 (target)" if pid == target else "")
#             fh.write(f"{pid}\t{GENE_NAMES.get(pid, pid)}{tag}\n")
#     print(f"[viz] Protein list saved -> {out_path}")


# def prune_to_top_edges(G, max_edges, must_keep_nodes=None):
#     """
#     Keep the top-N edges by weight.

#     Edges touching must_keep_nodes (source / target) are always included first,
#     then the remaining budget is filled with the next highest-confidence edges.
#     After pruning, any node with degree == 0 that is NOT in must_keep_nodes
#     is removed so there are no floating isolated boxes.
#     """
#     must_keep_nodes = must_keep_nodes or set()

#     all_edges = sorted(G.edges(data=True),
#                        key=lambda e: e[2].get("weight", 0), reverse=True)

#     pinned = [(u, v, d) for u, v, d in all_edges
#               if u in must_keep_nodes or v in must_keep_nodes]
#     rest   = [(u, v, d) for u, v, d in all_edges
#               if u not in must_keep_nodes and v not in must_keep_nodes]

#     budget = max(max_edges, len(pinned))
#     chosen = pinned + rest[:budget - len(pinned)]

#     G_pruned = nx.DiGraph()
#     G_pruned.add_nodes_from(G.nodes())
#     G_pruned.add_edges_from((u, v, d) for u, v, d in chosen)

#     # Drop isolated nodes that are not the source or target
#     isolated = [n for n in list(G_pruned.nodes())
#                 if G_pruned.degree(n) == 0 and n not in must_keep_nodes]
#     G_pruned.remove_nodes_from(isolated)

#     print(f"[viz] Pruned: {G.number_of_edges()} -> {G_pruned.number_of_edges()} edges | "
#           f"{G.number_of_nodes()} -> {G_pruned.number_of_nodes()} nodes "
#           f"({len(isolated)} isolated nodes removed)")
#     return G_pruned


# def kegg_style_layout(G):
#     """Spread nodes evenly within their assigned layer band."""
#     layer_nodes = defaultdict(list)
#     for node in G.nodes():
#         layer_nodes[LAYER_MAP.get(node, 3)].append(node)

#     pos      = {}
#     n_layers = 7
#     layer_h  = 1.0 / (n_layers - 1)

#     for layer, nodes in layer_nodes.items():
#         nodes = sorted(nodes, key=lambda n: GENE_NAMES.get(n, n))
#         n = len(nodes)
#         for i, node in enumerate(nodes):
#             pos[node] = ((i + 0.5) / n, 1.0 - layer * layer_h)
#     return pos


# # ── Drawing ────────────────────────────────────────────────────────────────────

# def draw_kegg_style(G, source, target, conf_threshold, out_path):
#     if G.number_of_nodes() == 0:
#         print("[viz] Nothing to draw.")
#         return

#     pos = kegg_style_layout(G)

#     FIG_W, FIG_H = 26, 18
#     fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
#     ax.set_facecolor("#f8f9fa")
#     fig.patch.set_facecolor("#f8f9fa")

#     # ── layer band backgrounds ─────────────────────────────────────────────────
#     layer_labels = {
#         0: "Growth Factor Receptors (RTKs)",
#         1: "Signal Adaptors / RAS GTPases",
#         2: "PI3K Complex / PTEN",
#         3: "AKT Hub / PDK1 / AMPK",
#         4: "mTOR Complex / TSC / RAF-MEK",
#         5: "Downstream Effectors",
#         6: "Nuclear / Apoptotic Targets",
#     }
#     n_layers    = 7
#     layer_h     = 1.0 / (n_layers - 1)
#     band_colors = ["#eaf4fb","#eafaf1","#fef9e7","#e8f8f5",
#                    "#f4ecf7","#fdfefe","#fdedec"]

#     for i in range(n_layers):
#         y_center = 1.0 - i * layer_h
#         ax.add_patch(mpatches.FancyBboxPatch(
#             (-0.01, y_center - layer_h * 0.48), 1.02, layer_h * 0.96,
#             boxstyle="round,pad=0.005",
#             linewidth=0.5, edgecolor="#cccccc",
#             facecolor=band_colors[i], alpha=0.55,
#             transform=ax.transData, zorder=0))
#         ax.text(1.035, y_center, layer_labels[i],
#                 fontsize=7.5, va="center", color="#555",
#                 style="italic", fontfamily="DejaVu Sans",
#                 transform=ax.transData)

#     # ── node box dimensions ────────────────────────────────────────────────────
#     BOX_W = 0.046
#     BOX_H = 0.026

#     # ── edges ─────────────────────────────────────────────────────────────────
#     edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
#     if edge_weights:
#         wmin, wmax = min(edge_weights), max(edge_weights)
#         wrange = wmax - wmin if wmax > wmin else 1.0

#         for (u, v), w in zip(G.edges(), edge_weights):
#             xu, yu = pos[u]
#             xv, yv = pos[v]
#             nw = (w - wmin) / wrange

#             r  = 0.15 + 0.65 * nw
#             g  = 0.25 - 0.15 * nw
#             b  = 0.85 - 0.65 * nw
#             lw = 1.0  + 2.0  * nw

#             # Offset endpoints to box edges so arrowhead is never hidden
#             # behind the target node patch.
#             # Guard: if source and target occupy the exact same position
#             # (same-layer edge), fall back to a pure arc with no offset.
#             dx, dy = xv - xu, yv - yu
#             dist   = np.hypot(dx, dy)
#             if dist < 1e-6:
#                 # Same-position nodes: draw a looping arc, skip offset
#                 sx, sy, ex, ey = xu, yu, xv, yv
#                 conn_style = "arc3,rad=0.35"
#             else:
#                 sx = xu + (BOX_H / 2)          * dx / dist
#                 sy = yu + (BOX_H / 2)          * dy / dist
#                 ex = xv - (BOX_H / 2 + 0.008) * dx / dist
#                 ey = yv - (BOX_H / 2 + 0.008) * dy / dist
#                 conn_style = "arc3,rad=0.06"

#             ax.annotate(
#                 "", xy=(ex, ey), xytext=(sx, sy),
#                 arrowprops=dict(
#                     arrowstyle="-|>",
#                     color=(r, g, b, 0.80),
#                     lw=lw,
#                     mutation_scale=14,
#                     connectionstyle=conn_style),
#                 zorder=2)

#     # ── nodes ─────────────────────────────────────────────────────────────────
#     for node in G.nodes():
#         x, y  = pos[node]
#         gene  = GENE_NAMES.get(node, node)
#         group = PROTEIN_GROUP.get(node, "Other")
#         color = GROUP_COLORS.get(group, "#d0d0d0")

#         if node == source:
#             fc, ec, tc, lw = "#e74c3c", "#922b21", "white", 2.0
#         elif node == target:
#             fc, ec, tc, lw = "#27ae60", "#1a6b3c", "white", 2.0
#         else:
#             fc, ec, tc, lw = color, "#555555", "#111111", 0.8

#         ax.add_patch(FancyBboxPatch(
#             (x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
#             boxstyle="round,pad=0.003",
#             facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3))

#         ax.text(x, y, gene,
#                 ha="center", va="center",
#                 fontsize=5.8, fontweight="bold",
#                 color=tc, fontfamily="DejaVu Sans Mono",
#                 zorder=4, clip_on=True)

#     # ── confidence colorbar ────────────────────────────────────────────────────
#     from matplotlib.colors import LinearSegmentedColormap
#     cmap = LinearSegmentedColormap.from_list(
#         "conf", [(0.15, 0.25, 0.85), (0.80, 0.10, 0.20)], N=256)
#     sm = plt.cm.ScalarMappable(
#         cmap=cmap, norm=plt.Normalize(vmin=conf_threshold, vmax=1.0))
#     sm.set_array([])
#     cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.35])
#     cb = fig.colorbar(sm, cax=cbar_ax)
#     cb.set_label("Interaction\nConfidence", fontsize=8, labelpad=6)
#     cb.ax.tick_params(labelsize=7)

#     # ── legend — placed BELOW the graph so it never overlaps layer-0 nodes ────
#     shown_groups = sorted({PROTEIN_GROUP.get(n, "Other") for n in G.nodes()})
#     legend_patches = [
#         mpatches.Patch(facecolor=GROUP_COLORS.get(g, "#d0d0d0"),
#                        edgecolor="#555", linewidth=0.6,
#                        label=g.replace("_", " "))
#         for g in shown_groups
#     ] + [
#         mpatches.Patch(facecolor="#e74c3c", edgecolor="#922b21",
#                        linewidth=1.2, label="TP53 – Source"),
#         mpatches.Patch(facecolor="#27ae60", edgecolor="#1a6b3c",
#                        linewidth=1.2, label="AKT1 – Target"),
#     ]
#     ax.legend(handles=legend_patches,
#               loc="lower center",
#               bbox_to_anchor=(0.45, -0.10),   # centred below axes
#               fontsize=6.5, framealpha=0.9, ncol=5,
#               title="Functional Group", title_fontsize=7,
#               edgecolor="#aaa")

#     # ── title & axes ──────────────────────────────────────────────────────────
#     edge_note = (f" (top {G.number_of_edges()} by confidence)"
#                  if PRUNE_EDGES else "")
#     ax.set_title(
#         f"PPI Sub-Network: {PATHWAY_NAME}  ({KEGG_PATHWAY_ID})\n"
#         f"{G.number_of_nodes()} proteins  ·  "
#         f"{G.number_of_edges()} interactions{edge_note}"
#         f"  ·  Confidence ≥ {conf_threshold}  ·  Hierarchical layout",
#         fontsize=13, fontweight="bold", pad=16, fontfamily="DejaVu Sans")

#     ax.set_xlim(-0.03, 1.03)
#     ax.set_ylim(-0.07, 1.07)
#     ax.axis("off")

#     plt.tight_layout(rect=[0, 0.07, 0.91, 1])   # leave room below for legend
#     plt.savefig(out_path, dpi=180, bbox_inches="tight",
#                 facecolor=fig.get_facecolor())
#     plt.close()
#     print(f"[viz] Figure saved -> {out_path}")


# # ═══════════════════════════════════════════════════════════════════════════════
# # FULL INTERACTOME GRAPH BUILDER
# # Loads every edge in the PathLinker file (no protein filter, no conf cutoff).
# # Saves GraphML + GPickle for fast reload in shortest-path scripts.
# # ═══════════════════════════════════════════════════════════════════════════════

# PI3K_AKT_PROTEINS = set(GENE_NAMES.keys())   # reuse the map defined above


# def build_full_graph(filepath):
#     """
#     Parse the entire PathLinker interactome into a directed weighted DiGraph.

#     Edge attributes
#     ---------------
#     weight          : raw confidence score  (0 – 1, higher = more confident)
#     cost            : 1 - weight            (use as Dijkstra edge weight —
#                                              lower cost = more confident path)
#     neg_log_weight  : -log(weight)          (use for max-product / probabilistic
#                                              path calculations)

#     Node attributes
#     ---------------
#     gene_name   : human-readable symbol if in GENE_NAMES, else the UniProt ID
#     in_pathway  : True if the node belongs to the curated PI3K-AKT set
#     """
#     G = nx.DiGraph()
#     skipped = 0

#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 3:
#                 skipped += 1
#                 continue
#             src, dst = parts[0], parts[1]
#             try:
#                 conf = float(parts[2])
#             except ValueError:
#                 skipped += 1
#                 continue

#             G.add_edge(src, dst,
#                        weight=conf,
#                        cost=round(1.0 - conf, 6),
#                        neg_log_weight=round(-math.log(max(conf, 1e-9)), 6))

#     # Annotate every node
#     for node in G.nodes():
#         G.nodes[node]["gene_name"]  = GENE_NAMES.get(node, node)
#         G.nodes[node]["in_pathway"] = node in PI3K_AKT_PROTEINS

#     if skipped:
#         print(f"[full] Skipped {skipped} malformed lines.")

#     return G


# def save_full_graph(G, graphml_path, gpickle_path):
#     """Save graph in two formats: GraphML (portable) and GPickle (fast)."""
#     nx.write_graphml(G, graphml_path)
#     print(f"[full] GraphML saved  -> {graphml_path}")

#     with open(gpickle_path, "wb") as f:
#         pickle.dump(G, f, protocol=pickle.HIGHEST_PROTOCOL)
#     print(f"[full] GPickle saved  -> {gpickle_path}")


# def print_full_graph_summary(G):
#     n_nodes   = G.number_of_nodes()
#     n_edges   = G.number_of_edges()
#     n_pathway = sum(1 for n in G.nodes() if G.nodes[n].get("in_pathway"))
#     weights   = [d["weight"] for _, _, d in G.edges(data=True)]
#     print("\n── Full interactome graph ──────────────────────────────────")
#     print(f"  Nodes              : {n_nodes:,}")
#     print(f"  Edges              : {n_edges:,}")
#     print(f"  PI3K-AKT nodes     : {n_pathway}")
#     print(f"  Confidence range   : {min(weights):.4f} – {max(weights):.4f}"
#           f"  (mean {sum(weights)/len(weights):.4f})")
#     print(f"  Is directed        : {G.is_directed()}")
#     n_wcc = nx.number_weakly_connected_components(G)
#     print(f"  Weakly conn. comps : {n_wcc}")
#     print("────────────────────────────────────────────────────────────\n")


# def shortest_path_smoke_test(G, source=SOURCE_ID, target=TARGET_ID):
#     """
#     Dijkstra smoke test: TP53 → AKT1 via 'cost' weight.
#     Prints the path and total cost so you can verify the graph is usable.
#     """
#     src_name = GENE_NAMES.get(source, source)
#     tgt_name = GENE_NAMES.get(target, target)
#     if not G.has_node(source):
#         print(f"[full] Source {source} ({src_name}) not in graph.")
#         return
#     if not G.has_node(target):
#         print(f"[full] Target {target} ({tgt_name}) not in graph.")
#         return
#     try:
#         length, path = nx.single_source_dijkstra(G, source, target, weight="cost")
#         gene_path = [GENE_NAMES.get(p, p) for p in path]
#         print(f"[full] Shortest path  {src_name} → {tgt_name}")
#         print(f"       Hops : {len(path)-1}   Total cost : {length:.4f}")
#         print(f"       Path : {' → '.join(gene_path)}")
#     except nx.NetworkXNoPath:
#         print(f"[full] No directed path from {src_name} to {tgt_name}.")


# # ── entry point ───────────────────────────────────────────────────────────────
# if __name__ == "__main__":

#     # ── 1. PI3K-AKT sub-network visualization ─────────────────────────────────
#     print(f"[viz] Building curated UniProt set for {KEGG_PATHWAY_ID}...")
#     pathway_proteins         = get_pi3k_akt_uniprot_ids()
#     G_full, G_pathway, found = load_pathway_subgraph(
#                                    INTERACTOME, pathway_proteins,
#                                    conf_threshold=CONF_THRESHOLD)
#     save_pathway_proteins(found, SOURCE_ID, TARGET_ID, OUT_PATHWAY_TXT)

#     if G_pathway.number_of_nodes() == 0:
#         print(f"[viz] WARNING: 0 nodes. "
#               f"Try lowering CONF_THRESHOLD (now {CONF_THRESHOLD}).")
#     else:
#         if PRUNE_EDGES:
#             G_draw = prune_to_top_edges(
#                 G_pathway, MAX_EDGES,
#                 must_keep_nodes={SOURCE_ID, TARGET_ID})
#         else:
#             G_draw = G_pathway
#         draw_kegg_style(G_draw, SOURCE_ID, TARGET_ID, CONF_THRESHOLD, OUT_FIG)

#     # ── 2. Full interactome DiGraph (no filter) ────────────────────────────────
#     print("\n[full] Building complete interactome graph (all proteins, all edges)...")
#     G_interactome = build_full_graph(INTERACTOME)
#     print_full_graph_summary(G_interactome)
#     save_full_graph(G_interactome, OUT_GRAPHML, OUT_GPICKLE)
#     shortest_path_smoke_test(G_interactome)

#     print("\n[full] To reload in another script:")
#     print("       import pickle")
#     print(f"       with open('{OUT_GPICKLE}', 'rb') as f: G = pickle.load(f)")
#     print("       path = nx.shortest_path(G, 'P04637', 'P31749', weight='cost')")
#     print("\n[load_graph] All done.")

"""
load_graph.py
=============
Two things in one file:

1. KEGG-style PI3K-AKT sub-network visualization
   - Colored rectangular boxes grouped by functional category
   - Directed edges with visible arrowheads (offset to box edges)
   - Hierarchical layer layout, confidence colorbar, legend below graph
   - Top-N edge pruning with isolated-node removal

2. Full PathLinker interactome DiGraph builder  (build_full_graph)
   - Loads every edge in the file, no protein filter, no confidence cutoff
   - Stores weight / cost / neg_log_weight on every edge
   - Annotates nodes with gene_name and in_pathway flag
   - Saves as GraphML + GPickle for fast reload
   - Runs a TP53 → AKT1 shortest-path smoke test

Public API (importable from other scripts):
   from load_graph import load_graph
   G = load_graph("data/PathLinker_2018_human-ppi-weighted-cap0_75.txt")

   load_graph() is a clean wrapper around build_full_graph() that:
       - Returns a nx.DiGraph with edge attributes: weight, cost, neg_log_weight
       - Prints a summary of the loaded graph
       - Accepts any file path as argument

Fixes vs previous version:
  - Arrowheads now guaranteed on ALL edges (zero-length vector guard added)
  - Isolated nodes removed after pruning
  - Source/target edges always preserved
  - Legend moved below axes (no overlap with RTK row)
  - Zorder stack: bands(0) → edges(2) → boxes(3) → labels(4)
"""

import os
import math
import pickle
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from collections import defaultdict

# ── paths ─────────────────────────────────────────────────────────────────────
INTERACTOME     = os.path.join("data", "PathLinker_2018_human-ppi-weighted-cap0_75.txt")
OUT_FIG         = os.path.join("outputs", "figures", "graph_overview.png")
OUT_PATHWAY_TXT = os.path.join("outputs", "txt", "pathway_proteins.txt")
OUT_GRAPHML     = os.path.join("outputs", "graphs", "full_interactome.graphml")
OUT_GPICKLE     = os.path.join("outputs", "graphs", "full_interactome.gpickle")
os.makedirs(os.path.dirname(OUT_FIG),         exist_ok=True)
os.makedirs(os.path.dirname(OUT_PATHWAY_TXT), exist_ok=True)
os.makedirs(os.path.dirname(OUT_GRAPHML),     exist_ok=True)

KEGG_PATHWAY_ID = "hsa04151"
PATHWAY_NAME    = "PI3K-AKT Signaling Pathway"
SOURCE_ID       = "P04637"   # TP53
TARGET_ID       = "P31749"   # AKT1
CONF_THRESHOLD  = 0.75

# Edge pruning: keep top-N edges; isolated nodes are then removed automatically.
# Set PRUNE_EDGES = False to draw all edges (denser but noisier).
PRUNE_EDGES = True
MAX_EDGES   = 200

# ── UniProt → Gene name map ────────────────────────────────────────────────────
GENE_NAMES = {
    "P31749": "AKT1",    "P31751": "AKT2",    "Q9Y243": "AKT3",
    "P42336": "PIK3CA",  "P42338": "PIK3CB",  "O00329": "PIK3CD",  "P48736": "PIK3CG",
    "P27986": "PIK3R1",  "O00459": "PIK3R2",  "Q92569": "PIK3R3",  "Q8WYR1": "PIK3R5",
    "P60484": "PTEN",    "Q6ZVD8": "PHLPP1",  "O60346": "PHLPP2",
    "O15327": "INPP4B",  "Q92835": "INPP5D",
    "P42345": "MTOR",    "Q8TB45": "RICTOR",  "Q6R327": "MLST8",
    "Q8N122": "DEPTOR",  "Q9BPZ7": "RPTOR",
    "Q92574": "TSC1",    "P49815": "TSC2",    "Q9Y3Q8": "RHEB",
    "O15530": "PDPK1",
    "P00533": "EGFR",    "P04626": "ERBB2",   "P21860": "ERBB3",   "Q15303": "ERBB4",
    "P35968": "KDR",     "P17948": "FLT1",    "P36888": "FLT3",
    "P09619": "PDGFRB",  "P09110": "PDGFRA",  "P10721": "KIT",
    "P08069": "IGF1R",   "P06213": "INSR",
    "P04085": "PDGFA",   "P01127": "PDGFB",   "P15692": "VEGFA",
    "P01116": "KRAS",    "P01112": "HRAS",    "P01111": "NRAS",
    "P15056": "BRAF",    "P04049": "RAF1",
    "Q02750": "MAP2K1",  "P36507": "MAP2K2",
    "P27361": "MAPK3",   "P28482": "MAPK1",
    "P20936": "RASA1",   "Q99490": "NF1",
    "P29353": "SHC1",    "P62993": "GRB2",    "Q07889": "SOS1",    "Q9Y6I9": "SOS2",
    "P35568": "IRS1",    "P35570": "IRS2",
    "P23443": "RPS6KB1", "P23444": "RPS6KB2", "Q13541": "EIF4EBP1",
    "P49840": "GSK3A",   "P49841": "GSK3B",
    "Q12778": "FOXO1",   "O43524": "FOXO3",   "Q9UKT9": "FOXO4",
    "P04637": "TP53",    "Q00987": "MDM2",    "O15151": "MDM4",
    "P24385": "CCND1",   "P11802": "CDK4",    "P24941": "CDK2",
    "P38936": "CDKN1A",  "P46527": "CDKN1B",
    "Q07812": "BAX",     "Q07817": "BCL2L1",  "P10415": "BCL2",
    "P19838": "NFKB1",   "Q04206": "RELA",
    "O14920": "IKBKB",   "Q15653": "NFKBIA",
    "P17612": "PRKACA",  "P22694": "PRKACB",  "P16220": "CREB1",
    "Q13131": "PRKAA1",  "P54646": "PRKAA2",
    "P51812": "RPS6KA1", "Q15208": "RPS6KA2",
    "Q16665": "HIF1A",   "P27540": "ARNT",
    "P51955": "NEK6",    "Q8WYQ9": "THEM4",
}

# ── Functional group colors ────────────────────────────────────────────────────
GROUP_COLORS = {
    "RTK":        "#a8d8a8",
    "RAS":        "#f9c784",
    "PI3K":       "#a0c4ff",
    "PTEN":       "#ff9aa2",
    "AKT":        "#2196F3",
    "mTOR":       "#b5ead7",
    "TSC":        "#c9b1ff",
    "FOXO":       "#ffd6e7",
    "Cell_cycle": "#ffe0ac",
    "Apoptosis":  "#ff6b6b",
    "NFkB":       "#d4a5a5",
    "TP53":       "#e74c3c",
    "Adaptor":    "#e8e8e8",
    "Other":      "#d0d0d0",
}

PROTEIN_GROUP = {
    **{p: "RTK" for p in [
        "P00533","P04626","P21860","Q15303","P35968","P17948",
        "P36888","P09619","P09110","P10721","P08069","P06213",
        "P04085","P01127","P15692"]},
    **{p: "RAS" for p in [
        "P01116","P01112","P01111","P15056","P04049",
        "Q02750","P36507","P27361","P28482","P20936","Q99490"]},
    **{p: "PI3K" for p in [
        "P42336","P42338","O00329","P48736",
        "P27986","O00459","Q92569","Q8WYR1",
        "O15530","P35568","P35570"]},
    **{p: "PTEN" for p in ["P60484","Q6ZVD8","O60346","O15327","Q92835"]},
    **{p: "AKT"  for p in ["P31749","P31751","Q9Y243"]},
    **{p: "mTOR" for p in ["P42345","Q8TB45","Q6R327","Q8N122","Q9BPZ7"]},
    **{p: "TSC"  for p in ["Q92574","P49815","Q9Y3Q8"]},
    **{p: "Adaptor" for p in ["P29353","P62993","Q07889","Q9Y6I9"]},
    **{p: "FOXO" for p in [
        "Q12778","O43524","Q9UKT9","P49840","P49841",
        "P23443","P23444","Q13541",
        "P17612","P22694","P16220",
        "Q13131","P54646","P51812","Q15208",
        "Q16665","P27540","P51955","Q8WYQ9"]},
    **{p: "Cell_cycle" for p in ["P24385","P11802","P24941","P38936","P46527"]},
    **{p: "Apoptosis"  for p in ["Q07812","Q07817","P10415"]},
    **{p: "NFkB"       for p in ["P19838","Q04206","O14920","Q15653"]},
    **{p: "TP53"       for p in ["P04637","Q00987","O15151"]},
}

# ── Layer assignment ───────────────────────────────────────────────────────────
LAYER_MAP = {
    **{p: 0 for p in [
        "P00533","P04626","P21860","Q15303","P35968","P17948",
        "P36888","P09619","P09110","P10721","P08069","P06213",
        "P04085","P01127","P15692"]},
    **{p: 1 for p in [
        "P29353","P62993","Q07889","Q9Y6I9",
        "P01116","P01112","P01111","P20936","Q99490",
        "P35568","P35570"]},
    **{p: 2 for p in [
        "P42336","P42338","O00329","P48736",
        "P27986","O00459","Q92569","Q8WYR1",
        "P60484","O15327","Q92835","Q6ZVD8","O60346"]},
    **{p: 3 for p in ["P31749","P31751","Q9Y243","O15530","Q13131","P54646"]},
    **{p: 4 for p in [
        "P42345","Q8TB45","Q6R327","Q8N122","Q9BPZ7",
        "Q92574","P49815","Q9Y3Q8",
        "P15056","P04049","Q02750","P36507",
        "O14920","Q15653"]},
    **{p: 5 for p in [
        "P23443","P23444","Q13541",
        "Q12778","O43524","Q9UKT9",
        "P49840","P49841","P27361","P28482",
        "P24385","P11802","P24941",
        "P17612","P22694","P16220",
        "P51812","Q15208","Q16665","P27540",
        "P51955","Q8WYQ9"]},
    **{p: 6 for p in [
        "P04637","Q00987","O15151",
        "P38936","P46527",
        "Q07812","Q07817","P10415",
        "P19838","Q04206"]},
}


# ══════════════════════════════════════════════════════════════════════════════
# PUBLIC API — this is what all other scripts import
# ══════════════════════════════════════════════════════════════════════════════

def load_graph(filepath: str) -> nx.DiGraph:
    """
    Load the PathLinker interactome into a directed weighted graph.

    This is the main entry point used by all other scripts:
        from load_graph import load_graph
        G = load_graph("data/PathLinker_2018_human-ppi-weighted-cap0_75.txt")

    Internally calls build_full_graph() then prints a summary.

    Parameters
    ----------
    filepath : str
        Path to the tab-separated interactome file.

    Returns
    -------
    G : nx.DiGraph
        Directed graph where every edge carries:
            weight         — raw confidence score [0, 1]
            cost           — 1 - weight  (for Dijkstra cost minimization)
            neg_log_weight — -log(weight) (for probabilistic path scoring)
        Every node carries:
            gene_name  — human-readable gene symbol (or UniProt ID if unknown)
            in_pathway — True if the protein belongs to the PI3K-AKT set
    """
    print(f"\n[load_graph] Reading: {filepath}")
    G = build_full_graph(filepath)
    print_full_graph_summary(G)
    return G


# ── Helpers ────────────────────────────────────────────────────────────────────

def get_pi3k_akt_uniprot_ids():
    return set(GENE_NAMES.keys())


def load_pathway_subgraph(filepath, pathway_proteins, conf_threshold=0.75):
    G_pathway = nx.DiGraph()
    G_full    = nx.DiGraph()
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            src, dst, conf = parts[0], parts[1], float(parts[2])
            G_full.add_node(src)
            G_full.add_node(dst)
            if (src in pathway_proteins and dst in pathway_proteins
                    and conf >= conf_threshold):
                G_pathway.add_edge(src, dst, weight=conf, cost=1.0 - conf)
    found = pathway_proteins & set(G_full.nodes())
    print(f"[viz] Pathway proteins in interactome : {len(found)}")
    print(f"[viz] After conf >= {conf_threshold} : "
          f"{G_pathway.number_of_nodes()} nodes, "
          f"{G_pathway.number_of_edges()} edges")
    return G_full, G_pathway, found


def save_pathway_proteins(found_proteins, source, target, out_path):
    with open(out_path, "w") as fh:
        fh.write(f"Pathway : {PATHWAY_NAME} ({KEGG_PATHWAY_ID})\n")
        fh.write(f"Proteins found in interactome: {len(found_proteins)}\n\n")
        for pid in sorted(found_proteins):
            tag = (" <- TP53 (source)" if pid == source else
                   " <- AKT1 (target)" if pid == target else "")
            fh.write(f"{pid}\t{GENE_NAMES.get(pid, pid)}{tag}\n")
    print(f"[viz] Protein list saved -> {out_path}")


def prune_to_top_edges(G, max_edges, must_keep_nodes=None):
    """
    Keep the top-N edges by weight.

    Edges touching must_keep_nodes (source / target) are always included first,
    then the remaining budget is filled with the next highest-confidence edges.
    After pruning, any node with degree == 0 that is NOT in must_keep_nodes
    is removed so there are no floating isolated boxes.
    """
    must_keep_nodes = must_keep_nodes or set()

    all_edges = sorted(G.edges(data=True),
                       key=lambda e: e[2].get("weight", 0), reverse=True)

    pinned = [(u, v, d) for u, v, d in all_edges
              if u in must_keep_nodes or v in must_keep_nodes]
    rest   = [(u, v, d) for u, v, d in all_edges
              if u not in must_keep_nodes and v not in must_keep_nodes]

    budget = max(max_edges, len(pinned))
    chosen = pinned + rest[:budget - len(pinned)]

    G_pruned = nx.DiGraph()
    G_pruned.add_nodes_from(G.nodes())
    G_pruned.add_edges_from((u, v, d) for u, v, d in chosen)

    # Drop isolated nodes that are not the source or target
    isolated = [n for n in list(G_pruned.nodes())
                if G_pruned.degree(n) == 0 and n not in must_keep_nodes]
    G_pruned.remove_nodes_from(isolated)

    print(f"[viz] Pruned: {G.number_of_edges()} -> {G_pruned.number_of_edges()} edges | "
          f"{G.number_of_nodes()} -> {G_pruned.number_of_nodes()} nodes "
          f"({len(isolated)} isolated nodes removed)")
    return G_pruned


def kegg_style_layout(G):
    """Spread nodes evenly within their assigned layer band."""
    layer_nodes = defaultdict(list)
    for node in G.nodes():
        layer_nodes[LAYER_MAP.get(node, 3)].append(node)

    pos      = {}
    n_layers = 7
    layer_h  = 1.0 / (n_layers - 1)

    for layer, nodes in layer_nodes.items():
        nodes = sorted(nodes, key=lambda n: GENE_NAMES.get(n, n))
        n = len(nodes)
        for i, node in enumerate(nodes):
            pos[node] = ((i + 0.5) / n, 1.0 - layer * layer_h)
    return pos


# ── Drawing ────────────────────────────────────────────────────────────────────

def draw_kegg_style(G, source, target, conf_threshold, out_path):
    if G.number_of_nodes() == 0:
        print("[viz] Nothing to draw.")
        return

    pos = kegg_style_layout(G)

    FIG_W, FIG_H = 26, 18
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_facecolor("#f8f9fa")
    fig.patch.set_facecolor("#f8f9fa")

    # ── layer band backgrounds ─────────────────────────────────────────────────
    layer_labels = {
        0: "Growth Factor Receptors (RTKs)",
        1: "Signal Adaptors / RAS GTPases",
        2: "PI3K Complex / PTEN",
        3: "AKT Hub / PDK1 / AMPK",
        4: "mTOR Complex / TSC / RAF-MEK",
        5: "Downstream Effectors",
        6: "Nuclear / Apoptotic Targets",
    }
    n_layers    = 7
    layer_h     = 1.0 / (n_layers - 1)
    band_colors = ["#eaf4fb","#eafaf1","#fef9e7","#e8f8f5",
                   "#f4ecf7","#fdfefe","#fdedec"]

    for i in range(n_layers):
        y_center = 1.0 - i * layer_h
        ax.add_patch(mpatches.FancyBboxPatch(
            (-0.01, y_center - layer_h * 0.48), 1.02, layer_h * 0.96,
            boxstyle="round,pad=0.005",
            linewidth=0.5, edgecolor="#cccccc",
            facecolor=band_colors[i], alpha=0.55,
            transform=ax.transData, zorder=0))
        ax.text(1.035, y_center, layer_labels[i],
                fontsize=7.5, va="center", color="#555",
                style="italic", fontfamily="DejaVu Sans",
                transform=ax.transData)

    # ── node box dimensions ────────────────────────────────────────────────────
    BOX_W = 0.046
    BOX_H = 0.026

    # ── edges ─────────────────────────────────────────────────────────────────
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    if edge_weights:
        wmin, wmax = min(edge_weights), max(edge_weights)
        wrange = wmax - wmin if wmax > wmin else 1.0

        for (u, v), w in zip(G.edges(), edge_weights):
            xu, yu = pos[u]
            xv, yv = pos[v]
            nw = (w - wmin) / wrange

            r  = 0.15 + 0.65 * nw
            g  = 0.25 - 0.15 * nw
            b  = 0.85 - 0.65 * nw
            lw = 1.0  + 2.0  * nw

            dx, dy = xv - xu, yv - yu
            dist   = np.hypot(dx, dy)
            if dist < 1e-6:
                sx, sy, ex, ey = xu, yu, xv, yv
                conn_style = "arc3,rad=0.35"
            else:
                sx = xu + (BOX_H / 2)          * dx / dist
                sy = yu + (BOX_H / 2)          * dy / dist
                ex = xv - (BOX_H / 2 + 0.008) * dx / dist
                ey = yv - (BOX_H / 2 + 0.008) * dy / dist
                conn_style = "arc3,rad=0.06"

            ax.annotate(
                "", xy=(ex, ey), xytext=(sx, sy),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=(r, g, b, 0.80),
                    lw=lw,
                    mutation_scale=14,
                    connectionstyle=conn_style),
                zorder=2)

    # ── nodes ─────────────────────────────────────────────────────────────────
    for node in G.nodes():
        x, y  = pos[node]
        gene  = GENE_NAMES.get(node, node)
        group = PROTEIN_GROUP.get(node, "Other")
        color = GROUP_COLORS.get(group, "#d0d0d0")

        if node == source:
            fc, ec, tc, lw = "#e74c3c", "#922b21", "white", 2.0
        elif node == target:
            fc, ec, tc, lw = "#27ae60", "#1a6b3c", "white", 2.0
        else:
            fc, ec, tc, lw = color, "#555555", "#111111", 0.8

        ax.add_patch(FancyBboxPatch(
            (x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
            boxstyle="round,pad=0.003",
            facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3))

        ax.text(x, y, gene,
                ha="center", va="center",
                fontsize=5.8, fontweight="bold",
                color=tc, fontfamily="DejaVu Sans Mono",
                zorder=4, clip_on=True)

    # ── confidence colorbar ────────────────────────────────────────────────────
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "conf", [(0.15, 0.25, 0.85), (0.80, 0.10, 0.20)], N=256)
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=conf_threshold, vmax=1.0))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.92, 0.15, 0.012, 0.35])
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.set_label("Interaction\nConfidence", fontsize=8, labelpad=6)
    cb.ax.tick_params(labelsize=7)

    # ── legend ─────────────────────────────────────────────────────────────────
    shown_groups = sorted({PROTEIN_GROUP.get(n, "Other") for n in G.nodes()})
    legend_patches = [
        mpatches.Patch(facecolor=GROUP_COLORS.get(g, "#d0d0d0"),
                       edgecolor="#555", linewidth=0.6,
                       label=g.replace("_", " "))
        for g in shown_groups
    ] + [
        mpatches.Patch(facecolor="#e74c3c", edgecolor="#922b21",
                       linewidth=1.2, label="TP53 – Source"),
        mpatches.Patch(facecolor="#27ae60", edgecolor="#1a6b3c",
                       linewidth=1.2, label="AKT1 – Target"),
    ]
    ax.legend(handles=legend_patches,
              loc="lower center",
              bbox_to_anchor=(0.45, -0.10),
              fontsize=6.5, framealpha=0.9, ncol=5,
              title="Functional Group", title_fontsize=7,
              edgecolor="#aaa")

    # ── title & axes ──────────────────────────────────────────────────────────
    edge_note = (f" (top {G.number_of_edges()} by confidence)"
                 if PRUNE_EDGES else "")
    ax.set_title(
        f"PPI Sub-Network: {PATHWAY_NAME}  ({KEGG_PATHWAY_ID})\n"
        f"{G.number_of_nodes()} proteins  ·  "
        f"{G.number_of_edges()} interactions{edge_note}"
        f"  ·  Confidence ≥ {conf_threshold}  ·  Hierarchical layout",
        fontsize=13, fontweight="bold", pad=16, fontfamily="DejaVu Sans")

    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.07, 1.07)
    ax.axis("off")

    plt.tight_layout(rect=[0, 0.07, 0.91, 1])
    plt.savefig(out_path, dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[viz] Figure saved -> {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# FULL INTERACTOME GRAPH BUILDER
# ══════════════════════════════════════════════════════════════════════════════

PI3K_AKT_PROTEINS = set(GENE_NAMES.keys())


def build_full_graph(filepath: str) -> nx.DiGraph:
    """
    Parse the entire PathLinker interactome into a directed weighted DiGraph.

    Edge attributes
    ---------------
    weight         : raw confidence score [0, 1]  (higher = more confident)
    cost           : 1 - weight                   (lower = more confident,
                                                   use as Dijkstra weight)
    neg_log_weight : -log(weight)                 (for probabilistic scoring)

    Node attributes
    ---------------
    gene_name  : human-readable symbol if in GENE_NAMES, else UniProt ID
    in_pathway : True if the node belongs to the curated PI3K-AKT set
    """
    G       = nx.DiGraph()
    skipped = 0

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                skipped += 1
                continue
            src, dst = parts[0], parts[1]
            try:
                conf = float(parts[2])
            except ValueError:
                skipped += 1
                continue

            G.add_edge(src, dst,
                       weight=conf,
                       cost=round(1.0 - conf, 6),
                       neg_log_weight=round(-math.log(max(conf, 1e-9)), 6),
                       method=parts[3] if len(parts) > 3 else "unknown")  # ← add this

    # Annotate every node with gene name and pathway membership
    for node in G.nodes():
        G.nodes[node]["gene_name"]  = GENE_NAMES.get(node, node)
        G.nodes[node]["in_pathway"] = node in PI3K_AKT_PROTEINS

    if skipped:
        print(f"[load_graph] Skipped {skipped} malformed lines.")

    return G


def save_full_graph(G: nx.DiGraph,
                    graphml_path: str, gpickle_path: str):
    """Save graph in two formats: GraphML (portable) and GPickle (fast)."""
    nx.write_graphml(G, graphml_path)
    print(f"[load_graph] GraphML saved  -> {graphml_path}")

    with open(gpickle_path, "wb") as f:
        pickle.dump(G, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[load_graph] GPickle saved  -> {gpickle_path}")


def print_full_graph_summary(G: nx.DiGraph):
    """Print a concise summary of the loaded graph to the console."""
    n_nodes   = G.number_of_nodes()
    n_edges   = G.number_of_edges()
    n_pathway = sum(1 for n in G.nodes() if G.nodes[n].get("in_pathway"))
    weights   = [d["weight"] for _, _, d in G.edges(data=True)]
    print("── Graph summary ───────────────────────────────────────────")
    print(f"  Nodes              : {n_nodes:,}")
    print(f"  Edges              : {n_edges:,}")
    print(f"  PI3K-AKT nodes     : {n_pathway}")
    print(f"  Confidence range   : {min(weights):.4f} – {max(weights):.4f}"
          f"  (mean {sum(weights)/len(weights):.4f})")
    print(f"  Is directed        : {G.is_directed()}")
    n_wcc = nx.number_weakly_connected_components(G)
    print(f"  Weakly conn. comps : {n_wcc}")
    print("────────────────────────────────────────────────────────────\n")


def shortest_path_smoke_test(G: nx.DiGraph,
                              source=SOURCE_ID, target=TARGET_ID):
    """
    Quick Dijkstra smoke test: TP53 → AKT1 via 'cost' weight.
    Prints the path and total cost to verify the graph is usable.
    """
    src_name = GENE_NAMES.get(source, source)
    tgt_name = GENE_NAMES.get(target, target)
    if not G.has_node(source):
        print(f"[load_graph] Source {source} ({src_name}) not in graph.")
        return
    if not G.has_node(target):
        print(f"[load_graph] Target {target} ({tgt_name}) not in graph.")
        return
    try:
        length, path = nx.single_source_dijkstra(
            G, source, target, weight="cost")
        gene_path = [GENE_NAMES.get(p, p) for p in path]
        print(f"[load_graph] Smoke test  {src_name} → {tgt_name}")
        print(f"             Hops : {len(path)-1}   Total cost : {length:.4f}")
        print(f"             Path : {' → '.join(gene_path)}")
    except nx.NetworkXNoPath:
        print(f"[load_graph] No directed path from {src_name} to {tgt_name}.")


# ── entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # ── 1. PI3K-AKT sub-network visualization ─────────────────────────────────
    print(f"[viz] Building curated UniProt set for {KEGG_PATHWAY_ID}...")
    pathway_proteins         = get_pi3k_akt_uniprot_ids()
    G_full, G_pathway, found = load_pathway_subgraph(
                                   INTERACTOME, pathway_proteins,
                                   conf_threshold=CONF_THRESHOLD)
    save_pathway_proteins(found, SOURCE_ID, TARGET_ID, OUT_PATHWAY_TXT)

    if G_pathway.number_of_nodes() == 0:
        print(f"[viz] WARNING: 0 nodes. "
              f"Try lowering CONF_THRESHOLD (now {CONF_THRESHOLD}).")
    else:
        if PRUNE_EDGES:
            G_draw = prune_to_top_edges(
                G_pathway, MAX_EDGES,
                must_keep_nodes={SOURCE_ID, TARGET_ID})
        else:
            G_draw = G_pathway
        draw_kegg_style(G_draw, SOURCE_ID, TARGET_ID, CONF_THRESHOLD, OUT_FIG)

    # ── 2. Full interactome DiGraph ────────────────────────────────────────────
    print("\n[full] Building complete interactome graph...")
    G_interactome = load_graph(INTERACTOME)     # ← uses the public wrapper
    save_full_graph(G_interactome, OUT_GRAPHML, OUT_GPICKLE)
    shortest_path_smoke_test(G_interactome)

    print("\n[load_graph] To reload in another script:")
    print("       import pickle")
    print(f"       with open('{OUT_GPICKLE}', 'rb') as f: G = pickle.load(f)")
    print("       path = nx.shortest_path(G, 'P04637', 'P31749', weight='cost')")
    print("\n[load_graph] All done.")