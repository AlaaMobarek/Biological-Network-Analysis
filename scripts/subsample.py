# import networkx as nx

# # ── Load full interactome ──────────────────────────────────────────────────────
# G_full = nx.DiGraph()

# with open("data/PathLinker_2018_human-ppi-weighted-cap0_75.txt", "r") as f:
#     for line in f:
#         if line.startswith("#"):        # skip header
#             continue
#         parts = line.strip().split("\t")
#         if len(parts) < 3:             # skip malformed lines
#             continue
#         tail, head, weight = parts[0], parts[1], float(parts[2])
#         method = parts[3] if len(parts) > 3 else ""
#         G_full.add_edge(tail, head, weight=weight, method=method)

# print(f"Full graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

# # ── Define your seed proteins (p53 pathway) ───────────────────────────────────
# seed_proteins = {"P04637", "P06400", "Q09472", "P38936", "O15151", "Q8NI29"}

# # ── Expand to include direct neighbors of seeds ───────────────────────────────
# pathway_nodes = set(seed_proteins)

# for protein in seed_proteins:
#     if protein in G_full:
#         # Add all direct neighbors (in and out)
#         pathway_nodes.update(G_full.predecessors(protein))
#         pathway_nodes.update(G_full.successors(protein))

# # ── Extract subgraph ───────────────────────────────────────────────────────────
# G_sub = G_full.subgraph(pathway_nodes).copy()

# print(f"Subgraph: {G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")

# # ── Save subgraph to a new file ───────────────────────────────────────────────
# with open("data/interactome.txt", "w") as f:
#     f.write("#tail\thead\tedge_weight\tedge_type\n")
#     for tail, head, data in G_sub.edges(data=True):
#         f.write(f"{tail}\t{head}\t{data['weight']:.6e}\t{data.get('method','')}\n")

# print("Subgraph saved to data/interactome.txt")
import networkx as nx

# ── Load full interactome ──────────────────────────────────────────────────────
G_full = nx.DiGraph()

with open("data/PathLinker_2018_human-ppi-weighted-cap0_75.txt", "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        tail, head, weight = parts[0], parts[1], float(parts[2])
        method = parts[3] if len(parts) > 3 else ""
        G_full.add_edge(tail, head, weight=weight, method=method)

print(f"Full graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

# ── TP53 → DNA Repair Pathway Seed Proteins ───────────────────────────────────
seed_proteins = {
    "P04637",   # TP53    — master regulator
    "Q13315",   # ATM     — detects DSBs, activates TP53
    "Q13535",   # ATR     — detects SSBs / replication stress
    "O96017",   # CHEK2   — phosphorylates TP53
    "P38936",   # CDKN1A  — p21, CDK inhibitor
    "Q00987",   # MDM2    — negative regulator of TP53
    "Q09472",   # EP300   — acetylates and stabilises TP53
    "P24522",   # GADD45A — NER / BER effector
    "P38398",   # BRCA1   — homologous recombination repair
    "Q92466",   # DDB2    — nucleotide excision repair
    "Q07812",   # BAX     — pro-apoptotic
    "Q9BXH1",   # PUMA    — p53-upregulated apoptosis
}

# ── Expand to 1-hop neighbors ─────────────────────────────────────────────────
present = {p for p in seed_proteins if p in G_full}
pathway_nodes = set(present)
for protein in present:
    pathway_nodes.update(G_full.predecessors(protein))
    pathway_nodes.update(G_full.successors(protein))

G_sub = G_full.subgraph(pathway_nodes).copy()
print(f"Subgraph: {G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")

# ── Save — exact same shape as the original file ──────────────────────────────
with open("data/interactome.txt", "w") as f:
    f.write("#tail\thead\tedge_weight\tedge_type\n")          # header row
    for tail, head, data in G_sub.edges(data=True):
        f.write(
            f"{tail}\t"
            f"{head}\t"
            f"{data['weight']:.6e}\t"                         # e.g. 7.500000e-01
            f"{data.get('method', '')}\n"                     # e.g. MI:0006 (anti…)
        )

print("Saved → data/interactome.txt")