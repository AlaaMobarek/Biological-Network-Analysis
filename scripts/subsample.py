import networkx as nx

# ── Load full interactome ──────────────────────────────────────────────────────
G_full = nx.DiGraph()

with open("data/PathLinker_2018_human-ppi-weighted-cap0_75.txt", "r") as f:
    for line in f:
        if line.startswith("#"):        # skip header
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:             # skip malformed lines
            continue
        tail, head, weight = parts[0], parts[1], float(parts[2])
        method = parts[3] if len(parts) > 3 else ""
        G_full.add_edge(tail, head, weight=weight, method=method)

print(f"Full graph: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges()} edges")

# ── Define your seed proteins (p53 pathway) ───────────────────────────────────
seed_proteins = {"P04637", "P06400", "Q09472", "P38936", "O15151", "Q8NI29"}

# ── Expand to include direct neighbors of seeds ───────────────────────────────
pathway_nodes = set(seed_proteins)

for protein in seed_proteins:
    if protein in G_full:
        # Add all direct neighbors (in and out)
        pathway_nodes.update(G_full.predecessors(protein))
        pathway_nodes.update(G_full.successors(protein))

# ── Extract subgraph ───────────────────────────────────────────────────────────
G_sub = G_full.subgraph(pathway_nodes).copy()

print(f"Subgraph: {G_sub.number_of_nodes()} nodes, {G_sub.number_of_edges()} edges")

# ── Save subgraph to a new file ───────────────────────────────────────────────
with open("data/interactome.txt", "w") as f:
    f.write("#tail\thead\tedge_weight\tedge_type\n")
    for tail, head, data in G_sub.edges(data=True):
        f.write(f"{tail}\t{head}\t{data['weight']:.6e}\t{data.get('method','')}\n")

print("Subgraph saved to data/interactome.txt")