import networkx as nx

G = nx.DiGraph()

with open("interactome.txt") as f:
    for line in f:
        if line.startswith("#"):  # Skip header
            continue
        parts = line.strip().split("\t")
        tail, head, weight, method = parts
        # source, target, weight, method = line.strip().split()
        
        G.add_edge(
            source = tail,
            target = head,
            weight=float(weight),
            method=method
        )