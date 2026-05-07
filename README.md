# PPI Network Analysis
## Biological Network Analysis using NetworkX

Analyzes the **PathLinker 2018 Human PPI Weighted** interactome.

---

## Project Structure

```
project/
├── data/
│   └── interactome.txt          ← download from PathLinker
├── scripts/
│   ├── load_graph.py            ← graph construction (shared module)
│   ├── shortest_path.py         ← requirement 1: shortest paths
│   ├── neighbors.py             ← requirement 2: neighbor listing
│   ├── degree_analysis.py       ← requirement 3: degree histogram + ranking
│   ├── mapping.py               ← requirement 4: UniProt → gene name
│   └── adjacency_matrix.py      ← requirement 5: unweighted adjacency matrix
├── outputs/
│   ├── txt/
│   ├── figures/
│   └── matrices/
├── requirements.txt
└── README.md
```

---

## Setup

```bash
# 1. Download the interactome
#    PathLinker_2018_human-ppi-weighted-cap0_75.txt → rename to data/interactome.txt

# 2. Install dependencies
pip install -r requirements.txt

# 3. Create output directories
mkdir -p outputs/txt outputs/figures outputs/matrices
```

---

## Usage

### Requirement 1 — Shortest Paths between two proteins
```bash
python scripts/shortest_path.py \
    --source  P04637 \
    --target  P00533

# Output:
#   outputs/txt/shortest_paths.txt     ← all paths with weights & confidences
#   outputs/figures/shortest_path_subnetwork.png
```

### Requirement 2 — Neighbors of a protein
```bash
python scripts/neighbors.py --protein P04637

# Output:
#   outputs/txt/neighbors.txt   ← degree, all in/out neighbors with weights
```

### Requirement 3 — Degree analysis for a set of proteins
```bash
python scripts/degree_analysis.py \
    --proteins P04637 P00533 Q9Y243 O14763

# Output:
#   outputs/txt/degree_ranking.txt        ← ranked from most to least connected
#   outputs/figures/degree_histogram.png  ← per-protein bar chart + ranking chart
```

### Requirement 4 — UniProt ID → Gene Name mapping
```bash
# Single protein
python scripts/mapping.py --proteins P04637

# Multiple proteins
python scripts/mapping.py --proteins P04637 P00533 Q9Y243 O14763

# Output:
#   outputs/txt/uniprot_gene_map.txt
```

### Requirement 5 — Unweighted Adjacency Matrix
```bash
# Sub-matrix for specific proteins (recommended)
python scripts/adjacency_matrix.py \
    --proteins P04637 P00533 Q9Y243 O14763

# Full graph – sparse COO format (safe for large graphs)
python scripts/adjacency_matrix.py --full --sparse

# Output:
#   outputs/matrices/adjacency_matrix.csv
#   outputs/figures/graph_overview.png
```

---

## Interactome File Format

| Column | Content |
|--------|---------|
| 0 | Tail / source protein (UniProt ID) |
| 1 | Head / destination protein (UniProt ID) |
| 2 | Interaction confidence (0–1) |
| 3 | Detection method |

Edge weight stored as **1 / confidence** so Dijkstra's algorithm treats
high-confidence interactions as "shorter" (more probable) paths.

---

## Notes

- The interactome is a **directed** graph; in/out degrees are tracked separately.
- `mapping.py` uses the **UniProt REST API** (no key needed); requires internet access.
- For the full interactome (~50k+ nodes), use `--sparse` with `adjacency_matrix.py`.