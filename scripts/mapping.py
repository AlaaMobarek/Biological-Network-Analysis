"""
mapping.py
──────────
Member 4 — UniProt ID to Gene Name Mapping

Responsibilities:
    - Given one or more UniProt IDs, retrieve their gene names
    - Uses the UniProt REST API (no manual database needed)
    - Saves the conversion map to a text file

Usage (standalone):
    python scripts/mapping.py

Usage (from main.py):
    from scripts.mapping import map_uniprot_to_gene
    map_uniprot_to_gene(["P04637", "P00533", "P06400"])
"""

import os
import time
import requests

# ── Constants ──────────────────────────────────────────────────────────────────
TXT_DIR     = "Code/outputs/txt"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
os.makedirs(TXT_DIR, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Query UniProt API for a single protein
# ══════════════════════════════════════════════════════════════════════════════
def _query_uniprot(uniprot_id: str) -> dict:
    """
    Queries the UniProt REST API for a single UniProt accession ID
    and returns gene name, protein name, and organism.

    Parameters
    ----------
    uniprot_id : str — UniProt accession (e.g. "P04637")

    Returns
    -------
    dict with keys: uniprot_id, gene_name, protein_name, organism
    """
    params = {
        "query" : f"accession:{uniprot_id}",
        "fields": "accession,gene_names,protein_name,organism_name",
        "format": "json",
        "size"  : 1,
    }

    try:
        response = requests.get(UNIPROT_URL, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()

        if not data.get("results"):
            return {
                "uniprot_id"   : uniprot_id,
                "gene_name"    : "NOT FOUND",
                "protein_name" : "NOT FOUND",
                "organism"     : "NOT FOUND",
            }

        result       = data["results"][0]
        gene_names   = result.get("genes", [])
        protein_desc = result.get("proteinDescription", {})
        organism     = result.get("organism", {})

        # Extract primary gene name
        gene_name = "N/A"
        if gene_names:
            primary = gene_names[0].get("geneName", {})
            gene_name = primary.get("value", "N/A")

        # Extract recommended protein name
        protein_name = "N/A"
        rec_name = protein_desc.get("recommendedName", {})
        if rec_name:
            full_name = rec_name.get("fullName", {})
            protein_name = full_name.get("value", "N/A")
        else:
            # Fall back to submitted name
            sub_names = protein_desc.get("submissionNames", [])
            if sub_names:
                protein_name = sub_names[0].get("fullName", {}).get("value", "N/A")

        organism_name = organism.get("scientificName", "N/A")

        return {
            "uniprot_id"   : uniprot_id,
            "gene_name"    : gene_name,
            "protein_name" : protein_name,
            "organism"     : organism_name,
        }

    except requests.exceptions.RequestException as e:
        print(f"    ⚠ API error for {uniprot_id}: {e}")
        return {
            "uniprot_id"   : uniprot_id,
            "gene_name"    : "API ERROR",
            "protein_name" : "API ERROR",
            "organism"     : "API ERROR",
        }


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Map a list of UniProt IDs to gene names
# ══════════════════════════════════════════════════════════════════════════════
def map_uniprot_to_gene(uniprot_ids: list) -> dict:
    """
    Converts a list of UniProt IDs to gene names using the UniProt REST API.
    Saves the full mapping table to a text file.

    Parameters
    ----------
    uniprot_ids : list — list of UniProt accession strings

    Returns
    -------
    mapping : dict — { uniprot_id : gene_name }
    """
    print(f"\n[mapping] Mapping {len(uniprot_ids)} UniProt IDs to gene names...")
    print(f"  (Querying UniProt REST API — this may take a few seconds)\n")

    results = []
    mapping = {}

    for i, uid in enumerate(uniprot_ids, 1):
        print(f"  [{i}/{len(uniprot_ids)}] Querying: {uid}", end=" ... ")
        info = _query_uniprot(uid)
        results.append(info)
        mapping[uid] = info["gene_name"]
        print(f"→ {info['gene_name']}")

        # Be polite to the API — small delay between requests
        time.sleep(0.3)

    # ── Save to text file ─────────────────────────────────────────────────────
    out_txt = os.path.join(TXT_DIR, "uniprot_gene_map.txt")
    with open(out_txt, "w") as f:
        f.write(f"UniProt ID to Gene Name Mapping\n")
        f.write(f"{'='*75}\n")
        f.write(f"{'UniProt ID':<14} {'Gene Name':<14} "
                f"{'Protein Name':<35} {'Organism'}\n")
        f.write(f"{'─'*75}\n")
        for r in results:
            f.write(f"{r['uniprot_id']:<14} {r['gene_name']:<14} "
                    f"{r['protein_name']:<35} {r['organism']}\n")

    print(f"\n  Mapping saved → {out_txt}")

    # ── Print summary table ───────────────────────────────────────────────────
    print(f"\n  {'UniProt ID':<14} {'Gene Name':<14} {'Protein Name'}")
    print(f"  {'─'*60}")
    for r in results:
        print(f"  {r['uniprot_id']:<14} {r['gene_name']:<14} "
              f"{r['protein_name'][:35]}")

    print()
    return mapping


# ══════════════════════════════════════════════════════════════════════════════
# MAIN — runs when executed directly
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    # Change this list to any UniProt IDs you want to map
    uniprot_ids = [
        "P04637",   # TP53
        "P00533",   # EGFR
        "P06400",   # RB1
        "Q09472",   # EP300
        "P38936",   # CDKN1A
        "P12931",   # SRC
        "P27986",   # PIK3R1
        "P62993",   # GRB2
        "O14543",   # SOCS3
        "P43403",   # ZAP70
    ]

    mapping = map_uniprot_to_gene(uniprot_ids)
    print("✅ mapping.py complete!\n")