"""
mapping.py
──────────────────────────────
UniProt → Gene Name mapper
Uses the UniProt REST API to retrieve:
    - Gene name
    - Protein name
    - Organism name
"""

import os
import time
import requests

# Output folder and UniProt API endpoint
TXT_DIR     = "outputs/txt"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"

os.makedirs(TXT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════
# Internal helper — Query UniProt API
# ═══════════════════════════════════════════════════
def _query_uniprot(uniprot_id: str) -> dict:
    """
    Query UniProt for a single UniProt accession ID.

    Returns a dictionary containing:
        - gene_name
        - protein_name
        - organism
    """

    # API query parameters
    params = {
        "query": f"accession:{uniprot_id}",
        "fields": "accession,gene_names,protein_name,organism_name",
        "format": "json",
        "size": 1,
    }

    # Retry failed requests up to 3 times
    for attempt in range(3):

        try:
            response = requests.get(
                UNIPROT_URL,
                params=params,
                timeout=30
            )

            response.raise_for_status()
            data = response.json()

            # No matching protein found
            if not data.get("results"):

                return {
                    "uniprot_id": uniprot_id,
                    "gene_name": "NOT FOUND",
                    "protein_name": "NOT FOUND",
                    "organism": "NOT FOUND",
                }

            result = data["results"][0]

            gene_names  = result.get("genes", [])
            protein_desc = result.get("proteinDescription", {})
            organism     = result.get("organism", {})

            # Extract gene name
            gene_name = "N/A"

            if gene_names:
                gene_name = (
                    gene_names[0]
                    .get("geneName", {})
                    .get("value", "N/A")
                )

            # Extract protein name
            protein_name = "N/A"

            rec = protein_desc.get("recommendedName", {})

            if rec:
                protein_name = (
                    rec.get("fullName", {})
                    .get("value", "N/A")
                )

            return {
                "uniprot_id": uniprot_id,
                "gene_name": gene_name,
                "protein_name": protein_name,
                "organism": organism.get("scientificName", "N/A"),
            }

        # Retry if request fails
        except requests.exceptions.RequestException:

            print(
                f"    WARN retry {attempt+1}/3 "
                f"failed for {uniprot_id}"
            )

            time.sleep(2)

    # Returned if all retries fail
    return {
        "uniprot_id": uniprot_id,
        "gene_name": "API ERROR",
        "protein_name": "API ERROR",
        "organism": "API ERROR",
    }


# ═══════════════════════════════════════════════════
# Public mapping function
# ═══════════════════════════════════════════════════
def map_uniprot_to_gene(uniprot_ids: list) -> dict:
    """
    Map a list of UniProt IDs to gene names.

    Also saves a text report containing:
        - UniProt ID
        - Gene name
        - Protein name
    """

    print(f"\n[mapping] Mapping {len(uniprot_ids)} UniProt IDs...\n")

    results = []
    mapping = {}

    # Query each UniProt ID
    for i, uid in enumerate(uniprot_ids, 1):

        print(f"[{i}/{len(uniprot_ids)}] {uid}", end=" ... ")

        info = _query_uniprot(uid)

        results.append(info)

        # Store simple ID → gene mapping
        mapping[uid] = info["gene_name"]

        print(info["gene_name"])

        # Gentle API rate limiting
        time.sleep(0.2)

    # ── Save mapping report ─────────────────────────
    out_txt = os.path.join(
        TXT_DIR,
        "uniprot_gene_map.txt"
    )

    with open(out_txt, "w", encoding="utf-8") as f:

        f.write("UniProt ID to Gene Name Mapping\n")
        f.write("=" * 60 + "\n")

        f.write(
            f"{'UniProt ID':<15}"
            f"{'Gene Name':<15}"
            f"{'Protein Name'}\n"
        )

        f.write("-" * 60 + "\n")

        for r in results:

            f.write(
                f"{r['uniprot_id']:<15}"
                f"{r['gene_name']:<15}"
                f"{r['protein_name']}\n"
            )

    print(f"\nOK Saved -> {out_txt}")

    return mapping