"""
mapping.py
──────────
Converts UniProt protein IDs to their corresponding gene names using
the UniProt REST API (no API key required).

Output:
    outputs/txt/uniprot_gene_map.txt

Usage:
    # Single protein
    python scripts/mapping.py --proteins P04637

    # Multiple proteins
    python scripts/mapping.py --proteins P04637 P00533 Q9Y243 O14763

    # From a file (one UniProt ID per line)
    python scripts/mapping.py --proteins-file my_proteins.txt
"""

import os
import sys
import argparse
import time
import json

try:
    import requests
except ImportError:
    sys.exit("[ERROR] 'requests' package not found. Install it with: pip install requests")


# ── UniProt API helper ────────────────────────────────────────────────────────

UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"

def fetch_gene_names(uniprot_ids: list, batch_size: int = 100) -> dict:
    """
    Query the UniProt REST API to get gene names for a list of UniProt IDs.

    Returns a dict: uniprot_id → gene_name (or None if not found).
    Processes IDs in batches to avoid URL-length limits.
    """
    mapping = {uid: None for uid in uniprot_ids}

    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i : i + batch_size]
        query = " OR ".join(f"accession:{uid}" for uid in batch)

        params = {
            "query":  query,
            "fields": "accession,gene_names",
            "format": "json",
            "size":   batch_size,
        }

        print(f"[INFO] Querying UniProt for batch {i//batch_size + 1} "
              f"({len(batch)} proteins) ...")

        try:
            resp = requests.get(f"{UNIPROT_API_BASE}/search",
                                params=params, timeout=30)
            resp.raise_for_status()
        except requests.RequestException as exc:
            print(f"[WARN] UniProt API request failed: {exc}")
            print(f"       Proteins in this batch will be marked as 'API_ERROR'.")
            for uid in batch:
                mapping[uid] = "API_ERROR"
            continue

        data = resp.json()

        for entry in data.get("results", []):
            accession = entry.get("primaryAccession", "")
            genes     = entry.get("genes", [])

            # Pick the first recommended gene name, fall back to any synonym
            gene_name = None
            for g in genes:
                gn = g.get("geneName", {})
                if gn.get("value"):
                    gene_name = gn["value"]
                    break
            if gene_name is None:
                for g in genes:
                    for syns in g.get("synonyms", []):
                        if syns.get("value"):
                            gene_name = syns["value"]
                            break

            if accession in mapping:
                mapping[accession] = gene_name or "NO_GENE_NAME"

        # Mark proteins that returned no result
        for uid in batch:
            if mapping[uid] is None:
                mapping[uid] = "NOT_FOUND"

        time.sleep(0.25)   # be polite to the UniProt servers

    return mapping


# ── Text report ───────────────────────────────────────────────────────────────

def write_mapping_report(mapping: dict, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    found     = sum(1 for v in mapping.values()
                    if v not in (None, "NOT_FOUND", "API_ERROR", "NO_GENE_NAME"))
    not_found = len(mapping) - found

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("=" * 55 + "\n")
        fh.write("  UNIPROT ID -> GENE NAME CONVERSION MAP\n")
        fh.write("=" * 55 + "\n")
        fh.write(f"  Total proteins queried : {len(mapping)}\n")
        fh.write(f"  Gene names found       : {found}\n")
        fh.write(f"  Not found / errors     : {not_found}\n")
        fh.write("=" * 55 + "\n\n")

        fh.write(f"  {'UniProt ID':<20}  Gene Name\n")
        fh.write(f"  {'-'*20}  {'-'*30}\n")
        for uid, gene in sorted(mapping.items()):
            fh.write(f"  {uid:<20}  {gene or 'NOT_FOUND'}\n")

    print(f"[INFO] Mapping saved → {out_path}")


# ── Pretty print to console ───────────────────────────────────────────────────

def print_mapping(mapping: dict) -> None:
    print("\n-- UniProt -> Gene Name ------------------------------------------")
    print(f"  {'UniProt ID':<20}  Gene Name")
    print(f"  {'-'*20}  {'-'*30}")
    for uid, gene in sorted(mapping.items()):
        print(f"  {uid:<20}  {gene or 'NOT_FOUND'}")
    print()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert UniProt IDs to gene names via the UniProt REST API.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--proteins",      nargs="+",
                       help="One or more UniProt IDs")
    group.add_argument("--proteins-file", metavar="FILE",
                       help="File with one UniProt ID per line")
    parser.add_argument("--out-txt", default="outputs/txt/uniprot_gene_map.txt")
    args = parser.parse_args()

    if args.proteins:
        proteins = args.proteins
    else:
        if not os.path.isfile(args.proteins_file):
            sys.exit(f"[ERROR] File not found: {args.proteins_file}")
        with open(args.proteins_file) as fh:
            proteins = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

    print(f"[INFO] Proteins to map: {proteins}")

    mapping = fetch_gene_names(proteins)
    print_mapping(mapping)
    write_mapping_report(mapping, args.out_txt)

    print("[DONE] mapping.py finished.")


if __name__ == "__main__":
    main()