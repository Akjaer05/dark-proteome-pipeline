"""
Parse an InterProScan JSON result file into a flat TSV summary table.
One row per match location.
"""

import argparse
import csv
import json
import sys
from pathlib import Path


COLUMNS = [
    "protein_id",
    "database",
    "signature_accession",
    "signature_name",
    "signature_description",
    "type",
    "interpro_accession",
    "interpro_description",
    "start",
    "end",
    "evalue",
    "go_terms",
]


def parse_result(result: dict) -> list[dict]:
    protein_id = result["xref"][0]["id"] if result.get("xref") else "unknown"
    rows = []

    for match in result.get("matches", []):
        sig = match["signature"]
        entry = sig.get("entry") or {}
        lib = sig.get("signatureLibraryRelease", {}).get("library", "")

        go_terms = "|".join(
            f"{g['id']}({g['name']})" for g in entry.get("goXRefs", [])
        )

        evalue = match.get("evalue", "")

        for loc in match.get("locations", []):
            rows.append({
                "protein_id": protein_id,
                "database": lib,
                "signature_accession": sig.get("accession", ""),
                "signature_name": sig.get("name", ""),
                "signature_description": sig.get("description") or "",
                "type": sig.get("type", ""),
                "interpro_accession": entry.get("accession", ""),
                "interpro_description": entry.get("description") or "",
                "start": loc.get("start", ""),
                "end": loc.get("end", ""),
                "evalue": evalue,
                "go_terms": go_terms,
            })

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Parse InterProScan JSON into a TSV summary table."
    )
    parser.add_argument("json_file", type=Path, help="InterProScan JSON result file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory to write TSV (default: data/output)",
    )
    args = parser.parse_args()

    if not args.json_file.is_file():
        sys.exit(f"Error: file not found: {args.json_file}")

    data = json.loads(args.json_file.read_text())
    results = data.get("results", [])

    rows = []
    for result in results:
        rows.extend(parse_result(result))

    out_file = args.output_dir / (args.json_file.stem + "_summary.tsv")
    with out_file.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_file}")


if __name__ == "__main__":
    main()
