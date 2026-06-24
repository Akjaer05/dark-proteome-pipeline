"""
Submit a FASTA file to NCBI BLASTp (nr database) via the NCBI REST API,
poll until complete, and save the raw JSON result.

NCBI fair-use guidelines:
  - Do not submit more than one request every 10 seconds (enforced by POLL_INTERVAL).
  - Always include email and tool parameters so NCBI can contact you if needed.
  - Do not submit more than 20 searches at a time on weekdays 06:00–24:00 EST.
  - See: https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html
"""

import argparse
import io
import json
import os
import sys
import time
import zipfile
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
POLL_INTERVAL = 10  # seconds between status checks (NCBI asks ≥10 s)
TOOL = "dark-proteome-pipeline"
EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")


def submit_job(fasta_path: Path) -> tuple[str, int]:
    sequence = fasta_path.read_text()
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "HITLIST_SIZE": 10,
        "EMAIL": EMAIL,
        "TOOL": TOOL,
    }
    response = requests.post(BASE_URL, data=params)
    response.raise_for_status()

    rid, rtoe = None, 30
    for line in response.text.splitlines():
        if line.startswith("    RID = "):
            rid = line.split("=", 1)[1].strip()
        elif line.startswith("    RTOE = "):
            rtoe = int(line.split("=", 1)[1].strip())

    if not rid:
        sys.exit("Error: could not parse RID from NCBI response.")

    print(f"Submitted job: RID={rid}, estimated wait={rtoe}s")
    return rid, rtoe


def poll_until_done(rid: str, rtoe: int) -> str:
    # Wait the server's estimated time before first poll
    print(f"Waiting {rtoe}s before first status check...")
    time.sleep(rtoe)

    while True:
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
            "EMAIL": EMAIL,
            "TOOL": TOOL,
        }
        response = requests.get(BASE_URL, params=params)
        response.raise_for_status()

        status = "UNKNOWN"
        for line in response.text.splitlines():
            if "Status=" in line:
                status = line.strip().split("=", 1)[1].strip()
                break

        print(f"Status: {status}")
        if status in ("READY", "FAILED", "UNKNOWN"):
            return status
        time.sleep(POLL_INTERVAL)


def fetch_result(rid: str) -> dict:
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "JSON2",
        "DESCRIPTIONS": 10,
        "ALIGNMENTS": 10,
        "EMAIL": EMAIL,
        "TOOL": TOOL,
    }
    response = requests.get(BASE_URL, params=params)
    response.raise_for_status()

    # NCBI returns JSON2 results as a ZIP archive containing one or more JSON files.
    # The index file (<RID>.json) lists the actual result files (<RID>_1.json, etc.).
    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        result_files = [n for n in zf.namelist() if n.endswith(".json") and "_" in n]
        if not result_files:
            sys.exit("Error: no result files found inside NCBI ZIP response.")
        return json.loads(zf.read(result_files[0]))


def main():
    parser = argparse.ArgumentParser(
        description="Submit a FASTA file to NCBI BLASTp (nr) and save JSON results."
    )
    parser.add_argument("fasta", type=Path, help="Path to input FASTA file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory to write JSON results (default: data/output)",
    )
    args = parser.parse_args()

    if not args.fasta.is_file():
        sys.exit(f"Error: FASTA file not found: {args.fasta}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    rid, rtoe = submit_job(args.fasta)
    status = poll_until_done(rid, rtoe)

    if status != "READY":
        sys.exit(f"Job ended with status '{status}'. No results saved.")

    result = fetch_result(rid)
    out_file = args.output_dir / f"{args.fasta.stem}_blast_{rid}.json"
    out_file.write_text(json.dumps(result, indent=2))
    print(f"Results saved to: {out_file}")


if __name__ == "__main__":
    main()
