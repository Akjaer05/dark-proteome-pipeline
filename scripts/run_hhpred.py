"""
Submit a FASTA file to HHpred via the MPI Bioinformatics Toolkit REST API,
poll until complete, and save the raw .hhr results to data/output.

HHpred detects remote structural homologs by HMM-HMM comparison.

MPI Toolkit fair-use guidelines:
  - Do not submit large numbers of jobs in rapid succession.
  - Poll at most once every 10 seconds (enforced by POLL_INTERVAL).
  - See: https://toolkit.tuebingen.mpg.de
"""

import argparse
import os
import sys
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL = "https://toolkit.tuebingen.mpg.de/api/jobs"
POLL_INTERVAL = 10  # seconds between status checks

# Numeric status codes returned by the MPI Toolkit API
STATUS_LABELS = {
    1: "Prepared",
    2: "Queued",
    3: "Running",
    4: "Error",
    5: "Done",
    6: "Submitted",
    7: "Pending",
    9: "Deleted",
}
TERMINAL_STATUSES = {4, 5, 9}  # Error, Done, Deleted

# Available databases: PDB_mmCIF70 (default), PDB_mmCIF30, SCOPe, Pfam, UniRef30
DATABASES = ("PDB_mmCIF70", "PDB_mmCIF30", "SCOPe", "Pfam")


def submit_job(fasta_path: Path, db: str) -> str:
    sequence = fasta_path.read_text()
    payload = {
        "sequence": sequence,
        "db": [db],
    }
    response = requests.post(
        f"{BASE_URL}/?toolName=hhpred",
        json=payload,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    data = response.json()
    if not data.get("successful"):
        sys.exit(f"Error: job submission failed: {data.get('message')}")
    job_id = data["jobID"]
    print(f"Submitted job: {job_id}")
    return job_id


def poll_until_done(job_id: str) -> int:
    while True:
        response = requests.get(f"{BASE_URL}/{job_id}")
        response.raise_for_status()
        data = response.json()
        status_code = data.get("status")
        label = STATUS_LABELS.get(status_code, f"Unknown({status_code})")
        print(f"Status: {label}")
        if status_code in TERMINAL_STATUSES:
            return status_code
        time.sleep(POLL_INTERVAL)


def fetch_result(job_id: str) -> str:
    url = f"{BASE_URL}/{job_id}/results/files/{job_id}.hhr"
    response = requests.get(url)
    response.raise_for_status()
    return response.text


def main():
    parser = argparse.ArgumentParser(
        description="Submit a FASTA file to HHpred and save raw .hhr results."
    )
    parser.add_argument("fasta", type=Path, help="Path to input FASTA file")
    parser.add_argument(
        "--db",
        choices=DATABASES,
        default="PDB_mmCIF70",
        help="Target database (default: PDB_mmCIF70)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory to write results (default: data/output)",
    )
    args = parser.parse_args()

    if not args.fasta.is_file():
        sys.exit(f"Error: FASTA file not found: {args.fasta}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    job_id = submit_job(args.fasta, args.db)
    status_code = poll_until_done(job_id)

    if status_code != 5:
        label = STATUS_LABELS.get(status_code, str(status_code))
        sys.exit(f"Job ended with status '{label}'. No results saved.")

    result = fetch_result(job_id)
    out_file = args.output_dir / f"{args.fasta.stem}_hhpred_{job_id}.hhr"
    out_file.write_text(result)
    print(f"Results saved to: {out_file}")


if __name__ == "__main__":
    main()
