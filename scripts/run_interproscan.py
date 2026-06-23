"""
Submit a FASTA file to EBI InterProScan REST API, poll until complete,
and save the raw JSON result.

EBI fair-use guidelines:
  - Do not submit more than 30 jobs simultaneously.
  - Poll at most once every 3 seconds (enforced by POLL_INTERVAL).
  - Include a meaningful email address in the clientId parameter.
  - See: https://www.ebi.ac.uk/Tools/webservices/help/faq
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
POLL_INTERVAL = 10  # seconds between status checks (EBI asks ≥3 s)
EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")


def submit_job(fasta_path: Path) -> str:
    sequence = fasta_path.read_text()
    payload = {
        "email": EMAIL,
        "title": fasta_path.stem,
        "goterms": "true",
        "pathways": "true",
        "stype": "p",
        "sequence": sequence,
    }
    response = requests.post(f"{BASE_URL}/run", data=payload)
    response.raise_for_status()
    job_id = response.text.strip()
    print(f"Submitted job: {job_id}")
    return job_id


def poll_until_done(job_id: str) -> str:
    while True:
        response = requests.get(f"{BASE_URL}/status/{job_id}")
        response.raise_for_status()
        status = response.text.strip()
        print(f"Status: {status}")
        if status in ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
            return status
        time.sleep(POLL_INTERVAL)


def fetch_result(job_id: str) -> dict:
    response = requests.get(
        f"{BASE_URL}/result/{job_id}/json",
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    return response.json()


def main():
    parser = argparse.ArgumentParser(
        description="Submit a FASTA file to EBI InterProScan and save JSON results."
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

    job_id = submit_job(args.fasta)
    status = poll_until_done(job_id)

    if status != "FINISHED":
        sys.exit(f"Job ended with status '{status}'. No results saved.")

    result = fetch_result(job_id)
    out_file = args.output_dir / f"{args.fasta.stem}_{job_id}.json"
    out_file.write_text(json.dumps(result, indent=2))
    print(f"Results saved to: {out_file}")


if __name__ == "__main__":
    main()
