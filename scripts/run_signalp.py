"""
Submit a FASTA file to SignalP 6.0 via the DTU Health Tech webface2 REST API,
poll until complete, and save the raw prediction results to data/output.

DTU fair-use guidelines:
  - Poll at most once every 10 seconds (enforced by POLL_INTERVAL).
  - The service is free for academic use only.
  - See: https://services.healthtech.dtu.dk/services/SignalP-6.0/
"""

import argparse
import os
import re
import sys
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL = "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi"
POLL_INTERVAL = 10  # seconds between status checks
CONFIG_FILE = "/var/www/services/services/SignalP-6.0/webface.cf"
EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")

# The DTU web API exposes two organism classes. Bacterial inputs (gramneg,
# grampos, archaea) map to "Other", which covers all non-eukaryotic lineages.
ORGANISMS = ("eukarya", "other", "gramneg", "grampos", "archaea")
_ORGANISM_MAP = {
    "eukarya": "Eukarya",
    "other": "Other",
    "gramneg": "Other",
    "grampos": "Other",
    "archaea": "Other",
}


def submit_job(fasta_path: Path, organism: str) -> str:
    sequence = fasta_path.read_text()
    # DTU form requires multipart/form-data encoding
    files = {
        "configfile": (None, CONFIG_FILE),
        "fasta":      (None, sequence),
        "organism":   (None, _ORGANISM_MAP[organism]),
        "format":     (None, "short"),
        "mode":       (None, "fast"),
    }
    response = requests.post(BASE_URL, files=files, allow_redirects=True)
    response.raise_for_status()

    # After submission, DTU redirects to ?jobid=XXXX&wait=20 (or embeds it in the body)
    match = re.search(r"jobid=([A-Za-z0-9_\-]+)", response.url + response.text)
    if not match:
        sys.exit(
            "Error: could not parse job ID from DTU response.\n"
            f"Response URL: {response.url}\n"
            f"Response body (first 500 chars):\n{response.text[:500]}"
        )
    job_id = match.group(1)
    print(f"Submitted job: {job_id}")
    return job_id


def poll_until_done(job_id: str) -> str:
    # DTU exposes a JSON status endpoint used by their browser client
    while True:
        response = requests.get(BASE_URL, params={"ajax": "1", "jobid": job_id})
        response.raise_for_status()
        data = response.json()
        status = data.get("status", "unknown")
        print(f"Status: {status}")
        if status not in ("active", "queued"):
            return status
        time.sleep(POLL_INTERVAL)


def fetch_result(job_id: str) -> str:
    # Fetch the completed result page and find the CSV prediction summary link
    response = requests.get(BASE_URL, params={"jobid": job_id, "wait": "0"})
    response.raise_for_status()

    # DTU stores results under a predictable tmp path
    result_url = (
        f"https://services.healthtech.dtu.dk/services/SignalP-6.0/tmp"
        f"/{job_id}/prediction_results.txt"
    )
    result_response = requests.get(result_url)
    if result_response.status_code != 200:
        sys.exit(
            f"Error: could not fetch result file from {result_url}\n"
            f"HTTP status: {result_response.status_code}\n"
            f"Result page (first 500 chars):\n{response.text[:500]}"
        )
    return result_response.text


def main():
    parser = argparse.ArgumentParser(
        description="Submit a FASTA file to SignalP 6.0 and save prediction results."
    )
    parser.add_argument("fasta", type=Path, help="Path to input FASTA file")
    parser.add_argument(
        "--organism",
        choices=ORGANISMS,
        default="other",
        help="Organism type (default: other). Use 'gramneg' or 'grampos' for bacteria.",
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

    job_id = submit_job(args.fasta, args.organism)
    status = poll_until_done(job_id)

    if status != "FINISHED":
        sys.exit(f"Job ended with status '{status}'. No results saved.")

    result = fetch_result(job_id)
    out_file = args.output_dir / f"{args.fasta.stem}_signalp_{job_id}.txt"
    out_file.write_text(result)
    print(f"Results saved to: {out_file}")


if __name__ == "__main__":
    main()
