"""
Submit a PDB file to the FoldSeek web server API, poll until complete,
and save the raw results (tar.gz of per-database M8 alignment files) to data/output.

FoldSeek performs fast structural similarity searches against large
protein structure databases using 3Di structural alphabets.

Fair-use guidelines:
  - Do not flood the server; poll at most once every 10 seconds (POLL_INTERVAL).
  - The server will return status "RATELIMIT" if too many jobs are submitted.
  - See: https://search.foldseek.com
"""

import argparse
import os
import sys
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL = "https://search.foldseek.com/api"
POLL_INTERVAL = 10  # seconds between status checks
EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")

# Search mode: 3diaa (3Di + amino acid, recommended), 3di, tmalign
DEFAULT_MODE = "3diaa"

# Default databases requested: exact path identifiers from the FoldSeek API
DEFAULT_DATABASES = ["afdb50", "afdb-swissprot", "afdb-proteome", "BFVD"]

TERMINAL_STATUSES = {"COMPLETE", "ERROR", "FAILED", "UNKNOWN"}


def submit_job(pdb_path: Path, databases: list[str], mode: str) -> str:
    pdb_content = pdb_path.read_text()
    # Ensure content ends with newline (required by the server)
    if not pdb_content.endswith("\n"):
        pdb_content += "\n"

    # Multiple databases are sent as repeated 'database[]' fields
    data = [("q", pdb_content), ("mode", mode), ("email", EMAIL)]
    for db in databases:
        data.append(("database[]", db))

    response = requests.post(f"{BASE_URL}/ticket", data=data)
    response.raise_for_status()
    result = response.json()

    if result.get("status") == "RATELIMIT":
        sys.exit("Error: FoldSeek rate limit reached. Please wait before resubmitting.")
    if result.get("status") == "MAINTENANCE":
        sys.exit("Error: FoldSeek server is under maintenance. Try again later.")

    ticket = result.get("id")
    if not ticket:
        sys.exit(f"Error: no ticket ID in response: {result}")

    print(f"Submitted job: {ticket}")
    return ticket


def poll_until_done(ticket: str) -> str:
    while True:
        response = requests.get(f"{BASE_URL}/ticket/{ticket}")
        response.raise_for_status()
        status = response.json().get("status", "UNKNOWN")
        print(f"Status: {status}")
        if status in TERMINAL_STATUSES:
            return status
        time.sleep(POLL_INTERVAL)


def fetch_result(ticket: str) -> bytes:
    response = requests.get(f"{BASE_URL}/result/download/{ticket}")
    response.raise_for_status()
    return response.content


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Submit a PDB file to FoldSeek and save raw M8 results (tar.gz). "
            "Searches afdb50, afdb-swissprot, afdb-proteome, and BFVD by default."
        )
    )
    parser.add_argument("pdb", type=Path, help="Path to input PDB file")
    parser.add_argument(
        "--databases",
        nargs="+",
        default=DEFAULT_DATABASES,
        metavar="DB",
        help=(
            "Databases to search (default: afdb50 afdb-swissprot afdb-proteome BFVD). "
            "Available: afdb50, afdb-swissprot, afdb-proteome, BFVD, pdb100, cath50, "
            "mgnify_esm30, gmgcl_id"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["3diaa", "3di", "tmalign"],
        default=DEFAULT_MODE,
        help="Search mode (default: 3diaa — combined 3Di + amino acid)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/output"),
        help="Directory to write results (default: data/output)",
    )
    args = parser.parse_args()

    if not args.pdb.is_file():
        sys.exit(f"Error: PDB file not found: {args.pdb}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    ticket = submit_job(args.pdb, args.databases, args.mode)
    status = poll_until_done(ticket)

    if status != "COMPLETE":
        sys.exit(f"Job ended with status '{status}'. No results saved.")

    result = fetch_result(ticket)
    out_file = args.output_dir / f"{args.pdb.stem}_foldseek_{ticket}.tar.gz"
    out_file.write_bytes(result)
    print(f"Results saved to: {out_file}")
    print("Contains per-database M8 alignment files (alis_<db>.m8).")


if __name__ == "__main__":
    main()
