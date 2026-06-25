"""
Submit a FASTA file to EBI NCBI BLAST (blastp, UniProt Swiss-Prot) via the EBI
Job Dispatcher REST API, poll until complete, and save the raw XML result.

Uses the same EBI infrastructure as run_interproscan.py and run_phobius.py —
more reliable from cloud servers than direct NCBI QBLAST.

EBI fair-use guidelines:
  - Do not submit more than 30 jobs simultaneously.
  - Poll at most once every 3 seconds (enforced by POLL_INTERVAL).
  - Include a meaningful email address in all requests.
  - See: https://www.ebi.ac.uk/Tools/webservices/help/faq
"""

import argparse
import os
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

BASE_URL     = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
POLL_INTERVAL = 15   # seconds between status checks
EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")


def submit_job(fasta_path: Path, database: str, evalue: str) -> str:
    sequence = fasta_path.read_text()
    payload = {
        "email":    EMAIL,
        "program":  "blastp",
        "database": database,
        "sequence": sequence,
        "stype":    "protein",
        "scores":   10,
        "alignments": 10,
        "exp":      str(evalue),
    }
    for attempt in range(3):
        try:
            r = requests.post(f"{BASE_URL}/run", data=payload, timeout=60)
            r.raise_for_status()
            jid = r.text.strip()
            if jid:
                print(f"Submitted job: {jid}")
                return jid
        except Exception as exc:
            if attempt < 2:
                print(f"Submission attempt {attempt+1} failed: {exc} — retrying in 10s")
                time.sleep(10)
            else:
                sys.exit(f"Error: submission failed after 3 attempts: {exc}")
    sys.exit("Error: no job ID in EBI response.")


def poll_until_done(job_id: str, timeout: int = 720) -> str:
    print(f"Waiting 10s before first status check...")
    time.sleep(10)
    start    = time.time()
    deadline = start + timeout
    poll_n   = 0
    status   = "PENDING"
    while time.time() < deadline:
        poll_n += 1
        elapsed = int(time.time() - start)
        try:
            r = requests.get(f"{BASE_URL}/status/{job_id}", timeout=30)
            r.raise_for_status()
            status = r.text.strip()
            print(f"Poll #{poll_n} at {elapsed}s: {status}")
            if status == "FINISHED":
                return status
            if status in ("FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
                sys.exit(f"Job ended with status '{status}'.")
        except requests.exceptions.Timeout:
            print(f"Poll #{poll_n} timed out (30s), retrying")
        remaining = deadline - time.time()
        if remaining > 0:
            time.sleep(min(POLL_INTERVAL, remaining))
    sys.exit(f"Job did not complete within {timeout}s. Last status: {status}")


def fetch_result(job_id: str) -> str:
    r = requests.get(f"{BASE_URL}/result/{job_id}/xml", timeout=60)
    r.raise_for_status()
    return r.text


def summarise(xml_text: str) -> int:
    """Print a brief summary of hits and return the hit count."""
    try:
        import re as _re
        xml_clean = _re.sub(r'\s+xmlns(?::\w+)?="[^"]+"', '', xml_text)
        xml_clean = _re.sub(r'\s+xsi:\w+="[^"]+"', '', xml_clean)
        root = ET.fromstring(xml_clean)
        hits = list(root.iter("hit"))
        print(f"\n{len(hits)} hit(s) found:")
        for hit in hits[:5]:
            acc  = hit.get("ac") or hit.get("id") or "?"
            desc = (hit.get("description") or "")[:70]
            aln  = hit.find(".//alignment")
            ev   = f"  E={aln.findtext('expectation','?')}" if aln is not None else ""
            print(f"  {acc:<14} {desc}{ev}")
        if len(hits) > 5:
            print(f"  ... and {len(hits)-5} more")
        return len(hits)
    except Exception as exc:
        print(f"Warning: could not parse result XML: {exc}")
        return 0


def main():
    parser = argparse.ArgumentParser(
        description="Submit a FASTA to EBI NCBI BLAST (blastp) and save XML results."
    )
    parser.add_argument("fasta", type=Path, help="Path to input FASTA file")
    parser.add_argument(
        "--database", default="uniprotkb_bacteria",
        help="Database to search (default: uniprotkb_bacteria)"
    )
    # EBI ncbiblast only accepts exp in scientific notation (e.g. "1e-3", not "0.001")
    parser.add_argument(
        "--evalue", default="1e-3",
        help="E-value threshold in scientific notation (default: 1e-3)"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/output"),
        help="Directory to write XML results (default: data/output)",
    )
    args = parser.parse_args()

    if not args.fasta.is_file():
        sys.exit(f"Error: FASTA file not found: {args.fasta}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    job_id = submit_job(args.fasta, args.database, args.evalue)
    poll_until_done(job_id)

    result = fetch_result(job_id)
    out_file = args.output_dir / f"{args.fasta.stem}_blast_{job_id}.xml"
    out_file.write_text(result)
    print(f"\nResults saved to: {out_file}")
    summarise(result)


if __name__ == "__main__":
    main()
