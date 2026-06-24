"""
Master annotation pipeline for dark proteome proteins.

Runs all available tools in parallel and saves results to a named
subfolder under data/output/<protein_name>/.

  FASTA input -> InterProScan, NCBI BLASTp, Phobius
  PDB input   -> FoldSeek
  Both        -> all four tools simultaneously
"""

import argparse
import concurrent.futures
import io
import json
import os
import sys
import threading
import time
import zipfile
from datetime import datetime
from pathlib import Path

import requests
from dotenv import load_dotenv

load_dotenv()

EMAIL = os.environ.get("EBI_EMAIL")
if not EMAIL:
    sys.exit("Error: EBI_EMAIL is not set. Add it to your .env file.")

# ── Thread-safe logging ───────────────────────────────────────────────────────

_lock = threading.Lock()


def log(tag: str, msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    with _lock:
        print(f"[{ts}] [{tag:<13}] {msg}", flush=True)


# ── InterProScan ──────────────────────────────────────────────────────────────

_IPR = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"


def run_interproscan(fasta_path: Path, out_dir: Path, name: str) -> Path:
    TAG = "InterProScan"
    log(TAG, "Submitting...")

    r = requests.post(f"{_IPR}/run", data={
        "email": EMAIL, "title": name, "goterms": "true",
        "pathways": "true", "stype": "p",
        "sequence": fasta_path.read_text(),
    })
    r.raise_for_status()
    job_id = r.text.strip()
    log(TAG, f"Job: {job_id}")

    while True:
        status = requests.get(f"{_IPR}/status/{job_id}").text.strip()
        log(TAG, status)
        if status in ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
            break
        time.sleep(10)

    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")

    result = requests.get(
        f"{_IPR}/result/{job_id}/json", headers={"Accept": "application/json"}
    ).json()
    out = out_dir / f"interproscan_{job_id}.json"
    out.write_text(json.dumps(result, indent=2))
    log(TAG, f"Saved -> {out.name}")
    return out


# ── NCBI BLASTp ───────────────────────────────────────────────────────────────

_BLAST = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
_BLAST_TOOL = "dark-proteome-pipeline"


def run_blast(fasta_path: Path, out_dir: Path, name: str) -> Path:
    TAG = "BLASTp"
    log(TAG, "Submitting...")

    r = requests.post(_BLAST, data={
        "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
        "QUERY": fasta_path.read_text(), "FORMAT_TYPE": "JSON2",
        "HITLIST_SIZE": 10, "EMAIL": EMAIL, "TOOL": _BLAST_TOOL,
    })
    r.raise_for_status()

    rid, rtoe = None, 30
    for line in r.text.splitlines():
        if line.startswith("    RID = "):
            rid = line.split("=", 1)[1].strip()
        elif line.startswith("    RTOE = "):
            rtoe = int(line.split("=", 1)[1].strip())

    if not rid:
        raise RuntimeError("Could not parse RID from NCBI response")

    log(TAG, f"Job: {rid} (est. wait: {rtoe}s)")
    time.sleep(rtoe)

    while True:
        r = requests.get(_BLAST, params={
            "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo",
            "EMAIL": EMAIL, "TOOL": _BLAST_TOOL,
        })
        status = "UNKNOWN"
        for line in r.text.splitlines():
            if "Status=" in line:
                status = line.strip().split("=", 1)[1].strip()
                break
        log(TAG, status)
        if status in ("READY", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)

    if status != "READY":
        raise RuntimeError(f"Job ended with status '{status}'")

    r = requests.get(_BLAST, params={
        "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2",
        "DESCRIPTIONS": 10, "ALIGNMENTS": 10,
        "EMAIL": EMAIL, "TOOL": _BLAST_TOOL,
    })
    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        result_files = [n for n in zf.namelist() if n.endswith(".json") and "_" in n]
        result = json.loads(zf.read(result_files[0]))

    out = out_dir / f"blast_{rid}.json"
    out.write_text(json.dumps(result, indent=2))
    log(TAG, f"Saved -> {out.name}")
    return out


# ── Phobius ───────────────────────────────────────────────────────────────────

_PHOB = "https://www.ebi.ac.uk/Tools/services/rest/phobius"


def run_phobius(fasta_path: Path, out_dir: Path, name: str) -> Path:
    TAG = "Phobius"
    log(TAG, "Submitting...")

    r = requests.post(f"{_PHOB}/run", data={
        "email": EMAIL, "sequence": fasta_path.read_text(),
        "format": "short", "stype": "protein",
    })
    r.raise_for_status()
    job_id = r.text.strip()
    log(TAG, f"Job: {job_id}")

    while True:
        status = requests.get(f"{_PHOB}/status/{job_id}").text.strip()
        log(TAG, status)
        if status in ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
            break
        time.sleep(10)

    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")

    result = requests.get(f"{_PHOB}/result/{job_id}/out").text
    out = out_dir / f"phobius_{job_id}.txt"
    out.write_text(result)
    log(TAG, f"Saved -> {out.name}")
    return out


# ── HMMER hmmscan ─────────────────────────────────────────────────────────────

_HMMER = "https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan"


def run_hmmer(fasta_path: Path, out_dir: Path, name: str) -> Path:
    TAG = "HMMER"
    log(TAG, "Submitting...")

    r = requests.post(f"{_HMMER}/run", data={
        "email": EMAIL, "sequence": fasta_path.read_text(),
        "database": "pfam", "E": "1.0",
    })
    r.raise_for_status()
    job_id = r.text.strip()
    log(TAG, f"Job: {job_id}")

    while True:
        status = requests.get(f"{_HMMER}/status/{job_id}").text.strip()
        log(TAG, status)
        if status in ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
            break
        time.sleep(10)

    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")

    result = requests.get(f"{_HMMER}/result/{job_id}/out").text
    out = out_dir / f"hmmer_{job_id}.txt"
    out.write_text(result)
    log(TAG, f"Saved -> {out.name}")
    return out


# ── FoldSeek ──────────────────────────────────────────────────────────────────

_FS = "https://search.foldseek.com/api"
_FS_DBS = ["afdb50", "afdb-swissprot", "afdb-proteome", "BFVD"]
_FS_DONE = {"COMPLETE", "ERROR", "FAILED", "UNKNOWN"}


def run_foldseek(pdb_path: Path, out_dir: Path, name: str) -> Path:
    TAG = "FoldSeek"
    log(TAG, "Submitting...")

    pdb_content = pdb_path.read_text()
    if not pdb_content.endswith("\n"):
        pdb_content += "\n"

    data = [("q", pdb_content), ("mode", "3diaa"), ("email", EMAIL)]
    for db in _FS_DBS:
        data.append(("database[]", db))

    r = requests.post(f"{_FS}/ticket", data=data)
    r.raise_for_status()
    resp = r.json()

    if resp.get("status") in ("RATELIMIT", "MAINTENANCE"):
        raise RuntimeError(f"Server returned '{resp['status']}'")

    ticket = resp["id"]
    log(TAG, f"Job: {ticket}")

    while True:
        status = requests.get(f"{_FS}/ticket/{ticket}").json().get("status", "UNKNOWN")
        log(TAG, status)
        if status in _FS_DONE:
            break
        time.sleep(10)

    if status != "COMPLETE":
        raise RuntimeError(f"Job ended with status '{status}'")

    result = requests.get(f"{_FS}/result/download/{ticket}").content
    out = out_dir / f"foldseek_{ticket}.tar.gz"
    out.write_bytes(result)
    log(TAG, f"Saved -> {out.name}")
    return out


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the full dark proteome annotation pipeline for one protein. "
            "Provide --fasta, --pdb, or both. All applicable tools run in parallel."
        )
    )
    parser.add_argument("--fasta", type=Path,
                        help="FASTA file -> InterProScan, BLASTp, Phobius, HMMER")
    parser.add_argument("--pdb",   type=Path,
                        help="PDB file  -> FoldSeek")
    parser.add_argument("--name",
                        help="Protein name for the output folder "
                             "(default: stem of the FASTA or PDB filename)")
    parser.add_argument("--output-dir", type=Path, default=Path("data/output"),
                        help="Parent output directory (default: data/output)")
    args = parser.parse_args()

    if not args.fasta and not args.pdb:
        parser.error("Provide at least one of --fasta or --pdb.")
    if args.fasta and not args.fasta.is_file():
        parser.error(f"FASTA file not found: {args.fasta}")
    if args.pdb and not args.pdb.is_file():
        parser.error(f"PDB file not found: {args.pdb}")

    name = args.name or (args.fasta or args.pdb).stem
    out_dir = args.output_dir / name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nProtein : {name}")
    print(f"Output  : {out_dir}")
    if args.fasta:
        print(f"FASTA   : {args.fasta}")
    if args.pdb:
        print(f"PDB     : {args.pdb}")
    print()

    tasks: dict[str, tuple] = {}
    if args.fasta:
        tasks["InterProScan"] = (run_interproscan, args.fasta)
        tasks["BLASTp"]       = (run_blast,        args.fasta)
        tasks["Phobius"]      = (run_phobius,       args.fasta)
        tasks["HMMER"]        = (run_hmmer,         args.fasta)
    if args.pdb:
        tasks["FoldSeek"]     = (run_foldseek, args.pdb)

    results: dict[str, Path | None] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(tasks)) as pool:
        futures = {
            pool.submit(fn, path, out_dir, name): tool
            for tool, (fn, path) in tasks.items()
        }
        for future in concurrent.futures.as_completed(futures):
            tool = futures[future]
            try:
                results[tool] = future.result()
            except Exception as exc:
                log(tool, f"ERROR: {exc}")
                results[tool] = None

    passed = sum(1 for p in results.values() if p)
    print(f"\n{'─' * 52}")
    print(f"Pipeline complete  |  protein: {name}")
    print(f"Output folder: {out_dir}\n")
    for tool, path in results.items():
        if path:
            print(f"  [OK]     {tool:<15} {path.name}")
        else:
            print(f"  [FAILED] {tool}")
    print(f"\n{passed}/{len(tasks)} tools succeeded.")


if __name__ == "__main__":
    main()
