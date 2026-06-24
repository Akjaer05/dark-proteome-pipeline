"""
Dark Proteome Annotation Pipeline — Streamlit web interface.

Run:
    streamlit run dark_proteome_app.py
"""

import io
import json
import os
import tarfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
import streamlit as st
from dotenv import load_dotenv

load_dotenv()
EMAIL = os.environ.get("EBI_EMAIL", "")

# ── Page config ───────────────────────────────────────────────────────────────

st.set_page_config(page_title="Dark Proteome Pipeline", layout="wide")
st.title("Dark Proteome Annotation Pipeline")

if not EMAIL:
    st.error("EBI_EMAIL not set. Add `EBI_EMAIL=your@email.com` to .env and restart.")
    st.stop()

# ── Tool runners ──────────────────────────────────────────────────────────────

def _ebi_poll(base_url, job_id):
    terminal = {"FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"}
    while True:
        r = requests.get(f"{base_url}/status/{job_id}", timeout=30)
        r.raise_for_status()
        s = r.text.strip()
        if s in terminal:
            return s
        time.sleep(10)


def run_interproscan(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "goterms": "true", "pathways": "true", "stype": "p",
    }, timeout=30)
    r.raise_for_status()
    jid = r.text.strip()
    if _ebi_poll(url, jid) != "FINISHED":
        raise RuntimeError("Job did not finish successfully")
    r = requests.get(f"{url}/result/{jid}/json",
                     headers={"Accept": "application/json"}, timeout=60)
    r.raise_for_status()
    return r.json()


def run_blast(sequence):
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    r = requests.put(url, params={
        "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
        "QUERY": sequence, "FORMAT_TYPE": "JSON2",
        "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
    }, timeout=30)
    r.raise_for_status()
    rid, rtoe = None, 30
    for line in r.text.splitlines():
        if line.startswith("    RID = "):
            rid = line.split("=", 1)[1].strip()
        elif line.startswith("    RTOE = "):
            rtoe = int(line.split("=", 1)[1].strip())
    if not rid:
        raise RuntimeError("No RID in NCBI response")
    time.sleep(rtoe)
    while True:
        r = requests.get(url, params={
            "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo",
            "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
        }, timeout=30)
        status = "UNKNOWN"
        for line in r.text.splitlines():
            if "Status=" in line:
                status = line.strip().split("=", 1)[1].strip()
                break
        if status in ("READY", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)
    if status != "READY":
        raise RuntimeError(f"BLAST status: {status}")
    r = requests.get(url, params={
        "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2",
        "DESCRIPTIONS": 10, "ALIGNMENTS": 10,
        "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
    }, timeout=60)
    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        names = [n for n in zf.namelist() if n.endswith(".json")]
        return json.loads(zf.read(names[0]))


def run_phobius(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/phobius"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "format": "short", "stype": "protein",
    }, timeout=30)
    r.raise_for_status()
    jid = r.text.strip()
    if _ebi_poll(url, jid) != "FINISHED":
        raise RuntimeError("Job did not finish successfully")
    r = requests.get(f"{url}/result/{jid}/out", timeout=30)
    r.raise_for_status()
    return r.text


def run_hmmer(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "database": "pfam", "E": "1.0",
    }, timeout=30)
    r.raise_for_status()
    jid = r.text.strip()
    if _ebi_poll(url, jid) != "FINISHED":
        raise RuntimeError("Job did not finish successfully")
    r = requests.get(f"{url}/result/{jid}/out", timeout=30)
    r.raise_for_status()
    return r.text


def run_foldseek(pdb_text):
    url = "https://search.foldseek.com/api"
    if not pdb_text.endswith("\n"):
        pdb_text += "\n"
    data = [("q", pdb_text), ("mode", "3diaa"), ("email", EMAIL)]
    for db in ["afdb50", "afdb-swissprot", "afdb-proteome", "BFVD"]:
        data.append(("database[]", db))
    r = requests.post(f"{url}/ticket", data=data, timeout=60)
    r.raise_for_status()
    resp = r.json()
    if resp.get("status") in ("RATELIMIT", "MAINTENANCE"):
        raise RuntimeError(f"FoldSeek server: {resp['status']}")
    ticket = resp["id"]
    while True:
        s = requests.get(f"{url}/ticket/{ticket}", timeout=30).json().get("status", "UNKNOWN")
        if s in ("COMPLETE", "ERROR", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)
    if s != "COMPLETE":
        raise RuntimeError(f"FoldSeek status: {s}")
    r = requests.get(f"{url}/result/download/{ticket}", timeout=120)
    r.raise_for_status()
    return r.content  # raw tar.gz bytes


# ── UI ────────────────────────────────────────────────────────────────────────

fasta_file = st.file_uploader(
    "FASTA file — InterProScan · BLASTp · Phobius · HMMER",
    type=["fasta", "fa", "txt"],
)
pdb_file = st.file_uploader(
    "PDB file (optional) — FoldSeek structural search",
    type=["pdb"],
)

run = st.button(
    "Run Pipeline",
    type="primary",
    disabled=not (fasta_file or pdb_file),
)

if run:
    fasta_text = fasta_file.read().decode() if fasta_file else None
    pdb_text   = pdb_file.read().decode()   if pdb_file   else None

    tasks = {}
    if fasta_text:
        tasks["InterProScan"] = (run_interproscan, fasta_text)
        tasks["BLASTp"]       = (run_blast,        fasta_text)
        tasks["Phobius"]      = (run_phobius,       fasta_text)
        tasks["HMMER"]        = (run_hmmer,         fasta_text)
    if pdb_text:
        tasks["FoldSeek"]     = (run_foldseek,      pdb_text)

    results = {}
    with st.spinner(f"Running {len(tasks)} tools in parallel — typically 8–12 minutes…"):
        with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
            futures = {pool.submit(fn, arg): name for name, (fn, arg) in tasks.items()}
            for future in as_completed(futures):
                name = futures[future]
                try:
                    results[name] = {"ok": True, "data": future.result()}
                except Exception as exc:
                    results[name] = {"ok": False, "error": str(exc)}

    n_ok = sum(r["ok"] for r in results.values())
    if n_ok == len(results):
        st.success(f"All {n_ok} tools completed.")
    else:
        st.warning(f"{n_ok}/{len(results)} tools completed.")

    for name in tasks:  # show in submission order
        if name not in results:
            continue
        r = results[name]
        with st.expander(f"{'✓' if r['ok'] else '✗'}  {name}", expanded=r["ok"]):
            if not r["ok"]:
                st.error(r["error"])
                continue

            data = r["data"]

            if name in ("InterProScan", "BLASTp"):
                st.json(data)

            elif name == "FoldSeek":
                tar_bytes = data
                with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tar:
                    m8s = [m for m in tar.getmembers() if m.name.endswith(".m8")]
                    if not m8s:
                        st.info("No result files found.")
                    else:
                        db_tabs = st.tabs([
                            m.name.split("alis_")[-1].replace(".m8", "") for m in m8s
                        ])
                        for tab, member in zip(db_tabs, m8s):
                            with tab:
                                f = tar.extractfile(member)
                                content = f.read().decode("utf-8", errors="replace") if f else ""
                                hit_count = len([l for l in content.splitlines() if l.strip()])
                                st.caption(f"{hit_count} hits")
                                st.text(content)

            else:  # Phobius, HMMER — plain text
                st.text(data)
