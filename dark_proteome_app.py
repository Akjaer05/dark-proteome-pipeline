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

import pandas as pd
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


# ── Result display ────────────────────────────────────────────────────────────

def show_interproscan(data: dict) -> None:
    rows = []
    for res in data.get("results", []):
        for match in res.get("matches", []):
            sig   = match.get("signature", {})
            entry = sig.get("entry") or {}
            lib   = sig.get("signatureLibraryRelease", {}).get("library", "")
            desc  = entry.get("description") or sig.get("description") or ""
            for loc in match.get("locations", []):
                rows.append({
                    "Database":    lib,
                    "Accession":   sig.get("accession", ""),
                    "Name":        sig.get("name", ""),
                    "Description": desc,
                    "Start":       loc.get("start", ""),
                    "End":         loc.get("end", ""),
                    "E-value":     match.get("evalue", ""),
                })
    if not rows:
        st.info("No domain matches found.")
        return
    df = pd.DataFrame(rows)
    c1, c2 = st.columns(2)
    c1.metric("Matches", len(df))
    c2.metric("Databases hit", df["Database"].nunique())
    st.dataframe(df, use_container_width=True, hide_index=True)


def show_blast(data: dict) -> None:
    try:
        bo2  = data["BlastOutput2"]
        bo2  = bo2[0] if isinstance(bo2, list) else bo2
        hits = bo2["report"]["results"]["search"]["hits"]
    except (KeyError, IndexError, TypeError):
        st.info("No BLAST hits found.")
        return
    if not hits:
        st.info("No BLAST hits found.")
        return
    rows = []
    for hit in hits:
        desc      = (hit.get("description") or [{}])[0]
        hsp       = (hit.get("hsps")        or [{}])[0]
        align_len = hsp.get("align_len") or 1
        rows.append({
            "Accession":   desc.get("accession", ""),
            "Description": (desc.get("title") or "")[:80],
            "Organism":    desc.get("sciname", ""),
            "% Identity":  round((hsp.get("identity") or 0) / align_len * 100, 1),
            "E-value":     hsp.get("evalue", ""),
            "Bit score":   hsp.get("bit_score", ""),
        })
    df = pd.DataFrame(rows)
    st.metric("Hits returned", len(df))
    st.dataframe(df, use_container_width=True, hide_index=True)


def show_phobius(text: str) -> None:
    for line in text.splitlines():
        if not line.strip() or line.startswith("SEQENCE"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            c1, c2, c3 = st.columns(3)
            c1.metric("TM helices",    parts[1])
            c2.metric("Signal peptide", "Yes" if parts[2] == "1" else "No")
            c3.metric("Topology",       parts[3])
            st.caption("Topology codes: o = outside, i = inside, h = TM helix")
            return
    st.info("No Phobius results parsed.")


def show_hmmer(text: str) -> None:
    rows, in_hits, below = [], False, False
    for line in text.splitlines():
        if "Scores for complete sequence" in line:
            in_hits = True;  continue
        if not in_hits:
            continue
        if "Domain annotation" in line:
            break
        s = line.strip()
        if not s or s.startswith("E-value"):
            continue
        # Must check inclusion threshold before the "---" guard — the line
        # "------ inclusion threshold ------" starts with "---" and would
        # otherwise be skipped, leaving `below` permanently False.
        if "inclusion threshold" in s:
            below = True;  continue
        if s.startswith("---"):
            continue
        parts = s.split(None, 9)
        if len(parts) < 9:
            continue
        try:
            float(parts[0])
        except ValueError:
            continue
        rows.append({
            "Model (Pfam)": parts[8],
            "Description":  parts[9] if len(parts) > 9 else "",
            "E-value":      parts[0],
            "Score":        float(parts[1]),
            "N domains":    int(parts[7]),
            "Significant":  "Yes" if not below else "No",
        })
    if not rows:
        st.info("No Pfam domain matches found.")
        return
    df = pd.DataFrame(rows)
    above = int((df["Significant"] == "Yes").sum())
    c1, c2 = st.columns(2)
    c1.metric("Significant domains", above)
    c2.metric("Total (incl. below threshold)", len(df))
    st.dataframe(df, use_container_width=True, hide_index=True)


def show_foldseek(tar_bytes: bytes) -> None:
    cols = ["query","target","fident","alnlen","mismatch",
            "gapopen","qstart","qend","tstart","tend","prob","evalue","bits"]
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tar:
        m8s = [m for m in tar.getmembers()
               if m.name.endswith(".m8") and "_report" not in m.name]
        if not m8s:
            st.info("No result files in archive.")
            return
        db_names = [m.name.split("alis_")[-1].replace(".m8", "") for m in m8s]
        total = 0
        frames = {}
        for member, db in zip(m8s, db_names):
            f = tar.extractfile(member)
            if not f:
                continue
            lines = [l for l in f.read().decode("utf-8", errors="replace").splitlines() if l.strip()]
            total += len(lines)
            if not lines:
                frames[db] = pd.DataFrame()
                continue
            rows = []
            for line in lines:
                p = line.split("\t")
                rows.append({c: (p[i] if i < len(p) else "") for i, c in enumerate(cols)})
            df = pd.DataFrame(rows)
            for c in ("fident", "prob", "evalue", "bits"):
                df[c] = pd.to_numeric(df[c], errors="coerce")
            df["target"] = df["target"].str[:60]
            df = df[["target","fident","prob","evalue","bits","qstart","qend"]].copy()
            df.columns = ["Target","Identity %","Probability","E-value","Bits","Q.Start","Q.End"]
            frames[db] = df

        st.metric("Total hits across all databases", total)
        db_tabs = st.tabs(db_names)
        for tab, db in zip(db_tabs, db_names):
            with tab:
                df = frames.get(db, pd.DataFrame())
                st.caption(f"{len(df)} hits")
                if not df.empty:
                    st.dataframe(df, use_container_width=True, hide_index=True)
                else:
                    st.info("No hits.")


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
            if   name == "InterProScan": show_interproscan(data)
            elif name == "BLASTp":       show_blast(data)
            elif name == "Phobius":      show_phobius(data)
            elif name == "HMMER":        show_hmmer(data)
            elif name == "FoldSeek":     show_foldseek(data)
