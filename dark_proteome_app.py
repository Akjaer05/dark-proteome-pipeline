"""
Streamlit web interface for the Dark Proteome Annotation Pipeline.

Launch with:
    streamlit run dark_proteome_app.py
"""

import io
import json
import os
import tarfile
import threading
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
import streamlit as st
from dotenv import load_dotenv

load_dotenv()

EMAIL = os.environ.get("EBI_EMAIL", "")
EBI_POLL_INTERVAL = 10  # seconds

# ── Module-level shared job state ─────────────────────────────────────────────
# Each session writes its job results here; keyed by UUID so sessions don't clash.

_JOBS: dict = {}
_JOBS_LOCK = threading.Lock()


def _job_set(job_id: str, key: str, value) -> None:
    with _JOBS_LOCK:
        _JOBS[job_id][key] = value


def _tool_update(job_id: str, tool: str, **kwargs) -> None:
    with _JOBS_LOCK:
        _JOBS[job_id]["tools"][tool].update(kwargs)


# ── EBI REST polling helper ───────────────────────────────────────────────────

def _ebi_poll(base_url: str, job_id: str) -> str:
    terminal = {"FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"}
    while True:
        r = requests.get(f"{base_url}/status/{job_id}", timeout=30)
        r.raise_for_status()
        status = r.text.strip()
        if status in terminal:
            return status
        time.sleep(EBI_POLL_INTERVAL)


# ── Tool runners (raise RuntimeError on failure, no sys.exit) ─────────────────

def _run_interproscan(sequence: str) -> dict:
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "goterms": "true", "pathways": "true", "stype": "p",
    }, timeout=30)
    r.raise_for_status()
    job_id = r.text.strip()
    status = _ebi_poll(url, job_id)
    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")
    r = requests.get(f"{url}/result/{job_id}/json",
                     headers={"Accept": "application/json"}, timeout=60)
    r.raise_for_status()
    return r.json()


def _run_blast(sequence: str) -> dict:
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
        raise RuntimeError("Could not parse RID from NCBI response")
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
        raise RuntimeError(f"BLAST ended with status '{status}'")
    r = requests.get(url, params={
        "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2",
        "DESCRIPTIONS": 10, "ALIGNMENTS": 10,
        "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
    }, timeout=60)
    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        names = [n for n in zf.namelist() if n.endswith(".json")]
        return json.loads(zf.read(names[0]))


def _run_phobius(sequence: str) -> str:
    url = "https://www.ebi.ac.uk/Tools/services/rest/phobius"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "format": "short", "stype": "protein",
    }, timeout=30)
    r.raise_for_status()
    job_id = r.text.strip()
    status = _ebi_poll(url, job_id)
    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")
    r = requests.get(f"{url}/result/{job_id}/out", timeout=30)
    r.raise_for_status()
    return r.text


def _run_hmmer(sequence: str) -> str:
    url = "https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan"
    r = requests.post(f"{url}/run", data={
        "email": EMAIL, "sequence": sequence,
        "database": "pfam", "E": "1.0",
    }, timeout=30)
    r.raise_for_status()
    job_id = r.text.strip()
    status = _ebi_poll(url, job_id)
    if status != "FINISHED":
        raise RuntimeError(f"Job ended with status '{status}'")
    r = requests.get(f"{url}/result/{job_id}/out", timeout=30)
    r.raise_for_status()
    return r.text


def _run_foldseek(pdb_text: str) -> bytes:
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
        raise RuntimeError(f"Server returned '{resp['status']}'")
    ticket = resp["id"]
    while True:
        status = requests.get(f"{url}/ticket/{ticket}", timeout=30).json().get("status", "UNKNOWN")
        if status in ("COMPLETE", "ERROR", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)
    if status != "COMPLETE":
        raise RuntimeError(f"Job ended with status '{status}'")
    r = requests.get(f"{url}/result/download/{ticket}", timeout=120)
    r.raise_for_status()
    return r.content


# ── Result parsers ────────────────────────────────────────────────────────────

def _parse_interproscan(data: dict) -> pd.DataFrame:
    rows = []
    for res in data.get("results", []):
        for match in res.get("matches", []):
            sig = match.get("signature", {})
            lib = sig.get("signatureLibraryRelease", {}).get("library", "")
            entry = sig.get("entry") or {}
            go = "; ".join(g.get("id", "") for g in entry.get("goXRefs", []))
            for loc in match.get("locations", []):
                rows.append({
                    "Database": lib,
                    "Accession": sig.get("accession", ""),
                    "Name": sig.get("name", ""),
                    "IPR Entry": entry.get("accession", ""),
                    "Description": entry.get("description") or sig.get("description") or "",
                    "Start": loc.get("start", ""),
                    "End": loc.get("end", ""),
                    "E-value": match.get("evalue", ""),
                    "GO Terms": go,
                })
    return pd.DataFrame(rows)


def _parse_blast(data: dict) -> pd.DataFrame:
    try:
        bo2 = data["BlastOutput2"]
        if isinstance(bo2, list):
            bo2 = bo2[0]
        hits = bo2["report"]["results"]["search"]["hits"]
    except (KeyError, IndexError, TypeError):
        return pd.DataFrame()
    rows = []
    for hit in hits:
        desc = (hit.get("description") or [{}])[0]
        hsp = (hit.get("hsps") or [{}])[0]
        align_len = hsp.get("align_len") or 1
        rows.append({
            "Accession": desc.get("accession", ""),
            "Title": (desc.get("title") or "")[:100],
            "Species": desc.get("sciname", ""),
            "E-value": hsp.get("evalue", ""),
            "Bit score": hsp.get("bit_score", ""),
            "% Identity": round((hsp.get("identity") or 0) / align_len * 100, 1),
            "Align len": align_len,
            "Q.start": hsp.get("query_from", ""),
            "Q.end": hsp.get("query_to", ""),
        })
    return pd.DataFrame(rows)


def _parse_phobius(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("SEQENCE"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            rows.append({
                "Protein ID": parts[0],
                "TM helices": int(parts[1]),
                "Signal peptide": "Yes" if parts[2] == "1" else "No",
                "Topology": parts[3],
            })
    return pd.DataFrame(rows)


def _parse_hmmer(text: str) -> pd.DataFrame:
    rows = []
    in_hits = False
    below_threshold = False
    for line in text.splitlines():
        if "Scores for complete sequence" in line:
            in_hits = True
            continue
        if not in_hits:
            continue
        if "Domain annotation" in line:
            break
        stripped = line.strip()
        if not stripped or stripped.startswith("E-value") or stripped.startswith("---"):
            continue
        if "inclusion threshold" in stripped:
            below_threshold = True
            continue
        # split(None, 9) → at most 10 parts; index 8 = model, 9 = rest of description
        parts = stripped.split(None, 9)
        if len(parts) < 9:
            continue
        try:
            float(parts[0])
        except ValueError:
            continue
        rows.append({
            "Model (Pfam)": parts[8],
            "Description": parts[9] if len(parts) > 9 else "",
            "E-value": parts[0],
            "Score": float(parts[1]),
            "Bias": float(parts[2]),
            "N domains": int(parts[7]),
            "Included": "Yes" if not below_threshold else "No (below threshold)",
        })
    return pd.DataFrame(rows)


def _parse_foldseek(tar_bytes: bytes) -> dict:
    base_cols = [
        "query", "target", "fident", "alnlen", "mismatch",
        "gapopen", "qstart", "qend", "tstart", "tend",
        "prob", "evalue", "bits",
    ]
    display_cols = [
        "Target", "Identity %", "Probability", "E-value", "Bits",
        "Q.Start", "Q.End", "T.Start", "T.End",
    ]
    result = {}
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".m8"):
                continue
            stem = Path(member.name).stem  # e.g. "alis_afdb50"
            db = stem[5:] if stem.startswith("alis_") else stem
            f = tar.extractfile(member)
            if not f:
                continue
            content = f.read().decode("utf-8", errors="replace")
            if not content.strip():
                result[db] = pd.DataFrame(columns=display_cols)
                continue
            rows = []
            for line in content.splitlines():
                parts = line.split("\t")
                rows.append({c: (parts[i] if i < len(parts) else "") for i, c in enumerate(base_cols)})
            df = pd.DataFrame(rows)
            for c in ("fident", "prob", "evalue", "bits"):
                if c in df.columns:
                    df[c] = pd.to_numeric(df[c], errors="coerce")
            df = df[["target", "fident", "prob", "evalue", "bits",
                      "qstart", "qend", "tstart", "tend"]].copy()
            df.columns = display_cols
            result[db] = df
    return result


# ── Result renderers ──────────────────────────────────────────────────────────

def _render_interproscan(result: dict) -> None:
    df = _parse_interproscan(result)
    if df.empty:
        st.info("No domain matches found.")
    else:
        c1, c2 = st.columns(2)
        c1.metric("Total matches", len(df))
        c2.metric("Databases with hits", df["Database"].nunique())
        st.dataframe(df, use_container_width=True, hide_index=True)
    with st.expander("Raw JSON"):
        st.json(result)


def _render_blast(result: dict) -> None:
    df = _parse_blast(result)
    if df.empty:
        st.info("No BLAST hits found.")
    else:
        st.metric("Top hits returned", len(df))
        st.dataframe(df, use_container_width=True, hide_index=True)
    with st.expander("Raw JSON"):
        st.json(result)


def _render_phobius(result: str) -> None:
    df = _parse_phobius(result)
    if df.empty:
        st.info("No Phobius results.")
    else:
        for _, row in df.iterrows():
            c1, c2, c3 = st.columns(3)
            c1.metric("TM helices", row["TM helices"])
            c2.metric("Signal peptide", row["Signal peptide"])
            c3.metric("Topology code", row["Topology"])
        st.caption("Topology: o = non-cytoplasmic side, i = cytoplasmic, h = TM helix region")
    with st.expander("Raw output"):
        st.text(result)


def _render_hmmer(result: str) -> None:
    df = _parse_hmmer(result)
    if df.empty:
        st.info("No Pfam domain matches found.")
    else:
        above = int((df["Included"] == "Yes").sum())
        c1, c2 = st.columns(2)
        c1.metric("Domains above threshold", above)
        c2.metric("Total matches (incl. below threshold)", len(df))
        st.dataframe(df, use_container_width=True, hide_index=True)
    with st.expander("Raw output"):
        st.text(result)


def _render_foldseek(result: bytes) -> None:
    dbs = _parse_foldseek(result)
    if not dbs:
        st.info("No FoldSeek results.")
        return
    total = sum(len(df) for df in dbs.values())
    st.metric("Total structural hits across all databases", total)
    sub_tabs = st.tabs(list(dbs.keys()))
    for tab, (db_name, df) in zip(sub_tabs, dbs.items()):
        with tab:
            st.caption(f"{len(df)} hits in **{db_name}**")
            if not df.empty:
                st.dataframe(df, use_container_width=True, hide_index=True)
            else:
                st.info("No hits in this database.")


_RENDERERS = {
    "InterProScan": _render_interproscan,
    "BLASTp": _render_blast,
    "Phobius": _render_phobius,
    "HMMER": _render_hmmer,
    "FoldSeek": _render_foldseek,
}

_STATUS_LABELS = {
    "QUEUED": "Queued",
    "RUNNING": "Running...",
    "FINISHED": "Done",
    "FAILED": "Failed",
}


# ── Background job runner ─────────────────────────────────────────────────────

def _job_runner(job_id: str, fasta_text, pdb_text) -> None:
    tasks = {}
    if fasta_text:
        tasks["InterProScan"] = (_run_interproscan, fasta_text)
        tasks["BLASTp"]       = (_run_blast,        fasta_text)
        tasks["Phobius"]      = (_run_phobius,       fasta_text)
        tasks["HMMER"]        = (_run_hmmer,         fasta_text)
    if pdb_text:
        tasks["FoldSeek"]     = (_run_foldseek,      pdb_text)

    with _JOBS_LOCK:
        _JOBS[job_id]["tools"] = {
            name: {"status": "QUEUED", "result": None, "error": None}
            for name in tasks
        }

    def run_one(name, fn, arg):
        _tool_update(job_id, name, status="RUNNING")
        try:
            _tool_update(job_id, name, status="FINISHED", result=fn(arg))
        except Exception as exc:
            _tool_update(job_id, name, status="FAILED", error=str(exc))

    with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
        futures = {pool.submit(run_one, n, fn, arg): n for n, (fn, arg) in tasks.items()}
        for _ in as_completed(futures):
            pass

    _job_set(job_id, "done", True)


# ── Streamlit page config (must be first Streamlit call) ─────────────────────

st.set_page_config(
    page_title="Dark Proteome Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)


# ── Main app ──────────────────────────────────────────────────────────────────

def main():
    st.title("Dark Proteome Annotation Pipeline")
    st.caption(
        "Upload a protein sequence (FASTA) and optionally its predicted structure (PDB) "
        "to run five annotation tools in parallel: InterProScan, BLASTp, Phobius, HMMER, and FoldSeek."
    )

    if not EMAIL:
        st.error(
            "EBI_EMAIL is not set. Create a `.env` file in the project root with "
            "`EBI_EMAIL=your@email.com`, then restart the app."
        )
        return

    # Silently drop any job ID left over from a previous server session.
    # Do this before rendering anything so the form and job views never clash.
    stored_id = st.session_state.get("job_id")
    if stored_id and stored_id not in _JOBS:
        st.session_state.pop("job_id", None)
        st.session_state.pop("job_name", None)
        stored_id = None

    # ── Upload widgets (shown only when no job is active) ────────────────────
    # NOTE: st.file_uploader must NOT be inside st.form — Streamlit clears
    # uploader widget state on form submission, making the file appear as None
    # in the submit handler. Plain widgets + st.button avoid this entirely.
    if not stored_id:
        col_fasta, col_pdb = st.columns(2)
        with col_fasta:
            fasta_upload = st.file_uploader(
                "FASTA file — InterProScan · BLASTp · Phobius · HMMER",
                type=["fasta", "fa", "txt"],
            )
        with col_pdb:
            pdb_upload = st.file_uploader(
                "PDB file (optional) — FoldSeek structural search",
                type=["pdb"],
            )
        protein_name = st.text_input(
            "Protein name",
            placeholder="e.g. SLINKY  (defaults to filename stem if left blank)",
        )

        if st.button("Run Pipeline", type="primary", use_container_width=True):
            if not fasta_upload and not pdb_upload:
                st.error("Please upload at least one file.")
            else:
                fasta_text = fasta_upload.read().decode("utf-8") if fasta_upload else None
                pdb_text = pdb_upload.read().decode("utf-8") if pdb_upload else None
                name = protein_name.strip() or (
                    Path(fasta_upload.name).stem if fasta_upload else Path(pdb_upload.name).stem
                )

                job_id = str(uuid.uuid4())
                with _JOBS_LOCK:
                    _JOBS[job_id] = {"done": False, "tools": {}, "name": name}

                threading.Thread(
                    target=_job_runner, args=(job_id, fasta_text, pdb_text), daemon=True
                ).start()

                st.session_state["job_id"] = job_id
                st.session_state["job_name"] = name
                st.rerun()

        return  # nothing more to show until a job is running

    # ── Job status and results ────────────────────────────────────────────────
    job_id = stored_id
    job = _JOBS[job_id]  # guaranteed to exist by the stale-check above

    with _JOBS_LOCK:
        done = job["done"]
        name = job.get("name", "")
        tools = {k: dict(v) for k, v in job["tools"].items()}

    # Header
    hdr_col, btn_col = st.columns([6, 1])
    with hdr_col:
        st.subheader(f"{'Completed' if done else 'Running'}: {name}")
    with btn_col:
        if st.button("New search", use_container_width=True):
            st.session_state.pop("job_id", None)
            st.session_state.pop("job_name", None)
            st.rerun()

    # Per-tool status cards
    if tools:
        cols = st.columns(len(tools))
        for col, (tool_name, state) in zip(cols, tools.items()):
            col.metric(
                label=tool_name,
                value=_STATUS_LABELS.get(state["status"], state["status"]),
            )

    # Inline error messages for failed tools
    for tool_name, state in tools.items():
        if state["status"] == "FAILED":
            st.error(f"**{tool_name}** failed: {state.get('error', 'Unknown error')}")

    # Results tabs — one per finished tool, appearing as they complete
    finished = {n: s for n, s in tools.items() if s["status"] == "FINISHED"}
    if finished:
        st.divider()
        tabs = st.tabs(list(finished.keys()))
        for tab, (tool_name, state) in zip(tabs, finished.items()):
            with tab:
                renderer = _RENDERERS.get(tool_name)
                if renderer:
                    renderer(state["result"])

    # Keep refreshing until the job is done
    if not done:
        time.sleep(3)
        st.rerun()


main()
