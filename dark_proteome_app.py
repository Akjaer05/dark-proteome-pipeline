"""
DarkProteome — Professional scientific annotation tool for the microbial dark proteome.
Run:  streamlit run dark_proteome_app.py
"""

import base64
import hashlib
import html as _html
import io
import json
import math
import os
import pathlib
import pickle
import tarfile
import tempfile
import time
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
try:
    from Bio.PDB.DSSP import DSSP as _DSSP
except Exception:
    _DSSP = None
from dotenv import load_dotenv

load_dotenv()
EMAIL  = os.environ.get("EBI_EMAIL", "")
GITHUB = "https://github.com/Akjaer05/dark-proteome-pipeline"

# ── Result persistence helpers ─────────────────────────────────────────────────

_CACHE_DIR = pathlib.Path(tempfile.gettempdir())


def _cache_key(fasta_text, pdb_text) -> str:
    payload = ((fasta_text or "") + "|||" + (pdb_text or "")).encode()
    return hashlib.md5(payload).hexdigest()[:14]


def _cache_path(fasta_text, pdb_text) -> pathlib.Path:
    return _CACHE_DIR / f"dpp_{_cache_key(fasta_text, pdb_text)}.pkl"


def _save_cache(fasta_text, pdb_text, results: dict, active_tools: list) -> None:
    try:
        payload = {
            "results":      results,
            "fasta_text":   fasta_text,
            "pdb_text":     pdb_text,
            "active_tools": active_tools,
        }
        _cache_path(fasta_text, pdb_text).write_bytes(pickle.dumps(payload))
    except Exception:
        pass  # cache failure is never fatal


def _load_cache(fasta_text, pdb_text):
    try:
        p = _cache_path(fasta_text, pdb_text)
        if p.exists():
            return pickle.loads(p.read_bytes())
    except Exception:
        pass
    return None


def _extract_protein_name(fasta_text: str) -> str:
    if not fasta_text:
        return "protein"
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            first = line[1:].strip().split()[0] if line[1:].strip() else ""
            return first or "protein"
    return "protein"


def _results_to_json_bytes(results: dict, protein_name: str = "protein") -> bytes:
    """Serialise pipeline results to JSON bytes for download.
    FoldSeek tar.gz bytes are base64-encoded; everything else is JSON-native."""
    out: dict = {}
    for tool, r in results.items():
        entry: dict = {"ok": r["ok"]}
        if not r["ok"]:
            entry["error"] = r.get("error", "")
        else:
            d = r.get("data")
            if isinstance(d, bytes):
                entry["data"]     = base64.b64encode(d).decode()
                entry["encoding"] = "base64"
            else:
                try:
                    json.dumps(d)
                    entry["data"] = d
                except (TypeError, ValueError):
                    entry["data"]      = str(d)[:2000]
                    entry["truncated"] = True
        out[tool] = entry
    return json.dumps({
        "schema":      "dark-proteome-pipeline/v1",
        "protein":     protein_name,
        "exported_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results":     out,
    }, indent=2).encode()


# ── Page config ────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="DarkProteome",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Global CSS ─────────────────────────────────────────────────────────────────

st.markdown("""
<style>
/* ── Hide Streamlit chrome ─────────────────────────────── */
#MainMenu, footer { visibility: hidden; }
[data-testid="stDecoration"],
[data-testid="stStatusWidget"]  { display: none !important; }
[data-testid="stHeader"]        { background: transparent !important; height: 0 !important; }
[data-testid="stToolbar"]       { display: none !important; }

/* ── Base ──────────────────────────────────────────────── */
html, body, .stApp,
[data-testid="stAppViewContainer"],
[data-testid="stMain"], .main {
    background-color: #0a0e1a !important;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Inter', sans-serif !important;
    color: #94a3b8 !important;
}
.block-container { padding: 0 !important; max-width: 100% !important; }

/* ── Column layout — sidebar only on 2-column layouts ──────────────────────── */
[data-testid="stHorizontalBlock"] {
    gap: 0 !important;
    align-items: stretch !important;
}
[data-testid="stHorizontalBlock"]:has(> [data-testid="column"]:last-child:nth-child(2)) > [data-testid="column"]:first-child {
    background: #080c17 !important;
    border-right: 1px solid #1e2d4a !important;
    padding: 28px 20px 60px !important;
    min-height: 65vh;
}
[data-testid="stHorizontalBlock"]:has(> [data-testid="column"]:last-child:nth-child(2)) > [data-testid="column"]:last-child {
    padding: 22px 28px 60px !important;
}

/* ── Widget labels ──────────────────────────────────────── */
[data-testid="stWidgetLabel"] p,
[data-testid="stWidgetLabel"] {
    color: #475569 !important;
    font-size: 10px !important;
    font-weight: 700 !important;
    letter-spacing: 0.09em !important;
    text-transform: uppercase !important;
}

/* ── File uploader ──────────────────────────────────────── */
[data-testid="stFileUploaderDropzone"] {
    background: #0d1424 !important;
    border: 1px dashed #1e2d4a !important;
    border-radius: 8px !important;
    transition: border-color 0.15s, background 0.15s !important;
}
[data-testid="stFileUploaderDropzone"]:hover {
    border-color: #3b82f6 !important;
    background: rgba(59, 130, 246, 0.03) !important;
}
[data-testid="stFileUploaderDropzone"] p,
[data-testid="stFileUploaderDropzone"] span,
[data-testid="stFileUploaderDropzone"] small {
    color: #334155 !important;
    font-size: 12px !important;
}
[data-testid="stFileUploader"] { margin-bottom: 4px !important; }

/* ── Buttons ────────────────────────────────────────────── */
.stButton > button {
    background: linear-gradient(135deg, #3b82f6, #2563eb) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: 8px !important;
    font-weight: 600 !important;
    font-size: 13px !important;
    letter-spacing: 0.02em !important;
    padding: 10px 20px !important;
    width: 100% !important;
    box-shadow: 0 4px 14px rgba(59, 130, 246, 0.25) !important;
    transition: all 0.2s !important;
    margin-top: 6px !important;
}
.stButton > button:hover:not(:disabled) {
    box-shadow: 0 8px 28px rgba(59, 130, 246, 0.55) !important;
    transform: translateY(-2px) scale(1.025) !important;
}
@keyframes dpp-btn-glow {
    0%,100% { box-shadow: 0 4px 14px rgba(59,130,246,.25); }
    50%      { box-shadow: 0 4px 26px rgba(59,130,246,.52); }
}
.stButton > button[kind="primary"]:not(:disabled) {
    animation: dpp-btn-glow 3s ease-in-out infinite !important;
}
.stButton > button:disabled {
    background: #0d1424 !important;
    color: #334155 !important;
    box-shadow: none !important;
    border: 1px solid #1e2d4a !important;
}

/* ── Tabs ───────────────────────────────────────────────── */
.stTabs [data-baseweb="tab-list"] {
    background: #080c17 !important;
    border-bottom: 1px solid #1e2d4a !important;
    padding: 0 4px !important;
    gap: 0 !important;
}
.stTabs [data-baseweb="tab"] {
    color: #475569 !important;
    background: transparent !important;
    border: none !important;
    border-bottom: 2px solid transparent !important;
    padding: 10px 16px !important;
    font-size: 12px !important;
    font-weight: 500 !important;
    transition: color 0.15s !important;
}
.stTabs [data-baseweb="tab"]:hover { color: #94a3b8 !important; }
.stTabs [aria-selected="true"] {
    color: #3b82f6 !important;
    border-bottom: 2px solid #3b82f6 !important;
    background: rgba(59, 130, 246, 0.04) !important;
}
.stTabs [data-baseweb="tab-panel"] {
    background: transparent !important;
    padding: 14px 0 0 !important;
    animation: dpp-panel-reveal 0.42s cubic-bezier(0.25,0.46,0.45,0.94) both !important;
}
@keyframes dpp-panel-reveal {
    from { opacity:0; transform:translateY(14px); }
    to   { opacity:1; transform:translateY(0); }
}

/* ── Alerts ─────────────────────────────────────────────── */
[data-testid="stAlert"] {
    background: #080c17 !important;
    border: 1px solid #1e2d4a !important;
    border-radius: 8px !important;
    color: #94a3b8 !important;
}

/* ── Spinner ────────────────────────────────────────────── */
[data-testid="stSpinner"] > div > div {
    border-top-color: #3b82f6 !important;
}

/* ── Scrollbar ──────────────────────────────────────────── */
::-webkit-scrollbar { width: 5px; height: 5px; }
::-webkit-scrollbar-track { background: #0a0e1a; }
::-webkit-scrollbar-thumb { background: #1e2d4a; border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: #3b82f6; }

/* ── Dark result table ──────────────────────────────────── */
.dp-wrap {
    overflow-x: auto;
    border-radius: 8px;
    border: 1px solid #1e2d4a;
    margin: 10px 0 4px;
}
.dp-tbl {
    width: 100%;
    border-collapse: collapse;
    font-family: 'JetBrains Mono', 'SF Mono', 'Consolas', monospace;
    font-size: 11.5px;
}
.dp-tbl thead th {
    background: #080c17;
    color: #3b82f6;
    font-size: 9.5px;
    font-weight: 700;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    padding: 10px 14px;
    border-bottom: 1px solid #1e2d4a;
    white-space: nowrap;
    text-align: left;
}
.dp-tbl tbody td {
    color: #94a3b8;
    padding: 8px 14px;
    border-bottom: 1px solid #0d1424;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 320px;
}
.dp-tbl tbody tr:last-child td { border-bottom: none; }
.dp-tbl tbody tr:hover td {
    background: rgba(59, 130, 246, 0.04);
    color: #c8d9f0;
}

/* ── Metric cards ───────────────────────────────────────── */
.dp-cards { display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 14px; }
.dp-card {
    background: #080c17;
    border: 1px solid #1e2d4a;
    border-radius: 8px;
    padding: 12px 18px;
    min-width: 110px;
}
.dp-card-label {
    color: #475569;
    font-size: 9.5px;
    font-weight: 700;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    margin-bottom: 4px;
}
.dp-card-value { color: #f0f6ff; font-size: 21px; font-weight: 700; line-height: 1.1; }
.dp-card-sub   { color: #334155; font-size: 10.5px; margin-top: 3px; }

/* ── Tool status list ───────────────────────────────────── */
.dp-status-row {
    display: flex;
    align-items: center;
    gap: 9px;
    padding: 7px 0;
    border-bottom: 1px solid #0d1424;
}
.dp-status-row:last-child { border-bottom: none; }
.dp-dot { width: 7px; height: 7px; border-radius: 50%; flex-shrink: 0; }
.dp-tool-name          { color: #475569; font-size: 12px; font-weight: 500; flex: 1; }
.dp-tool-name.active   { color: #94a3b8; }
.dp-tool-name.err      { color: #ef4444; }
.dp-badge { font-size: 10.5px; font-weight: 600; }
</style>
""", unsafe_allow_html=True)

# ── Navbar ─────────────────────────────────────────────────────────────────────

st.markdown(f"""
<nav style="background:#080c17; border-bottom:1px solid #1e2d4a; padding:0 28px;
            height:54px; display:flex; align-items:center; justify-content:space-between;
            position:sticky; top:0; z-index:200;">
  <div style="display:flex; align-items:center; gap:11px;">
    <div style="background:linear-gradient(135deg,#3b82f6,#1d4ed8); width:32px; height:32px;
                border-radius:8px; display:flex; align-items:center; justify-content:center;
                box-shadow:0 3px 10px rgba(59,130,246,0.35);">
      <svg width="16" height="16" viewBox="0 0 24 24" fill="none"
           stroke="#fff" stroke-width="2.3" stroke-linecap="round" stroke-linejoin="round">
        <path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8
                 a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/>
        <polyline points="3.27 6.96 12 12.01 20.73 6.96"/>
        <line x1="12" y1="22.08" x2="12" y2="12"/>
      </svg>
    </div>
    <span style="font-weight:700; font-size:16px; color:#f0f6ff;
                 letter-spacing:-0.01em;">DarkProteome</span>
  </div>
  <a href="{GITHUB}" target="_blank"
     style="display:flex; align-items:center; gap:6px; color:#475569;
            text-decoration:none; font-size:12px; font-weight:500;"
     onmouseover="this.style.color='#94a3b8'" onmouseout="this.style.color='#475569'">
    <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
      <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387
               .599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416
               -.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729
               1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997
               .107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931
               0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0
               1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404
               1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23
               .653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221
               0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293
               c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386
               0-6.627-5.373-12-12-12z"/>
    </svg>
    GitHub
  </a>
</nav>
""", unsafe_allow_html=True)

# ── Tool runners ──────────────────────────────────────────────────────────────

def _ebi_poll(base_url, job_id):
    terminal = {"FINISHED", "FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"}
    while True:
        r = requests.get(f"{base_url}/status/{job_id}", timeout=60)
        r.raise_for_status()
        s = r.text.strip()
        if s in terminal:
            return s
        time.sleep(10)


def run_interproscan(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    try:
        r = requests.post(f"{url}/run", data={
            "email": EMAIL, "sequence": sequence,
            "goterms": "true", "pathways": "true", "stype": "p",
        }, timeout=60)
        r.raise_for_status()
        jid = r.text.strip()
        if _ebi_poll(url, jid) != "FINISHED":
            raise RuntimeError("Job did not finish successfully")
        r = requests.get(f"{url}/result/{jid}/json",
                         headers={"Accept": "application/json"}, timeout=120)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.Timeout:
        raise RuntimeError(
            "InterProScan timed out — EBI servers may be slow. Try again in a few minutes."
        )


def _parse_blast_xml(xml_text: str) -> dict:
    """Parse EBI NCBI BLAST XML (EBIApplicationResult schema) into JSON2-compatible dict."""
    try:
        # Strip XML namespaces so we can use plain tag names in findall/iter
        import re as _re
        xml_clean = _re.sub(r'\s+xmlns(?::\w+)?="[^"]+"', '', xml_text)
        xml_clean = _re.sub(r'\s+xsi:\w+="[^"]+"', '', xml_clean)
        root = ET.fromstring(xml_clean)
        hits = []
        for hit in root.iter("hit"):
            acc  = hit.get("ac") or hit.get("id") or ""
            desc = hit.get("description") or ""
            # Extract organism from UniProt "OS=Name OX=..." field in description
            sciname = ""
            if " OS=" in desc:
                os_part = desc.split(" OS=", 1)[1]
                sciname = os_part.split(" OX=")[0].split(" GN=")[0].strip()
            # First alignment only
            aln = hit.find(".//alignment")
            if aln is None:
                continue
            try:
                evalue    = float(aln.findtext("expectation", "0") or "0")
                bit_score = float(aln.findtext("bits", "0") or "0")
                # EBI gives identity as a percentage (e.g. 47.9); encode as count/100
                # so that show_blast's formula (identity/align_len*100) gives the right %
                identity_pct = float(aln.findtext("identity", "0") or "0")
            except (ValueError, TypeError):
                continue
            hits.append({
                "description": [{"accession": acc, "title": desc, "sciname": sciname}],
                "hsps": [{"evalue": evalue, "bit_score": bit_score,
                          "identity": identity_pct, "align_len": 100}],
            })
        return {"BlastOutput2": [{"report": {"results": {"search": {"hits": hits}}}}]}
    except Exception:
        return {"BlastOutput2": []}


def run_blast(sequence):
    # Uses EBI NCBI BLAST Job Dispatcher — same reliable infrastructure as
    # InterProScan and Phobius, avoids direct NCBI QBLAST connectivity issues.
    url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"

    # SUBMISSION — retry up to 3 times
    jid = None
    for attempt in range(3):
        try:
            r = requests.post(f"{url}/run", data={
                "email": EMAIL, "program": "blastp",
                "database": "uniprotkb_bacteria",
                "sequence": sequence, "stype": "protein",
                "scores": 10, "alignments": 10, "exp": "1e-3",
            }, timeout=60)
            r.raise_for_status()
            jid = r.text.strip()
            if jid:
                break
        except Exception as exc:
            if attempt < 2:
                time.sleep(10)
            else:
                raise RuntimeError(f"BLASTp submission failed after 3 attempts: {exc}")
    if not jid:
        raise RuntimeError("BLASTp: no job ID in EBI response")

    print(f"[BLAST] Submitted EBI job={jid}, waiting 30s before first poll")
    time.sleep(30)

    # POLLING — 1200s total budget, 30s between polls (EBI fair-use), 30s per request
    start    = time.time()
    deadline = start + 1200
    status   = "PENDING"
    poll_n   = 0
    while time.time() < deadline:
        poll_n += 1
        elapsed = int(time.time() - start)
        try:
            r = requests.get(f"{url}/status/{jid}", timeout=30)
            r.raise_for_status()
            status = r.text.strip()
            print(f"[BLAST] Poll #{poll_n} at {elapsed}s elapsed: status={status}")
            if status == "FINISHED":
                break
            if status in ("FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
                raise RuntimeError(f"BLASTp job ended with status: {status}")
        except requests.exceptions.Timeout:
            print(f"[BLAST] Poll #{poll_n} request timed out (30s), retrying")
        remaining = deadline - time.time()
        if remaining > 0:
            time.sleep(min(30, remaining))

    if status != "FINISHED":
        raise RuntimeError(
            "BLASTp did not complete within 20 minutes — EBI queue may be busy. "
            "Try again later."
        )

    # RESULT DOWNLOAD — XML converted to JSON2-compatible dict
    try:
        r = requests.get(f"{url}/result/{jid}/xml", timeout=60)
        r.raise_for_status()
        return _parse_blast_xml(r.text)
    except requests.exceptions.Timeout:
        raise RuntimeError("BLASTp result download timed out — try again in a few minutes.")


def run_phobius(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/phobius"
    try:
        r = requests.post(f"{url}/run", data={
            "email": EMAIL, "sequence": sequence,
            "format": "short", "stype": "protein",
        }, timeout=60)
        r.raise_for_status()
        jid = r.text.strip()
        if _ebi_poll(url, jid) != "FINISHED":
            raise RuntimeError("Job did not finish successfully")
        r = requests.get(f"{url}/result/{jid}/out", timeout=90)
        r.raise_for_status()
        return r.text
    except requests.exceptions.Timeout:
        raise RuntimeError(
            "Phobius timed out — EBI servers may be slow. Try again in a few minutes."
        )


def run_hmmer(sequence):
    url = "https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan"

    # SUBMISSION — retry up to 3 times
    jid = None
    for attempt in range(3):
        try:
            r = requests.post(f"{url}/run", data={
                "email": EMAIL, "sequence": sequence,
                "database": "pfam", "E": "1.0",
            }, timeout=60)
            r.raise_for_status()
            jid = r.text.strip()
            if jid:
                break
        except Exception as exc:
            if attempt < 2:
                time.sleep(10)
            else:
                raise RuntimeError(f"HMMER submission failed after 3 attempts: {exc}")
    if not jid:
        raise RuntimeError("HMMER: no job ID in EBI response")

    print(f"[HMMER] Submitted job ID={jid}, waiting 5s before first poll")
    time.sleep(5)

    # POLLING — 360s total budget, 20s per individual request, 15s between polls
    start = time.time()
    deadline = start + 360
    status = "PENDING"
    poll_n = 0
    while time.time() < deadline:
        poll_n += 1
        elapsed = int(time.time() - start)
        try:
            r = requests.get(f"{url}/status/{jid}", timeout=20)
            r.raise_for_status()
            status = r.text.strip()
            print(f"[HMMER] Poll #{poll_n} at {elapsed}s elapsed: status={status}")
            if status == "FINISHED":
                break
            if status in ("FAILURE", "ERROR", "NOT_FOUND", "CANCELLED"):
                raise RuntimeError(f"HMMER job ended with status: {status}")
        except requests.exceptions.Timeout:
            print(f"[HMMER] Poll #{poll_n} request timed out (20s), retrying")
        remaining = deadline - time.time()
        if remaining > 0:
            time.sleep(min(15, remaining))

    if status != "FINISHED":
        raise RuntimeError(
            "HMMER did not complete within 6 minutes — EBI servers may be under load. "
            "Try again in a few minutes."
        )

    # RESULT DOWNLOAD
    try:
        r = requests.get(f"{url}/result/{jid}/out", timeout=60)
        r.raise_for_status()
        return r.text
    except requests.exceptions.Timeout:
        raise RuntimeError("HMMER result download timed out — try again in a few minutes.")


def run_foldseek(pdb_text):
    url = "https://search.foldseek.com/api"
    if not pdb_text.endswith("\n"):
        pdb_text += "\n"
    data = [("q", pdb_text), ("mode", "3diaa"), ("email", EMAIL)]
    for db in ["afdb50", "afdb-swissprot", "afdb-proteome", "BFVD"]:
        data.append(("database[]", db))

    # Submission — retry up to 3 times on 503 (server overload)
    resp = None
    for attempt in range(3):
        try:
            r = requests.post(f"{url}/ticket", data=data, timeout=60)
            if r.status_code == 503:
                if attempt < 2:
                    time.sleep(30)
                    continue
                r.raise_for_status()
            r.raise_for_status()
            resp = r.json()
            break
        except requests.exceptions.Timeout:
            if attempt < 2:
                time.sleep(30)
            else:
                raise RuntimeError("FoldSeek submission timed out after 3 attempts")
    if resp is None:
        raise RuntimeError("FoldSeek submission failed after 3 attempts")
    if resp.get("status") in ("RATELIMIT", "MAINTENANCE"):
        raise RuntimeError(f"FoldSeek server: {resp['status']}")
    ticket = resp["id"]

    while True:
        try:
            s = requests.get(f"{url}/ticket/{ticket}", timeout=60).json().get("status", "UNKNOWN")
        except requests.exceptions.HTTPError as e:
            if e.response is not None and e.response.status_code == 503:
                time.sleep(30)
                continue
            raise
        if s in ("COMPLETE", "ERROR", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)
    if s != "COMPLETE":
        raise RuntimeError(f"FoldSeek status: {s}")

    # Result download — retry up to 3 times on 503
    for attempt in range(3):
        try:
            r = requests.get(f"{url}/result/download/{ticket}", timeout=120)
            if r.status_code == 503:
                if attempt < 2:
                    time.sleep(30)
                    continue
                r.raise_for_status()
            r.raise_for_status()
            break
        except requests.exceptions.Timeout:
            if attempt < 2:
                time.sleep(30)
            else:
                raise RuntimeError("FoldSeek result download timed out after 3 attempts")
    r.raise_for_status()
    return r.content


# ── HTML rendering helpers ─────────────────────────────────────────────────────

def _esc(v: object) -> str:
    return _html.escape(str(v))


def _html_table(df: pd.DataFrame, max_rows: int = 500) -> str:
    if df.empty:
        return '<p style="color:#334155;font-size:13px;padding:16px 0;">No results found.</p>'
    df = df.head(max_rows)
    ths = "".join(f"<th>{_esc(c)}</th>" for c in df.columns)
    body = ""
    for _, row in df.iterrows():
        tds = "".join(
            f"<td title='{_esc(v)}'>"
            f"{_esc(str(v)[:72] if len(str(v)) > 72 else v)}</td>"
            for v in row
        )
        body += f"<tr>{tds}</tr>"
    return (
        f'<div class="dp-wrap"><table class="dp-tbl">'
        f"<thead><tr>{ths}</tr></thead><tbody>{body}</tbody>"
        f"</table></div>"
    )


def _cards(*items) -> str:
    """items: (label, value, sub_or_None) tuples"""
    inner = "".join(
        f'<div class="dp-card">'
        f'<div class="dp-card-label">{_esc(lbl)}</div>'
        f'<div class="dp-card-value">{_esc(str(val))}</div>'
        f'{"<div class=dp-card-sub>" + _esc(sub) + "</div>" if sub else ""}'
        f"</div>"
        for lbl, val, sub in items
    )
    return f'<div class="dp-cards">{inner}</div>'


# ── Result display functions ───────────────────────────────────────────────────

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
    df = pd.DataFrame(rows)
    n_db = df["Database"].nunique() if not df.empty else 0
    st.markdown(_cards(
        ("Matches",       len(df), None),
        ("Databases hit", n_db,    None),
    ), unsafe_allow_html=True)
    st.markdown(_html_table(df), unsafe_allow_html=True)


def show_blast(data: dict) -> None:
    parse_err = None
    hits = []
    try:
        bo2  = data["BlastOutput2"]
        bo2  = bo2[0] if isinstance(bo2, list) else bo2
        hits = bo2["report"]["results"]["search"]["hits"]
    except (KeyError, IndexError, TypeError) as e:
        parse_err = e

    if parse_err is not None or not hits:
        # Visible red diagnostic box — always shown when 0 results
        _diag_lines = []
        if parse_err:
            _diag_lines.append(f"Parse error: {parse_err}")
        _bo2 = data.get("BlastOutput2", "KEY MISSING")
        if isinstance(_bo2, list):
            _diag_lines.append(f"BlastOutput2 list length: {len(_bo2)}")
            if _bo2:
                _diag_lines.append(f"First entry keys: {list(_bo2[0].keys()) if isinstance(_bo2[0], dict) else type(_bo2[0])}")
        else:
            _diag_lines.append(f"BlastOutput2 value: {_bo2!r}")
        _diag_lines.append(f"Hit count: {len(hits)}")
        st.error(
            "BLASTp returned 0 hits.\n\n"
            + "\n".join(_diag_lines)
        )
        with st.expander("Raw BLAST response (first 4000 chars)"):
            try:
                st.code(json.dumps(data, indent=2)[:4000], language="json")
            except Exception:
                st.code(str(data)[:4000])
        if not hits:
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
    st.markdown(_cards(("Top hits", len(df), "uniprotkb_bacteria")), unsafe_allow_html=True)
    st.markdown(_html_table(df), unsafe_allow_html=True)


def show_phobius(text: str) -> None:
    for line in text.splitlines():
        if not line.strip() or line.startswith("SEQENCE"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            sp = "Yes" if parts[2] == "1" else "No"
            st.markdown(_cards(
                ("TM helices",    parts[1], None),
                ("Signal peptide", sp,       None),
                ("Topology",      parts[3], "o = outside · i = inside · h = TM helix"),
            ), unsafe_allow_html=True)
            return
    st.markdown(
        '<p style="color:#334155;font-size:13px;padding:16px 0;">No results parsed.</p>',
        unsafe_allow_html=True,
    )


def show_hmmer(text: str) -> None:
    rows, in_hits, below = [], False, False
    found_scores_section = False
    for line in text.splitlines():
        if "Scores for complete sequence" in line:
            in_hits = True
            found_scores_section = True
            continue
        if not in_hits:
            continue
        if "Domain annotation" in line:
            break
        s = line.strip()
        if not s or s.startswith("E-value"):
            continue
        # Must check inclusion threshold BEFORE the "---" guard — the line
        # "------ inclusion threshold ------" starts with "---" and would
        # otherwise be skipped, leaving `below` permanently False.
        if "inclusion threshold" in s:
            below = True
            continue
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
    df = pd.DataFrame(rows)

    if df.empty:
        if not text:
            # Genuinely empty response — keep as error for debugging
            st.error("HMMER returned an empty response (zero bytes).")
            return
        # 0 rows after parsing = valid "no hits" result (common for dark proteome proteins)
        st.markdown(_cards(
            ("Significant domains", 0, "above inclusion threshold"),
            ("Total",               0, None),
        ), unsafe_allow_html=True)
        st.info(
            "No Pfam domain matches found for this protein — this is common for "
            "dark proteome proteins with no characterised family."
        )
        return

    above = int((df["Significant"] == "Yes").sum())
    st.markdown(_cards(
        ("Significant domains",            above,   "above inclusion threshold"),
        ("Total (incl. below threshold)",  len(df), None),
    ), unsafe_allow_html=True)
    st.markdown(_html_table(df), unsafe_allow_html=True)


def show_foldseek(tar_bytes: bytes) -> None:
    base_cols = [
        "query", "target", "fident", "alnlen", "mismatch",
        "gapopen", "qstart", "qend", "tstart", "tend", "prob", "evalue", "bits",
    ]
    disp_cols = ["Target", "Identity %", "Probability", "E-value", "Bits", "Q.Start", "Q.End"]

    frames: dict   = {}
    db_names: list = []
    with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:gz") as tar:
        m8s = [m for m in tar.getmembers()
               if m.name.endswith(".m8") and "_report" not in m.name]
        if not m8s:
            st.markdown(
                '<p style="color:#334155;font-size:13px;padding:16px 0;">'
                'No result files in archive.</p>',
                unsafe_allow_html=True,
            )
            return
        for member in m8s:
            db = member.name.split("alis_")[-1].replace(".m8", "")
            db_names.append(db)
            f = tar.extractfile(member)
            if not f:
                frames[db] = pd.DataFrame(columns=disp_cols)
                continue
            lines = [
                l for l in f.read().decode("utf-8", errors="replace").splitlines()
                if l.strip()
            ]
            if not lines:
                frames[db] = pd.DataFrame(columns=disp_cols)
                continue
            rws = []
            for line in lines:
                p = line.split("\t")
                rws.append({c: (p[i] if i < len(p) else "") for i, c in enumerate(base_cols)})
            df = pd.DataFrame(rws)
            for c in ("fident", "prob", "evalue", "bits"):
                df[c] = pd.to_numeric(df[c], errors="coerce")
            df["target"] = df["target"].str[:60]
            df = df[["target", "fident", "prob", "evalue", "bits", "qstart", "qend"]].copy()
            df.columns = disp_cols
            frames[db] = df

    total = sum(len(frames.get(db, [])) for db in db_names)
    st.markdown(_cards(
        ("Total hits",        total,          "across all databases"),
        ("Databases searched", len(db_names), None),
    ), unsafe_allow_html=True)

    db_tabs = st.tabs(db_names)
    for tab, db in zip(db_tabs, db_names):
        with tab:
            df = frames.get(db, pd.DataFrame())
            st.markdown(
                f'<p style="color:#334155;font-size:11px;padding:4px 0 2px;">'
                f'{len(df)} hits</p>',
                unsafe_allow_html=True,
            )
            st.markdown(_html_table(df), unsafe_allow_html=True)


# ── Structure tab — 3D viewer template ────────────────────────────────────────
# Raw string so JS braces need no escaping; __SLOTS__ replaced at call time.

_VIEWER_TMPL = r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<style>
  html,body{margin:0;padding:0;background:#0a0e1a;overflow:hidden;}
  #viewer{width:100%;height:__HEIGHT__px;position:relative;}
  .brow{display:flex;gap:6px;padding:9px 12px;background:#080c17;
        border-top:1px solid #1e2d4a;align-items:center;}
  button{background:#0d1424;color:#64748b;border:1px solid #1e2d4a;border-radius:5px;
         padding:4px 11px;font-size:11px;cursor:pointer;
         font-family:-apple-system,sans-serif;transition:all .15s;}
  button:hover{background:#1e2d4a;color:#f0f6ff;}
  button.on{background:#3b82f6;color:#fff;border-color:#3b82f6;}
  .lbl{color:#334155;font-size:9.5px;font-weight:700;letter-spacing:.1em;
       text-transform:uppercase;margin-right:2px;}
</style>
</head>
<body>
<div id="viewer"></div>
<div class="brow">
  <span class="lbl">Colour</span>
  <button id="b0" class="on" onclick="cPLDDT()">pLDDT</button>
  <button id="b1"            onclick="cChain()">Chain</button>
  <button id="b2"            onclick="cHydro()">Hydrophobic</button>
  <button id="b3"            onclick="cCharge()">Charge</button>
</div>
<script>
var viewer=$3Dmol.createViewer('viewer',{backgroundColor:'#0a0e1a'});
viewer.addModel(atob("__PDB_B64__"),'pdb');
function act(id){
  ['b0','b1','b2','b3'].forEach(function(b){
    document.getElementById(b).classList.remove('on');
  });
  document.getElementById(id).classList.add('on');
}
function cPLDDT(){
  act('b0');
  viewer.setStyle({},{cartoon:{colorfunc:function(a){
    return a.b>90?'#3b82f6':a.b>70?'#eab308':a.b>50?'#f97316':'#ef4444';
  }}});
  viewer.render();
}
function cChain(){
  act('b1');
  viewer.setStyle({},{cartoon:{colorscheme:'chain'}});
  viewer.render();
}
function cHydro(){
  act('b2');
  var hp=['ILE','LEU','VAL','PHE','TRP','MET','ALA','TYR','CYS','PRO'];
  viewer.setStyle({},{cartoon:{colorfunc:function(a){
    return hp.indexOf(a.resn)>=0?'#f97316':'#3b82f6';
  }}});
  viewer.render();
}
function cCharge(){
  act('b3');
  var p=['LYS','ARG','HIS'],n=['ASP','GLU'];
  viewer.setStyle({},{cartoon:{colorfunc:function(a){
    if(p.indexOf(a.resn)>=0) return '#3b82f6';
    if(n.indexOf(a.resn)>=0) return '#ef4444';
    return '#475569';
  }}});
  viewer.render();
}
cPLDDT();
viewer.zoomTo();
viewer.zoom(1.1);
viewer.render();
</script>
</body>
</html>"""


def _3dmol_html(pdb_text: str, height: int = 440) -> str:
    pdb_b64 = base64.b64encode(pdb_text.encode()).decode()
    return (_VIEWER_TMPL
            .replace("__PDB_B64__", pdb_b64)
            .replace("__HEIGHT__", str(height)))


# ── Structure tab — analysis helpers ──────────────────────────────────────────

def _parse_fasta_seq(fasta_text: str) -> str:
    return "".join(
        line.strip() for line in fasta_text.splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()


def _parse_pdb_plddts(pdb_text: str) -> list:
    """Extract per-residue pLDDT scores from B-factor column (Cα atoms, first model)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", io.StringIO(pdb_text))
    scores = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    scores.append(res["CA"].get_bfactor())
        break
    return scores


def _run_dssp_raw(pdb_text: str):
    """Return per-residue SS code list ('H', 'E', or 'C'), or None on total failure.

    Strategy (in order):
      1. pydssp  — pure-Python DSSP if the package is installed.
      2. Phi/psi dihedral heuristic via BioPython — always available
         (BioPython is a hard dependency); less accurate than DSSP but
         never fails on valid PDB input.
      3. mkdssp / dssp binary via BioPython DSSP wrapper.
    """
    # ── 1. pydssp (pure Python, optional) ─────────────────────────────────
    try:
        import pydssp
        raw    = pydssp.read_pdbtext(pdb_text)
        coords = raw[0] if isinstance(raw, (list, tuple)) else raw
        ss_arr = pydssp.assign(coords, out_type="c3")
        return [str(s) for s in ss_arr]
    except Exception:
        pass

    # ── 2. Phi/psi dihedral heuristic (always available via BioPython) ────
    try:
        from Bio.PDB.Polypeptide import PPBuilder
        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure("p", io.StringIO(pdb_text))
        ppb       = PPBuilder()
        ss_list: list[str] = []
        for pp in ppb.build_peptides(structure):
            for phi, psi in pp.get_phi_psi_list():
                if phi is None or psi is None:
                    ss_list.append("C")
                    continue
                phi_d = math.degrees(phi)
                psi_d = math.degrees(psi)
                # α-helix region of Ramachandran plot
                if -90 <= phi_d <= -30 and -77 <= psi_d <= -17:
                    ss_list.append("H")
                # β-strand region (both parallel and antiparallel)
                elif (-170 <= phi_d <= -50 and
                      (100 <= psi_d <= 180 or -180 <= psi_d <= -150)):
                    ss_list.append("E")
                else:
                    ss_list.append("C")
        if ss_list:
            return ss_list
    except Exception:
        pass

    # ── 3. mkdssp / dssp binary via BioPython ─────────────────────────────
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", io.StringIO(pdb_text))
    model = structure[0]
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        pdbio.save(f)
        tmppath = f.name
    try:
        for exe in ("mkdssp", "dssp"):
            try:
                if _DSSP is None:
                    continue
                dssp_obj = _DSSP(model, tmppath, dssp=exe)
                return [dssp_obj[k][2] for k in dssp_obj]
            except Exception:
                continue
        return None
    except Exception:
        return None
    finally:
        os.unlink(tmppath)


def _detect_disordered(scores: list, threshold: float = 50.0, min_len: int = 5) -> list:
    regions, in_r, start = [], False, 0
    for i, s in enumerate(scores):
        if s < threshold and not in_r:
            in_r, start = True, i + 1
        elif s >= threshold and in_r:
            in_r = False
            if (i - start) >= min_len:
                regions.append((start, i))
    if in_r and (len(scores) - start) >= min_len:
        regions.append((start, len(scores)))
    return regions


def _detect_low_complexity(seq: str, window: int = 20, threshold: float = 0.65) -> list:
    if len(seq) < window:
        return []
    max_ent = math.log2(min(window, 20))
    regions, in_r, start = [], False, 0
    for i in range(len(seq) - window + 1):
        w = seq[i : i + window]
        cnt: dict = {}
        for aa in w:
            cnt[aa] = cnt.get(aa, 0) + 1
        norm = (-sum((c / window) * math.log2(c / window) for c in cnt.values()) / max_ent
                if max_ent else 1)
        if norm < threshold and not in_r:
            in_r, start = True, i + 1
        elif norm >= threshold and in_r:
            in_r = False
            regions.append((start, i + window))
    if in_r:
        regions.append((start, len(seq)))
    return regions


# ── Structure tab — HTML sub-components ───────────────────────────────────────

def _plddt_bar_html(scores: list) -> str:
    n = len(scores)
    bands = [
        (">90 Very high",   "#3b82f6", sum(1 for s in scores if s > 90)),
        ("70–90 Confident", "#eab308", sum(1 for s in scores if 70 < s <= 90)),
        ("50–70 Low",       "#f97316", sum(1 for s in scores if 50 < s <= 70)),
        ("<50 Very low",    "#ef4444", sum(1 for s in scores if s <= 50)),
    ]
    segs = "".join(
        f'<div style="flex:{c};background:{col};height:100%;min-width:3px;" '
        f'title="{lbl}: {round(c/n*100,1) if n else 0}%"></div>'
        for lbl, col, c in bands if c > 0
    )
    legend = "".join(
        f'<span style="display:inline-flex;align-items:center;gap:4px;">'
        f'<span style="width:9px;height:9px;border-radius:2px;background:{col};'
        f'display:inline-block;flex-shrink:0;"></span>'
        f'<span style="color:#475569;font-size:10.5px;">'
        f'{lbl} ({round(c/n*100,1) if n else 0}%)</span></span>'
        for lbl, col, c in bands
    )
    return (
        '<div style="margin-top:14px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        'text-transform:uppercase;margin:0 0 8px;">pLDDT Confidence</p>'
        '<div style="display:flex;height:12px;border-radius:4px;overflow:hidden;'
        f'border:1px solid #1e2d4a;">{segs}</div>'
        '<div style="display:flex;flex-wrap:wrap;gap:10px;margin-top:8px;">'
        f'{legend}</div></div>'
    )


def _ss_bars_html(ss: dict) -> str:
    items = [
        ("Helix",       "#3b82f6", ss.get("H", 0)),
        ("Strand",      "#22c55e", ss.get("E", 0)),
        ("Loop / coil", "#475569", ss.get("C", 0)),
    ]
    rows = "".join(
        f'<div style="margin-bottom:8px;">'
        f'<div style="display:flex;justify-content:space-between;margin-bottom:3px;">'
        f'<span style="color:#94a3b8;font-size:12px;">{lbl}</span>'
        f'<span style="color:#f0f6ff;font-size:12px;font-weight:600;">{pct}%</span></div>'
        f'<div style="background:#0d1424;border-radius:3px;height:5px;">'
        f'<div style="background:{col};width:{pct}%;height:100%;border-radius:3px;"></div>'
        f'</div></div>'
        for lbl, col, pct in items
    )
    return (
        '<div style="margin-top:18px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        f'text-transform:uppercase;margin:0 0 10px;">Secondary Structure</p>{rows}</div>'
    )


def _aa_comp_html(aa_comp: dict) -> str:
    top10 = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)[:10]
    mx = top10[0][1] if top10 else 1
    rows = "".join(
        f'<div style="display:flex;align-items:center;gap:7px;margin-bottom:4px;">'
        f'<span style="color:#64748b;font-size:11px;font-family:monospace;'
        f'width:14px;text-align:center;flex-shrink:0;">{aa}</span>'
        f'<div style="flex:1;background:#0d1424;border-radius:2px;height:13px;">'
        f'<div style="background:#3b82f6;width:{round(frac/mx*100,1)}%;height:100%;'
        f'border-radius:2px;opacity:0.75;"></div></div>'
        f'<span style="color:#64748b;font-size:10.5px;width:34px;text-align:right;'
        f'flex-shrink:0;">{round(frac*100,1)}%</span></div>'
        for aa, frac in top10
    )
    return (
        '<div style="margin-top:18px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        f'text-transform:uppercase;margin:0 0 10px;">Amino Acid Composition (top 10)</p>'
        f'{rows}</div>'
    )


def _features_html(scores: list, seq: str, phobius_text) -> str:
    feats = []
    if scores:
        dis = _detect_disordered(scores)
        if dis:
            ex = ", ".join(f"{s}–{e}" for s, e in dis[:3])
            if len(dis) > 3:
                ex += f" (+{len(dis)-3} more)"
            feats.append(("#f97316", "Disordered regions (pLDDT < 50)",
                           f"{len(dis)} region(s): residues {ex}"))
    if seq:
        lc = _detect_low_complexity(seq)
        if lc:
            ex = ", ".join(f"{s}–{e}" for s, e in lc[:3])
            if len(lc) > 3:
                ex += f" (+{len(lc)-3} more)"
            feats.append(("#eab308", "Low-complexity regions",
                           f"{len(lc)} region(s): residues {ex}"))
    if phobius_text:
        for line in phobius_text.splitlines():
            if not line.strip() or line.startswith("SEQENCE"):
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[1].isdigit():
                tm = int(parts[1])
                if tm > 0:
                    feats.append(("#22c55e", "Transmembrane helices (Phobius)",
                                   f"{tm} TM helix{'es' if tm > 1 else ''} predicted"))
                break
    if not feats:
        body = ('<p style="color:#334155;font-size:12px;padding:2px 0;">'
                'No notable structural features detected.</p>')
    else:
        body = "".join(
            f'<div style="background:#080c17;border:1px solid #1e2d4a;border-radius:7px;'
            f'padding:10px 14px;margin-bottom:8px;border-left:3px solid {col};">'
            f'<div style="color:{col};font-size:11.5px;font-weight:600;margin-bottom:2px;">'
            f'{title}</div>'
            f'<div style="color:#475569;font-size:11px;">{_esc(detail)}</div></div>'
            for col, title, detail in feats
        )
    return (
        '<div style="margin-top:18px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        f'text-transform:uppercase;margin:0 0 10px;">Detected Structural Features</p>'
        f'{body}</div>'
    )


# ── Fold Type Classifier ──────────────────────────────────────────────────────

def _run_dssp_strands(pdb_text: str):
    """Return (strand_list, pct_dict) or ([], None).
    strand_list: [(start, end), ...] 1-based indices of beta-strand runs."""
    ss_list = _run_dssp_raw(pdb_text)
    if ss_list is None:
        return [], None
    strands: list = []
    in_s, s0 = False, 0
    for i, ss in enumerate(ss_list):
        beta = ss in ("E", "B")
        if beta and not in_s:
            in_s, s0 = True, i + 1
        elif not beta and in_s:
            in_s = False
            strands.append((s0, i))
    if in_s:
        strands.append((s0, len(ss_list)))
    counts = {"H": 0, "E": 0, "C": 0}
    for ss in ss_list:
        if ss in ("H", "G", "I"): counts["H"] += 1
        elif ss in ("E", "B"):    counts["E"] += 1
        else:                     counts["C"] += 1
    total = sum(counts.values())
    pct = {k: round(v / total * 100, 1) for k, v in counts.items()} if total else None
    return strands, pct


def _parse_pdb_ca_coords(pdb_text: str) -> list:
    """Return [(x, y, z)] for Cα atoms in chain/residue order (first model)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", io.StringIO(pdb_text))
    coords = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    v = res["CA"].get_vector()
                    coords.append((float(v[0]), float(v[1]), float(v[2])))
        break
    return coords


# Accessions that mark the C-terminal β-barrel translocator in autotransporters
_BARREL_ACCESSIONS = {"PS51208", "SSF103515", "PF00291"}

# Accessions that positively identify autotransporter passenger (β-helix) domain
_AT_ACCESSIONS = {"PF20696", "PF27415", "PF22399", "PS51208", "PF00291"}

# BLAST keywords for Type Va autotransporters (fha/FHA moved to _FHA_BLAST_KEYWORDS)
_AT_BLAST_KEYWORDS = [
    "biga", "autotransporter", "aida", "trimeric autotransporter", "type v secretion",
    "pertactin", "intimin", "adhesin involved in",
]

# RHS repeat solenoid — very large proteins, highly irregular spacing, N-terminal coiled-coil
_RHS_ACCESSIONS = {"PF03917", "PF05593", "PF19434", "PF13796"}
_RHS_BLAST_KEYWORDS = [
    "rhs repeat", "rhs element", "type vi secretion",
    "vgrg", "hemolysin-coregulated", "wxg100",
]

# FHA-type beta-helix — TPS/Two-Partner Secretion; filamentous hemagglutinin family
_FHA_ACCESSIONS = {"PF02413", "PF13312", "PF15609"}
_FHA_BLAST_KEYWORDS = [
    "filamentous hemagglutinin", "filamentous haemagglutinin", "fha",
    "two-partner secretion", "tps domain",
]

# FN3/Ig beta-sandwich — tandem immunoglobulin/fibronectin type III repeats
_FN3_IG_ACCESSIONS = {
    "PF00041", "PF00047", "PF07654", "PF13895", "PF13927", "PF09661", "PF13833",
}
_FN3_IG_BLAST_KEYWORDS = [
    "fibronectin type iii", "fibronectin-type iii", "immunoglobulin-like",
    "ig-like fold", "fn3 repeat",
]

# TPR solenoid — alpha-alpha superhelix; tetratricopeptide repeats; all-helix, no strands
_TPR_ACCESSIONS = {
    "PF00515", "PF07719", "PF07720", "PF07721",
    "PF13428", "PF13432", "PF13512", "PF14559",
}
_TPR_BLAST_KEYWORDS = [
    "tetratricopeptide", "tpr repeat", "sel1 repeat", "hal repeat",
]

# TIM barrel — 8-stranded parallel beta-barrel with alternating outer helices; 200–450 aa
_TIM_ACCESSIONS = {"SSF51351", "G3DSA:3.20.20.70"}
_TIM_BLAST_KEYWORDS = [
    "tim barrel", "triosephosphate isomerase", "aldolase", "enolase",
    "xylose isomerase", "glucoamylase", "phosphotriesterase",
]

# Rossmann fold — parallel beta-sheet alternating with helices; nucleotide/cofactor binding
_ROSSMANN_ACCESSIONS = {"SSF52374", "SSF51735", "G3DSA:3.40.50.720"}
_ROSSMANN_BLAST_KEYWORDS = [
    "rossmann", "nad-binding", "nad binding", "nadh", "nadph",
    "nucleotide-binding", "dehydrogenase", "reductase", "oxidoreductase",
]

# OB fold — small oligonucleotide/oligosaccharide-binding 5-stranded beta-barrel
_OB_ACCESSIONS = {"SSF50249", "PF01336"}
_OB_BLAST_KEYWORDS = [
    "oligonucleotide binding", "ob fold", "ob-fold",
    "cold shock", "ribosomal protein s1",
]

# Beta-propeller — radially symmetric WD40/Kelch/YWTD repeat; many short antiparallel strands
_PROPELLER_ACCESSIONS = {"PF00400", "PF01344", "PF01966", "SSF50978"}
_PROPELLER_BLAST_KEYWORDS = [
    "wd40", "wd repeat", "wd-40", "kelch repeat", "kelch-type",
    "beta-propeller", "ywtd",
]

# Lectin/jelly-roll beta-sandwich — ConA-type; carbohydrate-binding antiparallel beta-sheet
_LECTIN_ACCESSIONS = {"SSF49899", "PF00139", "PF00652", "PF14200"}
_LECTIN_BLAST_KEYWORDS = [
    "lectin", "carbohydrate-binding", "concanavalin", "galectin",
    "jelly roll", "jelly-roll", "sugar-binding",
]

# Coiled-coil — heptad-repeat helical bundle; >70% helix, <5% strand
_COIL_ACCESSIONS = {"PF05765", "SSF47370"}
_COIL_BLAST_KEYWORDS = [
    "coiled-coil", "coiled coil", "leucine zipper", "coil domain",
]

# Alpha-solenoid (HEAT/ARM) — stacked helical hairpin repeats; large, high helix content
_HEAT_ACCESSIONS = {"PF02985", "PF00514", "PF13190", "PF19012"}
_HEAT_BLAST_KEYWORDS = [
    "heat repeat", "heat-repeat", "armadillo repeat", "arm repeat",
    "importin", "exportin", "beta-catenin",
]

# Calycin/lipocalin — small 8-stranded antiparallel beta-barrel cup; lipid/ligand carrier
_CALYCIN_ACCESSIONS = {"PF00061", "SSF50814"}
_CALYCIN_BLAST_KEYWORDS = [
    "lipocalin", "fatty acid-binding", "fatty acid binding",
    "retinol-binding", "calycin", "carrier protein",
]


def _blast_top_title(blast_data) -> str:
    """Return the lowercase title of the top BLAST hit, or ''."""
    if not blast_data:
        return ""
    for result in blast_data.get("BlastOutput2", []):
        search = result.get("report", {}).get("results", {}).get("search", {})
        hits   = search.get("hits", [])
        if hits:
            return hits[0].get("description", [{}])[0].get("title", "").lower()
    return ""


def _ipr_hits_for(ipr_data, accession_set: set) -> list:
    """Return list of (acc, lib, name) for each InterProScan hit in accession_set."""
    if not ipr_data:
        return []
    found, seen = [], set()
    for res in ipr_data.get("results", []):
        for match in res.get("matches", []):
            sig = match.get("signature", {})
            acc = sig.get("accession", "")
            if acc in accession_set and acc not in seen:
                seen.add(acc)
                name = sig.get("name", "") or sig.get("description", "")
                lib  = sig.get("signatureLibraryRelease", {}).get("library", "")
                found.append((acc, lib, name))
    return found


def _find_barrel_cutoff(ipr_data) -> tuple:
    """Return (cutoff_residue, hit_label) for the earliest C-terminal β-barrel hit."""
    if not ipr_data:
        return 0, None
    cutoff, hit_label = None, None
    for res in ipr_data.get("results", []):
        for match in res.get("matches", []):
            sig = match.get("signature", {})
            acc = sig.get("accession", "")
            if acc in _BARREL_ACCESSIONS:
                name = sig.get("name", "") or sig.get("description", "")
                lib  = sig.get("signatureLibraryRelease", {}).get("library", "")
                for loc in match.get("locations", []):
                    s = int(loc.get("start", 0))
                    if s > 0 and (cutoff is None or s < cutoff):
                        cutoff = s
                        hit_label = f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})"
    return (cutoff or 0), hit_label


def _score_autotransporter(ipr_data, blast_data, seq: str,
                            passenger_gap_sd: float) -> tuple:
    """Return (score, evidence_list) for the Autotransporter β-helix fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _AT_ACCESSIONS):
        score = min(100, score + 25)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
    top = _blast_top_title(blast_data)
    for kw in _AT_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 15)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    if seq and seq.rstrip()[-1].upper() == "F":
        score = min(100, score + 10)
        hits.append("C-terminal Phe (BAM complex recognition signal)")
    if 0 < passenger_gap_sd < 8:
        score = min(100, score + 15)
        hits.append(f"Passenger strand SD {passenger_gap_sd:.1f} < 8 (regular β-helix rungs)")
    return min(100, score), hits


def _score_rhs(ipr_data, blast_data, seq: str, gap_sd: float,
               parallel_pct: int, antiparallel_pct: int,
               first_strand_start: int) -> tuple:
    """Return (score, evidence_list) for the RHS repeat solenoid fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _RHS_ACCESSIONS):
        score = min(100, score + 30)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 60:
            break
    if seq and len(seq) > 1500:
        score = min(100, score + 20)
        hits.append(f"Large protein ({len(seq)} aa > 1500 aa RHS threshold)")
    if gap_sd > 20:
        score = min(100, score + 15)
        hits.append(f"Highly irregular strand spacing (SD={gap_sd:.1f} > 20)")
    if parallel_pct < 60 and antiparallel_pct < 60:
        score = min(100, score + 10)
        hits.append(f"Mixed strand orientation ({parallel_pct}% par, {antiparallel_pct}% antipar)")
    if first_strand_start > 200:
        score = min(100, score + 10)
        hits.append(f"N-terminal non-strand region (first strand at residue {first_strand_start})")
    top = _blast_top_title(blast_data)
    for kw in _RHS_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 15)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_fha(ipr_data, blast_data, seq: str, gap_sd: float) -> tuple:
    """Return (score, evidence_list) for the FHA-type β-helix (TPS) fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _FHA_ACCESSIONS):
        score = min(100, score + 25)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 75:
            break
    if seq and len(seq) > 2000:
        score = min(100, score + 15)
        hits.append(f"Very long protein ({len(seq)} aa > 2000 aa FHA threshold)")
    if 0 < gap_sd < 8:
        score = min(100, score + 10)
        hits.append(f"Semi-regular strand spacing (SD={gap_sd:.1f} < 8, like AT β-helix)")
    top = _blast_top_title(blast_data)
    for kw in _FHA_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 15)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_fn3_ig(ipr_data, blast_data, n_strands: int,
                  antiparallel_pct: int) -> tuple:
    """Return (score, evidence_list) for the FN3/Ig tandem-repeat sandwich fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _FN3_IG_ACCESSIONS)[:4]:
        score = min(100, score + 20)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
    if n_strands >= 6:
        score = min(100, score + 10)
        hits.append(f"{n_strands} beta strands — consistent with tandem Ig/FN3 domains")
    if antiparallel_pct > 50:
        score = min(100, score + 10)
        hits.append(f"{antiparallel_pct}% antiparallel — Ig/FN3 beta-sheets are antiparallel")
    top = _blast_top_title(blast_data)
    for kw in _FN3_IG_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 10)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_tpr(ipr_data, blast_data, ss_pct) -> tuple:
    """Return (score, evidence_list) for the TPR alpha-alpha solenoid fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _TPR_ACCESSIONS):
        score = min(100, score + 25)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 75:
            break
    if ss_pct and ss_pct.get("H", 0) > 40:
        score = min(100, score + 15)
        hits.append(f"Predominantly helical ({ss_pct['H']}% helix by DSSP)")
    if ss_pct and ss_pct.get("E", 100) < 10:
        score = min(100, score + 10)
        hits.append(f"Very low beta content ({ss_pct.get('E', 0)}% strand by DSSP)")
    top = _blast_top_title(blast_data)
    for kw in _TPR_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 15)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_tim_barrel(ipr_data, blast_data, seq, n_strands,
                      parallel_pct, mean_gap, ss_pct) -> tuple:
    """Return (score, evidence_list) for the TIM barrel fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _TIM_ACCESSIONS):
        score = min(100, score + 35)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 70:
            break
    if 6 <= n_strands <= 10:
        score = min(100, score + 20)
        hits.append(f"{n_strands} beta strands — TIM barrel has 8 ±2 parallel strands")
    if parallel_pct > 60:
        score = min(100, score + 15)
        hits.append(f"{parallel_pct}% parallel — TIM barrel beta-strands are parallel")
    if 10 <= mean_gap <= 45:
        score = min(100, score + 10)
        hits.append(f"Mean inter-strand gap {mean_gap:.0f} aa — consistent with helix-containing loop")
    if seq and 200 <= len(seq) <= 450:
        score = min(100, score + 10)
        hits.append(f"Protein size {len(seq)} aa — within TIM barrel range (200–450 aa)")
    if ss_pct and ss_pct.get("H", 0) > 20:
        score = min(100, score + 5)
        hits.append(f"Significant helix content ({ss_pct['H']}%) — TIM barrel has 8 outer helices")
    top = _blast_top_title(blast_data)
    for kw in _TIM_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 15)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_rossmann(ipr_data, blast_data, ss_pct, parallel_pct) -> tuple:
    """Return (score, evidence_list) for the Rossmann fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _ROSSMANN_ACCESSIONS):
        score = min(100, score + 35)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 70:
            break
    if ss_pct and ss_pct.get("H", 0) > 25 and ss_pct.get("E", 0) > 15:
        score = min(100, score + 15)
        hits.append(f"Mixed α/β (H={ss_pct['H']}%, E={ss_pct['E']}%) — Rossmann fold is α/β")
    if parallel_pct > 50:
        score = min(100, score + 10)
        hits.append(f"{parallel_pct}% parallel strands — Rossmann fold has parallel beta-sheet")
    top = _blast_top_title(blast_data)
    for kw in _ROSSMANN_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_ob_fold(ipr_data, blast_data, seq, n_strands, ss_pct) -> tuple:
    """Return (score, evidence_list) for the OB fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _OB_ACCESSIONS):
        score = min(100, score + 40)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 80:
            break
    if seq and len(seq) < 200:
        score = min(100, score + 20)
        hits.append(f"Small protein ({len(seq)} aa < 200 aa — OB folds are small domains)")
    if ss_pct and ss_pct.get("E", 0) > 30:
        score = min(100, score + 10)
        hits.append(f"High beta content ({ss_pct['E']}% strand) — OB fold is predominantly beta")
    if 4 <= n_strands <= 7:
        score = min(100, score + 10)
        hits.append(f"{n_strands} beta strands — OB fold typically has 5 beta-strands")
    top = _blast_top_title(blast_data)
    for kw in _OB_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_beta_propeller(ipr_data, blast_data, n_strands, ss_pct) -> tuple:
    """Return (score, evidence_list) for the beta-propeller fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _PROPELLER_ACCESSIONS):
        score = min(100, score + 30)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 60:
            break
    if n_strands >= 16:
        score = min(100, score + 20)
        hits.append(f"{n_strands} beta strands — propellers have many short blades (≥16 for 4-blade)")
    if ss_pct and ss_pct.get("E", 0) > 25 and ss_pct.get("H", 0) < 20:
        score = min(100, score + 15)
        hits.append(f"Predominantly beta (E={ss_pct['E']}%, H={ss_pct['H']}%)")
    top = _blast_top_title(blast_data)
    for kw in _PROPELLER_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_lectin(ipr_data, blast_data, ss_pct, antiparallel_pct) -> tuple:
    """Return (score, evidence_list) for the lectin/jelly-roll beta-sandwich fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _LECTIN_ACCESSIONS):
        score = min(100, score + 35)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 70:
            break
    if antiparallel_pct > 50:
        score = min(100, score + 15)
        hits.append(f"{antiparallel_pct}% antiparallel — lectin jelly-roll is antiparallel beta")
    if ss_pct and ss_pct.get("E", 0) > 25 and ss_pct.get("H", 0) < 20:
        score = min(100, score + 10)
        hits.append(f"Predominantly beta (E={ss_pct['E']}%, H={ss_pct['H']}%)")
    top = _blast_top_title(blast_data)
    for kw in _LECTIN_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_coiled_coil(ipr_data, blast_data, ss_pct) -> tuple:
    """Return (score, evidence_list) for the coiled-coil fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _COIL_ACCESSIONS):
        score = min(100, score + 35)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 70:
            break
    if ss_pct and ss_pct.get("H", 0) > 70:
        score = min(100, score + 25)
        hits.append(f"Predominantly helical ({ss_pct['H']}% helix > 70% coiled-coil threshold)")
    if ss_pct and ss_pct.get("E", 100) < 5:
        score = min(100, score + 15)
        hits.append(f"Negligible beta content ({ss_pct.get('E', 0)}% strand by DSSP)")
    top = _blast_top_title(blast_data)
    for kw in _COIL_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_heat_arm(ipr_data, blast_data, seq, ss_pct) -> tuple:
    """Return (score, evidence_list) for the HEAT/ARM alpha-solenoid fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _HEAT_ACCESSIONS):
        score = min(100, score + 30)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 60:
            break
    if ss_pct and ss_pct.get("H", 0) > 60:
        score = min(100, score + 20)
        hits.append(f"Predominantly helical ({ss_pct['H']}% helix > 60% HEAT/ARM threshold)")
    if seq and len(seq) > 500:
        score = min(100, score + 15)
        hits.append(f"Large protein ({len(seq)} aa > 500 aa — HEAT/ARM solenoids are typically large)")
    if ss_pct and ss_pct.get("E", 100) < 15:
        score = min(100, score + 10)
        hits.append(f"Low beta content ({ss_pct.get('E', 0)}% strand) — HEAT/ARM have few strands")
    top = _blast_top_title(blast_data)
    for kw in _HEAT_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _score_calycin(ipr_data, blast_data, seq, n_strands, antiparallel_pct) -> tuple:
    """Return (score, evidence_list) for the calycin/lipocalin beta-barrel cup fold."""
    score, hits = 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _CALYCIN_ACCESSIONS):
        score = min(100, score + 45)
        hits.append(f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})")
        if score >= 90:
            break
    if seq and len(seq) < 200:
        score = min(100, score + 15)
        hits.append(f"Small protein ({len(seq)} aa < 200 aa — lipocalins are compact domains)")
    if antiparallel_pct > 60:
        score = min(100, score + 10)
        hits.append(f"{antiparallel_pct}% antiparallel — lipocalin barrel is antiparallel")
    if 6 <= n_strands <= 10:
        score = min(100, score + 10)
        hits.append(f"{n_strands} beta strands — lipocalin has 8 antiparallel strands")
    top = _blast_top_title(blast_data)
    for kw in _CALYCIN_BLAST_KEYWORDS:
        if kw in top:
            score = min(100, score + 20)
            hits.append(f"BLAST top hit: '{kw}'")
            break
    return min(100, score), hits


def _analyze_fold_type(strands: list, ca_coords: list, plddt: list, seq: str,
                       dssp_available: bool = True,
                       ipr_data=None, blast_data=None, ss_pct=None) -> dict:
    rtx_hits = _rtx_matches(seq) if seq else []
    n_total  = len(strands)

    # Helix-dominant folds can be scored without any beta strands
    tpr_score,  tpr_hits  = _score_tpr(ipr_data, blast_data, ss_pct)
    coil_score, coil_hits = _score_coiled_coil(ipr_data, blast_data, ss_pct)
    heat_score, heat_hits = _score_heat_arm(ipr_data, blast_data, seq, ss_pct)

    if n_total < 2:
        helix_scores = {
            "TPR solenoid":    tpr_score,
            "Coiled-coil":     coil_score,
            "HEAT/ARM solenoid": heat_score,
        }
        best_helix       = max(helix_scores, key=helix_scores.get)
        best_helix_score = helix_scores[best_helix]

        if best_helix_score >= 25:
            plddt_mean = round(sum(plddt) / len(plddt), 1) if plddt else 0
            _zero_fs: dict = {k: 0 for k in [
                "Beta-solenoid", "RHS solenoid", "Autotransporter β-helix", "FHA β-helix",
                "FN3/Ig sandwich", "Lectin/jelly-roll", "Beta-barrel", "Beta-propeller",
                "TIM barrel", "OB fold", "Calycin/lipocalin", "Rossmann fold",
                "Coiled-coil", "HEAT/ARM solenoid", "TPR solenoid",
            ]}
            _zero_fs.update({"Coiled-coil": coil_score,
                              "HEAT/ARM solenoid": heat_score,
                              "TPR solenoid": tpr_score})
            _best_hits = {"TPR solenoid": tpr_hits,
                          "Coiled-coil": coil_hits,
                          "HEAT/ARM solenoid": heat_hits}[best_helix]
            _criteria = {
                "All-helix architecture":    (ss_pct.get("E", 100) < 10) if ss_pct else True,
                f"{best_helix} domain hits": best_helix_score >= 25,
                "High pLDDT in repeat core": plddt_mean > 70,
                "No/few beta strands":       n_total == 0,
            }
            _conclusion = {
                "TPR solenoid":      "Evidence is consistent with a TPR-type alpha-alpha superhelix.",
                "Coiled-coil":       "Evidence is consistent with a heptad-repeat coiled-coil bundle.",
                "HEAT/ARM solenoid": "Evidence is consistent with a HEAT/ARM stacked-helix solenoid.",
            }
            _reasoning = [
                f"DSSP found {n_total} beta strand(s) — consistent with an all-helix architecture.",
                f"SS composition: {ss_pct or 'unavailable'}.",
                f"{best_helix} score {best_helix_score}/100: "
                + ("; ".join(_best_hits) or "no specific evidence") + ".",
                f"Best fold match: {best_helix} ({best_helix_score}/100). "
                + _conclusion.get(best_helix, "FoldSeek recommended to confirm."),
            ]
            return {
                "available": True, "n_strands": n_total, "n_passenger": 0,
                "strands": [], "barrel_cutoff": 0, "barrel_hit": None,
                "lengths": [], "mean_len": 0.0, "gaps": [], "mean_gap": 0.0,
                "gap_sd": 0.0, "angles": [], "mean_angle": 0.0, "angle_sd": 0.0,
                "parallel_pct": 0, "antiparallel_pct": 0,
                "rtx_hits": rtx_hits, "rtx_on_strands": False,
                "strand_plddt_mean": plddt_mean,
                "best_fold": best_helix, "criteria": _criteria,
                "reasoning": _reasoning, "fold_scores": _zero_fs,
                "dssp_available": dssp_available,
            }
        return {"available": False, "n_strands": n_total,
                "dssp_available": dssp_available, "rtx_hits": rtx_hits}

    # ── Exclude C-terminal β-barrel translocator from passenger geometry ──────
    barrel_cutoff, barrel_hit = _find_barrel_cutoff(ipr_data)
    if barrel_cutoff > 0:
        work_strands = [(s, e) for s, e in strands if s < barrel_cutoff]
        work_ca      = ca_coords[:barrel_cutoff - 1] if ca_coords else []
        work_plddt   = plddt[:barrel_cutoff - 1]     if plddt     else plddt
    else:
        work_strands, work_ca, work_plddt = strands, ca_coords, plddt

    if len(work_strands) < 2:
        work_strands, work_ca, work_plddt = strands, ca_coords, plddt
        barrel_cutoff, barrel_hit = 0, None

    n                  = len(work_strands)
    first_strand_start = work_strands[0][0] if work_strands else 0
    lengths            = [e - s + 1 for s, e in work_strands]
    mean_len           = sum(lengths) / n
    gaps               = [work_strands[i+1][0] - work_strands[i][1] - 1 for i in range(n - 1)]
    mean_gap           = sum(gaps) / len(gaps) if gaps else 0
    gap_sd             = (math.sqrt(sum((g - mean_gap)**2 for g in gaps) / len(gaps))
                          if len(gaps) > 1 else 0)

    def _norm(v):
        m = math.sqrt(sum(x*x for x in v))
        return tuple(x/m for x in v) if m > 0 else (0.0, 0.0, 0.0)

    strand_dirs = []
    for s, e in work_strands:
        s0, e0 = s - 1, e - 1
        if 0 <= s0 < len(work_ca) and 0 <= e0 < len(work_ca):
            v = tuple(work_ca[e0][j] - work_ca[s0][j] for j in range(3))
            strand_dirs.append(_norm(v))
        else:
            strand_dirs.append(None)

    angles = []
    for i in range(len(strand_dirs) - 1):
        v1, v2 = strand_dirs[i], strand_dirs[i+1]
        if v1 and v2:
            cos_a = max(-1.0, min(1.0, sum(a*b for a, b in zip(v1, v2))))
            angles.append(math.degrees(math.acos(cos_a)))

    mean_angle       = sum(angles) / len(angles) if angles else 0
    angle_sd         = (math.sqrt(sum((a - mean_angle)**2 for a in angles) / len(angles))
                        if len(angles) > 1 else 0)
    parallel_pct     = round(sum(1 for a in angles if a < 45)  / len(angles) * 100) if angles else 0
    antiparallel_pct = round(sum(1 for a in angles if a > 135) / len(angles) * 100) if angles else 0

    rtx_on_strands = any(
        not (re < s - 3 or rs > e + 3)
        for rs, re in rtx_hits
        for s, e in work_strands
    )

    strand_pl = [work_plddt[i]
                 for s, e in work_strands
                 for i in range(s - 1, min(e, len(work_plddt)))]
    strand_plddt_mean = round(sum(strand_pl) / len(strand_pl), 1) if strand_pl else 0

    # ── Score all 15 fold types ───────────────────────────────────────────────
    at_score,   at_hits   = _score_autotransporter(ipr_data, blast_data, seq, gap_sd)
    rhs_score,  rhs_hits  = _score_rhs(ipr_data, blast_data, seq, gap_sd,
                                        parallel_pct, antiparallel_pct, first_strand_start)
    fha_score,  fha_hits  = _score_fha(ipr_data, blast_data, seq, gap_sd)
    fn3_ig_score, fn3_ig_hits = _score_fn3_ig(ipr_data, blast_data, n_total, antiparallel_pct)
    tim_score,  tim_hits  = _score_tim_barrel(ipr_data, blast_data, seq, n_total,
                                               parallel_pct, mean_gap, ss_pct)
    ros_score,  ros_hits  = _score_rossmann(ipr_data, blast_data, ss_pct, parallel_pct)
    ob_score,   ob_hits   = _score_ob_fold(ipr_data, blast_data, seq, n_total, ss_pct)
    prop_score, prop_hits = _score_beta_propeller(ipr_data, blast_data, n_total, ss_pct)
    lec_score,  lec_hits  = _score_lectin(ipr_data, blast_data, ss_pct, antiparallel_pct)
    cal_score,  cal_hits  = _score_calycin(ipr_data, blast_data, seq, n_total, antiparallel_pct)

    sol  = min(30, n * 4)
    sol += max(0, int(30 - gap_sd * 6))
    sol += round(max(parallel_pct, antiparallel_pct) * 0.2) if angles else 0
    sol += 15 if rtx_on_strands else 0
    sol += 5  if strand_plddt_mean > 70 else 0
    sol  = min(100, max(0, sol))

    barrel  = 30 if 8 <= n_total <= 24 else 0
    barrel += int(antiparallel_pct * 0.25)
    barrel += 20 if 1 <= mean_gap <= 5 else 0
    barrel += max(0, int(20 - gap_sd * 3))
    barrel += 15 if 5 <= mean_len <= 12 else 0
    barrel  = min(100, max(0, barrel))

    fold_scores = {
        "Beta-solenoid":           sol,
        "RHS solenoid":            rhs_score,
        "Autotransporter β-helix": at_score,
        "FHA β-helix":             fha_score,
        "FN3/Ig sandwich":         fn3_ig_score,
        "Lectin/jelly-roll":       lec_score,
        "Beta-barrel":             barrel,
        "Beta-propeller":          prop_score,
        "TIM barrel":              tim_score,
        "OB fold":                 ob_score,
        "Calycin/lipocalin":       cal_score,
        "Rossmann fold":           ros_score,
        "Coiled-coil":             coil_score,
        "HEAT/ARM solenoid":       heat_score,
        "TPR solenoid":            tpr_score,
    }
    best_fold = max(fold_scores, key=fold_scores.get)

    # ── Context-sensitive criteria (4 cards, tailored to best fold) ──────────
    _crit_defaults: dict = {
        "Regular strand spacing":    gap_sd < 4.0,
        "Parallel orientation":      parallel_pct > 60,
        "High pLDDT in repeat core": strand_plddt_mean > 70,
        "RTX sequence motif":        rtx_on_strands,
    }
    _crit_overrides: dict = {
        "Beta-solenoid": {**_crit_defaults},
        "RHS solenoid": {
            "RHS repeat domain hits":       rhs_score >= 25,
            "Large protein (>1500 aa)":     bool(seq and len(seq) > 1500),
            "Highly irregular spacing":     gap_sd > 20,
            "N-terminal non-strand region": first_strand_start > 200,
        },
        "Autotransporter β-helix": {
            "Regular strand spacing":       gap_sd < 4.0,
            "Parallel orientation":         parallel_pct > 60,
            "High pLDDT in repeat core":    strand_plddt_mean > 70,
            "Autotransporter domain hits":  at_score >= 25,
        },
        "FHA β-helix": {
            "TPS/ESPR domain hits":         fha_score >= 25,
            "Very long protein (>2000 aa)": bool(seq and len(seq) > 2000),
            "Semi-regular spacing (SD<8)":  0 < gap_sd < 8,
            "High pLDDT in repeat core":    strand_plddt_mean > 70,
        },
        "FN3/Ig sandwich": {
            "Multiple Ig/FN3 domain hits":  fn3_ig_score >= 20,
            "Antiparallel orientation":     antiparallel_pct > 50,
            "High pLDDT in repeat core":    strand_plddt_mean > 70,
            "≥6 beta strands":              n_total >= 6,
        },
        "Lectin/jelly-roll": {
            "Lectin/ConA domain hits":      lec_score >= 25,
            "Antiparallel orientation":     antiparallel_pct > 50,
            "Predominantly beta":           bool(ss_pct and ss_pct.get("E", 0) > 25),
            "High pLDDT":                   strand_plddt_mean > 70,
        },
        "Beta-barrel": {
            "8–24 beta strands":            8 <= n_total <= 24,
            "Antiparallel orientation":     antiparallel_pct > 50,
            "Short strands (5–12 aa)":      5 <= mean_len <= 12,
            "Tight inter-strand gaps":      1 <= mean_gap <= 5,
        },
        "Beta-propeller": {
            "WD40/Kelch domain hits":       prop_score >= 20,
            "≥16 beta strands":             n_total >= 16,
            "Predominantly beta":           bool(ss_pct and ss_pct.get("E", 0) > 25),
            "Low helix content":            bool(ss_pct and ss_pct.get("H", 0) < 20),
        },
        "TIM barrel": {
            "TIM barrel IPR hits":          tim_score >= 20,
            "8 ±2 parallel strands":        6 <= n_total <= 10,
            "Helix-containing gaps (10–45 aa)": 10 <= mean_gap <= 45,
            "Protein size 200–450 aa":      bool(seq and 200 <= len(seq) <= 450),
        },
        "OB fold": {
            "OB fold IPR hits":             ob_score >= 30,
            "Small protein (<200 aa)":      bool(seq and len(seq) < 200),
            "Predominantly beta":           bool(ss_pct and ss_pct.get("E", 0) > 30),
            "4–7 beta strands":             4 <= n_total <= 7,
        },
        "Calycin/lipocalin": {
            "Lipocalin IPR hits (PF00061)": cal_score >= 30,
            "Small protein (<200 aa)":      bool(seq and len(seq) < 200),
            "Antiparallel orientation":     antiparallel_pct > 60,
            "6–10 beta strands":            6 <= n_total <= 10,
        },
        "Rossmann fold": {
            "Rossmann/NAD-binding IPR hits": ros_score >= 25,
            "Parallel beta-sheet":          parallel_pct > 50,
            "Mixed α/β content":            bool(ss_pct and ss_pct.get("H", 0) > 25
                                                  and ss_pct.get("E", 0) > 15),
            "Nucleotide-binding keyword":   any(kw in _blast_top_title(blast_data)
                                                for kw in _ROSSMANN_BLAST_KEYWORDS),
        },
        "Coiled-coil": {
            ">70% helix by DSSP":           bool(ss_pct and ss_pct.get("H", 0) > 70),
            "<5% strand by DSSP":           bool(ss_pct and ss_pct.get("E", 100) < 5),
            "Coiled-coil IPR hits":         coil_score >= 25,
            "High pLDDT":                   strand_plddt_mean > 70,
        },
        "HEAT/ARM solenoid": {
            "HEAT/ARM domain hits":         heat_score >= 25,
            ">60% helix by DSSP":           bool(ss_pct and ss_pct.get("H", 0) > 60),
            "Large protein (>500 aa)":      bool(seq and len(seq) > 500),
            "Low beta content (<15%)":      bool(ss_pct and ss_pct.get("E", 100) < 15),
        },
        "TPR solenoid": {
            "TPR domain hits":              tpr_score >= 25,
            ">40% helix by DSSP":           bool(ss_pct and ss_pct.get("H", 0) > 40),
            "<10% strand by DSSP":          bool(ss_pct and ss_pct.get("E", 100) < 10),
            "High pLDDT in repeat core":    strand_plddt_mean > 70,
        },
    }
    criteria = _crit_overrides.get(best_fold, _crit_defaults)

    # ── Reasoning chain ───────────────────────────────────────────────────────
    spacing_q = ("regular" if gap_sd < 4 else
                 "moderately regular" if gap_sd < 8 else
                 "highly irregular" if gap_sd > 20 else "irregular")
    orient_q  = (f"{parallel_pct}% parallel" if parallel_pct > 60
                 else f"{antiparallel_pct}% antiparallel" if antiparallel_pct > 60
                 else "mixed orientation")
    reasoning: list = []

    if barrel_cutoff > 0:
        reasoning.append(
            f"β-barrel translocator detected ({barrel_hit}) at residue {barrel_cutoff}. "
            f"Passenger domain (res 1–{barrel_cutoff - 1}) used for geometry statistics."
        )
    if ss_pct:
        reasoning.append(
            f"Global SS composition: {ss_pct.get('H', 0)}% helix, "
            f"{ss_pct.get('E', 0)}% strand, {ss_pct.get('C', 0)}% coil."
        )
    reasoning.append(
        f"DSSP: {n_total} beta strands ({n} in passenger domain), "
        f"mean {round(mean_len, 1)} aa (range {min(lengths)}–{max(lengths)} aa)."
    )
    reasoning.append(
        f"Inter-strand spacing: mean {round(mean_gap, 1)} aa, SD {round(gap_sd, 1)} — {spacing_q}. "
        f"Orientation: {orient_q}."
    )
    if first_strand_start > 200:
        reasoning.append(
            f"First strand at residue {first_strand_start} — substantial N-terminal "
            f"non-strand region (RHS coiled-coil/WXG100, signal domain, or unstructured linker)."
        )
    # Evidence lines only for folds with non-zero scores
    _evidence_map = [
        ("Autotransporter β-helix", at_score,   at_hits),
        ("RHS solenoid",            rhs_score,  rhs_hits),
        ("FHA β-helix",             fha_score,  fha_hits),
        ("FN3/Ig sandwich",         fn3_ig_score, fn3_ig_hits),
        ("Lectin/jelly-roll",       lec_score,  lec_hits),
        ("TIM barrel",              tim_score,  tim_hits),
        ("Rossmann fold",           ros_score,  ros_hits),
        ("OB fold",                 ob_score,   ob_hits),
        ("Beta-propeller",          prop_score, prop_hits),
        ("Calycin/lipocalin",       cal_score,  cal_hits),
        ("Coiled-coil",             coil_score, coil_hits),
        ("HEAT/ARM solenoid",       heat_score, heat_hits),
        ("TPR solenoid",            tpr_score,  tpr_hits),
    ]
    for fname, fscore, fhits in _evidence_map:
        if fhits:
            reasoning.append(f"{fname} evidence ({fscore}/100): " + "; ".join(fhits) + ".")
    if rtx_on_strands:
        reasoning.append("RTX nonapeptide (GGXGXDXUX) detected at strand edges.")
    reasoning.append(
        f"Mean pLDDT over strand residues: {strand_plddt_mean} "
        f"({'high' if strand_plddt_mean > 70 else 'moderate/low'})."
    )
    _conclusion_map = {
        "Beta-solenoid":           "RTX-type parallel beta-solenoid.",
        "RHS solenoid":            "RHS-type solenoid with N-terminal coiled-coil/WXG100 domain.",
        "Autotransporter β-helix": "Autotransporter right-handed β-helix (Type Va).",
        "FHA β-helix":             "FHA-type β-helix; Two-Partner Secretion (TPS).",
        "FN3/Ig sandwich":         "Tandem FN3/Ig-like beta-sandwich repeat domains.",
        "Lectin/jelly-roll":       "Lectin/jelly-roll ConA-type antiparallel beta-sandwich.",
        "Beta-barrel":             "Transmembrane beta-barrel.",
        "Beta-propeller":          "Radially symmetric WD40/Kelch-type beta-propeller.",
        "TIM barrel":              "8-stranded parallel TIM barrel (beta/alpha)8 fold.",
        "OB fold":                 "OB-fold oligonucleotide/oligosaccharide-binding barrel.",
        "Calycin/lipocalin":       "Calycin/lipocalin 8-stranded antiparallel beta-barrel cup.",
        "Rossmann fold":           "Rossmann-fold nucleotide/cofactor-binding alpha/beta domain.",
        "Coiled-coil":             "Heptad-repeat coiled-coil helical bundle.",
        "HEAT/ARM solenoid":       "HEAT/ARM stacked-helix alpha-solenoid.",
        "TPR solenoid":            "TPR-type alpha-alpha superhelix with no beta strands.",
    }
    best_sc = fold_scores[best_fold]
    reasoning.append(
        f"Best fold match: {best_fold} ({best_sc}/100). "
        + _conclusion_map.get(best_fold, "FoldSeek recommended to confirm assignment.")
    )

    return {
        "available":          True,
        "n_strands":          n_total,
        "n_passenger":        n,
        "strands":            work_strands,
        "barrel_cutoff":      barrel_cutoff,
        "barrel_hit":         barrel_hit,
        "lengths":            lengths,
        "mean_len":           round(mean_len, 1),
        "gaps":               gaps,
        "mean_gap":           round(mean_gap, 1),
        "gap_sd":             round(gap_sd, 1),
        "angles":             [round(a, 1) for a in angles],
        "mean_angle":         round(mean_angle, 1),
        "angle_sd":           round(angle_sd, 1),
        "parallel_pct":       parallel_pct,
        "antiparallel_pct":   antiparallel_pct,
        "rtx_hits":           rtx_hits,
        "rtx_on_strands":     rtx_on_strands,
        "strand_plddt_mean":  strand_plddt_mean,
        "best_fold":          best_fold,
        "criteria":           criteria,
        "reasoning":          reasoning,
        "fold_scores":        fold_scores,
        "dssp_available":     dssp_available,
    }


def _solenoid_strand_viz_html(analysis: dict, seq_len: int) -> str:
    strands  = analysis["strands"]
    angles   = analysis.get("angles", [])
    rtx_hits = analysis.get("rtx_hits", [])
    if not strands:
        return ""

    W, H   = 660, 70
    LM, RM = 8, 8
    TW     = W - LM - RM
    Y      = 40
    AH     = 14

    def sx(pos):
        return LM + (pos - 1) / max(seq_len - 1, 1) * TW

    palette = ["#3b82f6", "#22c55e"]
    els: list[str] = []

    # Backbone rail
    els.append(
        f'<line x1="{LM}" y1="{Y}" x2="{LM+TW}" y2="{Y}" '
        f'stroke="#1e2d4a" stroke-width="1.5"/>'
    )

    # Orientation arcs between consecutive strands
    for i, angle in enumerate(angles):
        if i + 1 >= len(strands):
            break
        x1  = sx(strands[i][1])
        x2  = sx(strands[i+1][0])
        mid = (x1 + x2) / 2
        is_par = angle < 45
        col = "#3b82f6" if is_par else "#f97316"
        arc_y = Y - 18
        label = "Parallel" if is_par else "Antiparallel"
        els.append(
            f'<path d="M{x1:.1f},{Y-AH//2} Q{mid:.1f},{arc_y} {x2:.1f},{Y-AH//2}" '
            f'fill="none" stroke="{col}" stroke-width="1.2" stroke-dasharray="3,2" opacity="0.7">'
            f'<title>{label} ({angle:.0f}°)</title></path>'
        )

    # Strand arrow bars
    for i, (s, e) in enumerate(strands):
        col = palette[i % 2]
        x1, x2 = sx(s), sx(e)
        w = max(x2 - x1, 3.0)
        els.append(
            f'<rect x="{x1:.1f}" y="{Y-AH//2}" width="{w:.1f}" height="{AH}" '
            f'rx="2" fill="{col}" opacity="0.85">'
            f'<title>Strand {i+1}: residues {s}–{e} ({e-s+1} aa)</title></rect>'
        )
        if w > 8:
            ax = min(x1 + w, LM + TW - 1)
            ah2 = AH // 2
            els.append(
                f'<polygon points="{ax:.1f},{Y-ah2} {min(ax+5, LM+TW):.1f},{Y} {ax:.1f},{Y+ah2}" '
                f'fill="{col}" opacity="0.85"/>'
            )
        if w > 15:
            els.append(
                f'<text x="{(x1+x2)/2:.1f}" y="{Y+4.5}" text-anchor="middle" '
                f'font-size="7" fill="#f0f6ff" font-family="monospace" '
                f'pointer-events="none">{i+1}</text>'
            )

    # RTX motif diamonds above the backbone
    for rs, re in rtx_hits:
        mx = sx((rs + re) / 2)
        els.append(
            f'<polygon points="{mx:.1f},{Y-AH-8} {mx+4:.1f},{Y-AH-3} '
            f'{mx:.1f},{Y-AH+2} {mx-4:.1f},{Y-AH-3}" fill="#a855f7" opacity="0.9">'
            f'<title>RTX motif {rs}–{re}</title></polygon>'
        )

    # Ruler ticks
    tick_step = max(10, round(seq_len / 8 / 10) * 10)
    for t in list(range(1, seq_len + 1, tick_step)) + [seq_len]:
        tx     = sx(t)
        anchor = "start" if t == 1 else ("end" if t == seq_len else "middle")
        els.append(
            f'<line x1="{tx:.1f}" y1="{Y+AH//2+2}" x2="{tx:.1f}" y2="{Y+AH//2+6}" '
            f'stroke="#1e2d4a" stroke-width="1"/>'
        )
        els.append(
            f'<text x="{tx:.1f}" y="{Y+AH//2+15}" text-anchor="{anchor}" '
            f'font-size="8" fill="#334155" font-family="monospace">{t}</text>'
        )

    # Legend row
    lx = LM
    ly = 10
    for lbl, col in [("Strand A", "#3b82f6"), ("Strand B", "#22c55e")]:
        els.append(f'<rect x="{lx}" y="{ly-7}" width="10" height="8" rx="1" fill="{col}"/>')
        els.append(
            f'<text x="{lx+13}" y="{ly}" font-size="7.5" fill="#475569" '
            f'font-family="-apple-system,sans-serif">{lbl}</text>'
        )
        lx += 65
    for lbl, col in [("Parallel", "#3b82f6"), ("Antiparallel", "#f97316")]:
        els.append(
            f'<line x1="{lx}" y1="{ly-4}" x2="{lx+14}" y2="{ly-4}" '
            f'stroke="{col}" stroke-width="1.3" stroke-dasharray="3,2"/>'
        )
        els.append(
            f'<text x="{lx+17}" y="{ly}" font-size="7.5" fill="#475569" '
            f'font-family="-apple-system,sans-serif">{lbl}</text>'
        )
        lx += 80
    if rtx_hits:
        els.append(
            f'<polygon points="{lx},{ly-7} {lx+4},{ly-3} {lx},{ly+1} {lx-4},{ly-3}" '
            f'fill="#a855f7"/>'
        )
        els.append(
            f'<text x="{lx+7}" y="{ly}" font-size="7.5" fill="#475569" '
            f'font-family="-apple-system,sans-serif">RTX motif</text>'
        )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="100%" style="display:block;">'
        f'<rect width="{W}" height="{H}" fill="#080c17"/>'
        f'{"".join(els)}</svg>'
    )


def _fold_scores_html(fold_scores: dict, best: str) -> str:
    rows = ""
    for fold, score in fold_scores.items():
        is_best = fold == best
        col     = "#3b82f6" if is_best else "#334155"
        lc      = "#f0f6ff" if is_best else "#64748b"
        glow    = f"box-shadow:0 0 6px {col}66;" if is_best else ""
        rows += (
            f'<div style="margin-bottom:10px;">'
            f'<div style="display:flex;justify-content:space-between;margin-bottom:3px;">'
            f'<span style="color:{lc};font-size:12px;font-weight:{"600" if is_best else "400"};">'
            f'{_esc(fold)}{" ✓" if is_best else ""}</span>'
            f'<span style="color:{lc};font-size:12px;font-weight:600;">{score}/100</span></div>'
            f'<div style="background:#0d1424;border-radius:3px;height:6px;">'
            f'<div style="background:{col};width:{score}%;height:100%;border-radius:3px;{glow}">'
            f'</div></div></div>'
        )
    return (
        '<div style="margin-top:16px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        f'text-transform:uppercase;margin:0 0 10px;">Fold Type Score</p>{rows}</div>'
    )


def _criteria_cards_html(criteria: dict) -> str:
    items = ""
    for label, passed in criteria.items():
        icon = "✓" if passed else "✗"
        col  = "#22c55e" if passed else "#ef4444"
        bg   = "rgba(34,197,94,0.07)"  if passed else "rgba(239,68,68,0.07)"
        bdr  = "rgba(34,197,94,0.2)"   if passed else "rgba(239,68,68,0.2)"
        items += (
            f'<div style="background:{bg};border:1px solid {bdr};border-radius:7px;'
            f'padding:8px 12px;display:flex;align-items:center;gap:8px;">'
            f'<span style="color:{col};font-size:14px;font-weight:700;flex-shrink:0;">{icon}</span>'
            f'<span style="color:#94a3b8;font-size:11px;">{_esc(label)}</span></div>'
        )
    return (
        '<div style="margin-top:16px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        'text-transform:uppercase;margin:0 0 10px;">Classification Criteria</p>'
        f'<div style="display:grid;grid-template-columns:1fr 1fr;gap:6px;">{items}</div></div>'
    )


def _reasoning_chain_html(reasoning: list) -> str:
    items = "".join(
        f'<li style="color:#64748b;font-size:11px;line-height:1.65;margin-bottom:7px;">'
        f'{_esc(s)}</li>'
        for s in reasoning
    )
    return (
        '<div style="margin-top:16px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        f'text-transform:uppercase;margin:0 0 10px;">Reasoning Chain</p>'
        f'<ol style="margin:0;padding-left:18px;">{items}</ol></div>'
    )


def _fold_classifier_html(analysis: dict, seq_len: int) -> str:
    if not analysis.get("available"):
        n        = analysis.get("n_strands", 0)
        dssp_ok  = analysis.get("dssp_available", True)
        rtx_hits = analysis.get("rtx_hits", [])

        if not dssp_ok:
            msg = ("Secondary structure assignment failed — could not assign "
                   "strands from this PDB file. Sequence-only results shown below.")
        elif n < 2:
            msg = (f"Only {n} beta strand(s) detected by DSSP — "
                   f"insufficient for solenoid scoring (need ≥2).")
        else:
            msg = "Beta-solenoid analysis unavailable."

        # RTX motif box — shown even when DSSP is missing
        rtx_html = ""
        if rtx_hits:
            positions = ", ".join(f"{s}–{e}" for s, e in rtx_hits[:6])
            extra = f" (+{len(rtx_hits)-6} more)" if len(rtx_hits) > 6 else ""
            rtx_html = (
                '<div style="background:rgba(168,85,247,0.08);border:1px solid rgba(168,85,247,0.25);'
                'border-radius:7px;padding:10px 14px;margin-top:10px;">'
                '<div style="color:#a855f7;font-size:11.5px;font-weight:600;margin-bottom:4px;">'
                f'✓ {len(rtx_hits)} RTX nonapeptide motif(s) detected</div>'
                '<div style="color:#94a3b8;font-size:11px;">'
                f'Positions: {_esc(positions + extra)}</div>'
                '<div style="color:#64748b;font-size:10.5px;margin-top:3px;">'
                'GGXGXDXUX pattern — characteristic of RTX-family beta-solenoid toxins.</div>'
                '</div>'
            )
        else:
            rtx_html = (
                '<div style="background:rgba(71,85,105,0.08);border:1px solid rgba(71,85,105,0.2);'
                'border-radius:7px;padding:10px 14px;margin-top:10px;">'
                '<div style="color:#475569;font-size:11.5px;">'
                '✗ No RTX nonapeptide motif (GGXGXDXUX) detected in sequence.</div></div>'
            )

        return (
            '<div style="margin-top:16px;">'
            '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
            'text-transform:uppercase;margin:0 0 8px;">Fold Type Classifier</p>'
            f'<p style="color:#64748b;font-size:11.5px;margin:0 0 2px;">{_esc(msg)}</p>'
            f'{rtx_html}</div>'
        )

    best     = analysis["best_fold"]
    best_sc  = analysis["fold_scores"][best]
    conf_col = "#22c55e" if best_sc > 60 else "#f59e0b" if best_sc > 35 else "#ef4444"
    conf_lbl = ("High" if best_sc > 60 else "Moderate" if best_sc > 35 else "Low") + f" confidence · {best_sc}/100"

    header = (
        '<div style="margin-top:16px;">'
        '<div style="display:flex;align-items:center;justify-content:space-between;'
        'flex-wrap:wrap;gap:6px;margin-bottom:10px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        'text-transform:uppercase;margin:0;">Fold Type Classifier</p>'
        f'<span style="background:{conf_col}22;border:1px solid {conf_col}55;'
        f'border-radius:999px;padding:2px 9px;font-size:10px;color:{conf_col};font-weight:600;">'
        f'{conf_lbl} · Best: {_esc(best)}</span></div>'
    )

    viz      = _solenoid_strand_viz_html(analysis, seq_len)
    scores_h = _fold_scores_html(analysis["fold_scores"], best)
    crit_h   = _criteria_cards_html(analysis["criteria"])
    chain_h  = _reasoning_chain_html(analysis["reasoning"])

    return (
        header
        + f'<div style="background:#080c17;border:1px solid #1e2d4a;border-radius:8px;'
          f'padding:10px 8px 4px;overflow-x:auto;">{viz}</div>'
        + scores_h + crit_h + chain_h + "</div>"
    )


# ── Structure tab — main display function ─────────────────────────────────────

def show_structure(pdb_text: str, fasta_text, phobius_text=None) -> None:
    plddt   = _parse_pdb_plddts(pdb_text)
    mean_pl = round(sum(plddt) / len(plddt), 1) if plddt else 0
    seq     = _parse_fasta_seq(fasta_text) if fasta_text else ""

    seq_props = None
    if seq:
        try:
            ana = ProteinAnalysis(seq)
            seq_props = {
                "length":  len(seq),
                "mw_kda":  round(ana.molecular_weight() / 1000, 1),
                "pi":      round(ana.isoelectric_point(), 2),
                "aa_comp": ana.get_amino_acids_percent(),
            }
        except Exception:
            pass

    # Single DSSP run — yields strand segments + SS percentages
    # ss is None when the dssp/mkdssp binary is unavailable
    strands, ss = _run_dssp_strands(pdb_text)
    dssp_ok = ss is not None

    # Cα coordinates for crossing-angle analysis
    ca_coords = _parse_pdb_ca_coords(pdb_text) if strands else []

    # Pull InterProScan + BLAST results from session state if available
    _results    = st.session_state.get("results", {})
    _ipr_data   = (_results.get("InterProScan", {}) or {}).get("data")
    _blast_data = (_results.get("BLASTp", {}) or {}).get("data")

    # Beta-solenoid analysis
    sol_analysis = _analyze_fold_type(strands, ca_coords, plddt, seq,
                                     dssp_available=dssp_ok,
                                     ipr_data=_ipr_data,
                                     blast_data=_blast_data,
                                     ss_pct=ss)
    seq_len      = len(seq) or len(plddt) or 1

    col_v, col_a = st.columns([5, 4])

    with col_v:
        components.html(_3dmol_html(pdb_text, height=440), height=494, scrolling=False)
        if plddt:
            st.markdown(_plddt_bar_html(plddt), unsafe_allow_html=True)
        st.markdown(_fold_classifier_html(sol_analysis, seq_len), unsafe_allow_html=True)

    with col_a:
        if seq_props:
            st.markdown(_cards(
                ("Amino acids", seq_props["length"],           None),
                ("Mol. weight", f'{seq_props["mw_kda"]} kDa', None),
                ("pI",          seq_props["pi"],               None),
                ("Mean pLDDT",  mean_pl,                       None),
            ), unsafe_allow_html=True)
        else:
            st.markdown(_cards(("Mean pLDDT", mean_pl, None)), unsafe_allow_html=True)

        if ss:
            st.markdown(_ss_bars_html(ss), unsafe_allow_html=True)
        else:
            st.markdown(
                '<p style="color:#334155;font-size:11px;margin:14px 0 0;">'
                'Secondary structure unavailable — DSSP binary not found.</p>',
                unsafe_allow_html=True,
            )

        if seq_props:
            st.markdown(_aa_comp_html(seq_props["aa_comp"]), unsafe_allow_html=True)

        st.markdown(_features_html(plddt, seq, phobius_text), unsafe_allow_html=True)


# ── Domain Map tab ────────────────────────────────────────────────────────────

_DB_COLORS: dict = {
    "PFAM":             "#3b82f6",
    "PANTHER":          "#22c55e",
    "PROSITE_PATTERNS": "#a855f7",
    "PROSITE_PROFILES": "#a855f7",
    "SUPERFAMILY":      None,
    "GENE3D":           "#0891b2",
    "TIGRFAM":          "#f59e0b",
    "PRINTS":           "#ec4899",
    "PIRSF":            "#6366f1",
    "HAMAP":            "#14b8a6",
    "SFLD":             "#f97316",
    "SMART":            "#84cc16",
    "CDD":              "#fb7185",
}

_TW = 850    # track width (SVG units)
_LM = 135    # left margin
_RM = 15     # right margin
_TH = 20     # track height
_RH = 22     # row height per domain row


def _parse_ipr_domains(data) -> list:
    rows = []
    for res in data.get("results", []):
        for match in res.get("matches", []):
            sig   = match.get("signature", {})
            entry = sig.get("entry") or {}
            lib   = sig.get("signatureLibraryRelease", {}).get("library", "")
            desc  = entry.get("description") or sig.get("description") or ""
            name  = sig.get("name", "")
            acc   = sig.get("accession", "")
            ev    = match.get("evalue")
            for loc in match.get("locations", []):
                rows.append({
                    "db":    lib,
                    "acc":   acc,
                    "name":  name or acc,
                    "desc":  desc,
                    "start": int(loc.get("start", 0)),
                    "end":   int(loc.get("end", 0)),
                    "evalue": ev,
                })
    return sorted(rows, key=lambda d: d["start"])


def _dssp_segments(pdb_text: str) -> list:
    """Return [(start, end, type)] where type in H/E/C/D (D = pLDDT < 50)."""
    ss_list = _run_dssp_raw(pdb_text)
    plddts  = _parse_pdb_plddts(pdb_text)
    n = max(len(ss_list or []), len(plddts))
    if n == 0:
        return []
    def _t(i):
        bf = plddts[i] if i < len(plddts) else 100.0
        ss = ss_list[i] if ss_list and i < len(ss_list) else "C"
        if bf < 50:               return "D"
        if ss in ("H", "G", "I"): return "H"
        if ss in ("E", "B"):      return "E"
        return "C"
    segs, start, cur = [], 1, _t(0)
    for i in range(1, n):
        t = _t(i)
        if t != cur:
            segs.append((start, i, cur))
            start, cur = i + 1, t
    segs.append((start, n, cur))
    return segs


_SS_COLORS = {"H": "#3b82f6", "E": "#22c55e", "C": "#334155", "D": "#ef4444"}
_SS_LABELS = {"H": "Helix", "E": "Strand", "C": "Loop", "D": "Disordered"}


def _plddt_segments(scores: list) -> list:
    if not scores:
        return []
    def _col(s):
        return "#3b82f6" if s > 90 else "#eab308" if s > 70 else "#f97316" if s > 50 else "#ef4444"
    segs, start, cur = [], 1, _col(scores[0])
    for i, s in enumerate(scores[1:], start=2):
        c = _col(s)
        if c != cur:
            segs.append((start, i - 1, cur))
            start, cur = i, c
    segs.append((start, len(scores), cur))
    return segs


def _rtx_matches(seq: str) -> list:
    """Detect RTX nonapeptide GGXGXDXUX (U = hydrophobic: ILVCFYWM)."""
    hydro = set("ILVCFYWM")
    out = []
    for i in range(len(seq) - 8):
        w = seq[i : i + 9]
        if w[0] == "G" and w[1] == "G" and w[3] == "G" and w[5] == "D" and w[7] in hydro:
            out.append((i + 1, i + 9))
    return out


def _sliding_repeats(seq: str, min_unit: int = 8, max_unit: int = 50,
                     min_repeats: int = 3) -> list:
    """Detect tandem repeats by sliding-window unit similarity (≥60% identity)."""
    n = len(seq)
    found = []
    for unit in range(min_unit, min(max_unit + 1, n // min_repeats + 1)):
        i = 0
        while i + unit * min_repeats <= n:
            ref = seq[i : i + unit]
            count, j = 1, i + unit
            while j + unit <= n:
                cand = seq[j : j + unit]
                if sum(a == b for a, b in zip(ref, cand)) / unit >= 0.6:
                    count += 1
                    j += unit
                else:
                    break
            if count >= min_repeats:
                found.append((i + 1, j, count, unit))
                i = j
            else:
                i += 1
    found.sort()
    merged: list = []
    for r in found:
        if merged and r[0] <= merged[-1][1]:
            if (r[1] - r[0]) > (merged[-1][1] - merged[-1][0]):
                merged[-1] = r
        else:
            merged.append(r)
    return merged


def _svg_x(pos: int, seq_len: int) -> float:
    return _LM + (pos - 1) / max(seq_len - 1, 1) * _TW


def _svg_w(start: int, end: int, seq_len: int) -> float:
    return max((end - start + 1) / max(seq_len, 1) * _TW, 3.0)


def _svg_label(y: float, text: str) -> str:
    return (
        f'<text x="{_LM - 8}" y="{y + _TH / 2 + 4.5}" text-anchor="end" '
        f'font-family="-apple-system,monospace" font-size="9.5" font-weight="700" '
        f'fill="#3b82f6">{_esc(text)}</text>'
    )


def _domain_map_svg(seq_len: int, domains: list, ss_segs: list,
                    plddt_segs: list, repeats: list, rtx: list) -> str:
    if seq_len < 1:
        return ""
    els: list[str] = []
    y = 20

    # Track 1 — InterProScan domains
    if domains:
        rows: list[list] = []
        for d in domains:
            placed = False
            for row in rows:
                if all(d["start"] > r["end"] + 2 or d["end"] < r["start"] - 2 for r in row):
                    row.append(d)
                    placed = True
                    break
            if not placed:
                rows.append([d])

        els.append(_svg_label(y + (len(rows) * _RH) / 2 - _RH / 2, "InterProScan"))
        for ri, row in enumerate(rows):
            ry = y + ri * _RH
            els.append(
                f'<rect x="{_LM}" y="{ry + _RH/2 - 1}" width="{_TW}" height="2" '
                f'fill="#1e2d4a" rx="1"/>'
            )
            for d in row:
                col = _DB_COLORS.get(d["db"], "#475569")
                xd  = _svg_x(d["start"], seq_len)
                wd  = _svg_w(d["start"], d["end"], seq_len)
                ev  = f'{d["evalue"]:.2e}' if d["evalue"] is not None else "n/a"
                tip = f'{_esc(d["name"])} ({d["db"]})\nResidues {d["start"]}–{d["end"]}\nE-value: {ev}'
                fill_attr = f'fill="{col}" opacity="0.85"' if col else 'fill="none" stroke="#3b82f6" stroke-width="1.5"'
                els.append(
                    f'<rect x="{xd:.1f}" y="{ry + 3}" width="{wd:.1f}" '
                    f'height="{_RH - 6}" rx="3" {fill_attr}>'
                    f'<title>{tip}</title></rect>'
                )
                if wd > 40:
                    lbl = d["name"][:12] + ("…" if len(d["name"]) > 12 else "")
                    els.append(
                        f'<text x="{xd + wd/2:.1f}" y="{ry + _RH/2 + 3.5}" '
                        f'text-anchor="middle" font-size="7.5" fill="#f0f6ff" '
                        f'font-family="-apple-system,sans-serif" pointer-events="none">'
                        f'{_esc(lbl)}</text>'
                    )
        y += len(rows) * _RH + 12

    # Track 2 — Secondary structure
    if ss_segs:
        els.append(_svg_label(y, "Sec. structure"))
        els.append(f'<rect x="{_LM}" y="{y}" width="{_TW}" height="{_TH}" rx="3" fill="#0d1424"/>')
        for s, e, t in ss_segs:
            xd = _svg_x(s, seq_len)
            wd = _svg_w(s, e, seq_len)
            els.append(
                f'<rect x="{xd:.1f}" y="{y}" width="{wd:.1f}" height="{_TH}" '
                f'fill="{_SS_COLORS[t]}"><title>{_SS_LABELS[t]} ({s}–{e})</title></rect>'
            )
        y += _TH + 10

    # Track 3 — pLDDT
    if plddt_segs:
        els.append(_svg_label(y, "pLDDT"))
        els.append(f'<rect x="{_LM}" y="{y}" width="{_TW}" height="{_TH}" rx="3" fill="#0d1424"/>')
        for s, e, col in plddt_segs:
            xd = _svg_x(s, seq_len)
            wd = _svg_w(s, e, seq_len)
            els.append(
                f'<rect x="{xd:.1f}" y="{y}" width="{wd:.1f}" height="{_TH}" fill="{col}">'
                f'<title>pLDDT {s}–{e}</title></rect>'
            )
        y += _TH + 10

    # Track 4 — Repeats
    rtx_entries  = [(s, e, f"RTX motif ({s}–{e})", "#a855f7") for s, e in rtx]
    rep_entries  = [(s, e, f"{cnt}× repeat (unit ~{unit} aa)", "#f59e0b")
                    for s, e, cnt, unit in repeats]
    all_rep = rtx_entries + rep_entries
    if all_rep:
        els.append(_svg_label(y, "Repeats"))
        els.append(f'<rect x="{_LM}" y="{y}" width="{_TW}" height="{_TH}" rx="3" fill="#0d1424"/>')
        for s, e, label, col in all_rep:
            xd = _svg_x(s, seq_len)
            wd = _svg_w(s, e, seq_len)
            els.append(
                f'<rect x="{xd:.1f}" y="{y}" width="{wd:.1f}" height="{_TH}" '
                f'rx="2" fill="{col}" opacity="0.85"><title>{_esc(label)}</title></rect>'
            )
            if wd > 35:
                short = label[:14] + ("…" if len(label) > 14 else "")
                els.append(
                    f'<text x="{xd + wd/2:.1f}" y="{y + _TH/2 + 3.5}" '
                    f'text-anchor="middle" font-size="7.5" fill="#f0f6ff" '
                    f'font-family="-apple-system,sans-serif" pointer-events="none">'
                    f'{_esc(short)}</text>'
                )
        y += _TH + 12

    # Ruler
    tick_step = max(1, round(seq_len / 15 / 10) * 10)
    ticks = list(range(1, seq_len + 1, tick_step))
    if seq_len not in ticks:
        ticks.append(seq_len)
    els.append(
        f'<line x1="{_LM}" y1="{y}" x2="{_LM + _TW}" y2="{y}" '
        f'stroke="#1e2d4a" stroke-width="1"/>'
    )
    for t in ticks:
        tx = _svg_x(t, seq_len)
        els.append(
            f'<line x1="{tx:.1f}" y1="{y}" x2="{tx:.1f}" y2="{y + 5}" '
            f'stroke="#1e2d4a" stroke-width="1"/>'
        )
        anchor = "start" if t == 1 else ("end" if t == seq_len else "middle")
        els.append(
            f'<text x="{tx:.1f}" y="{y + 15}" text-anchor="{anchor}" '
            f'font-size="9" fill="#334155" font-family="-apple-system,monospace">{t}</text>'
        )

    total_h = y + 25
    vb_w = _LM + _TW + _RM
    bg = (
        f'<rect width="{vb_w}" height="{total_h}" fill="#080c17"/>'
        f'<rect x="{_LM}" y="0" width="{_TW}" height="{total_h}" fill="#0a0e1a"/>'
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {vb_w} {total_h}" '
        f'width="100%" style="display:block;">{bg}{"".join(els)}</svg>'
    )


def _chip(text: str, col: str = "#334155") -> str:
    return (
        f'<span style="display:inline-flex;align-items:center;'
        f'background:{col}22;border:1px solid {col}55;border-radius:999px;'
        f'padding:3px 10px;font-size:11px;color:{col};font-weight:500;'
        f'font-family:-apple-system,sans-serif;">{_esc(text)}</span>'
    )


def show_domain_map(ipr_data, pdb_text, fasta_text) -> None:
    domains    = _parse_ipr_domains(ipr_data) if ipr_data else []
    plddt      = _parse_pdb_plddts(pdb_text)  if pdb_text  else []
    ss_segs    = _dssp_segments(pdb_text)      if pdb_text  else []
    seq        = _parse_fasta_seq(fasta_text)  if fasta_text else ""
    rtx        = _rtx_matches(seq)             if seq else []
    repeats    = _sliding_repeats(seq)         if seq else []
    seq_len    = (
        len(seq)
        or (max(d["end"] for d in domains) if domains else 0)
        or len(plddt)
    )

    if seq_len == 0:
        st.markdown('<p style="color:#334155;">No sequence data available.</p>',
                    unsafe_allow_html=True)
        return

    plddt_segs = _plddt_segments(plddt)
    svg        = _domain_map_svg(seq_len, domains, ss_segs, plddt_segs, repeats, rtx)

    st.markdown(
        '<div style="background:#080c17;border:1px solid #1e2d4a;border-radius:10px;'
        f'padding:18px 16px 10px;overflow-x:auto;">{svg}</div>',
        unsafe_allow_html=True,
    )

    # Legend
    legend_items = [
        ("PFAM", "#3b82f6"), ("PANTHER", "#22c55e"), ("PROSITE", "#a855f7"),
        ("GENE3D", "#0891b2"), ("SUPERFAMILY", None), ("Other", "#475569"),
    ]
    legend_html = " ".join(
        f'<span style="display:inline-flex;align-items:center;gap:4px;margin-right:6px;">'
        f'<span style="width:8px;height:8px;border-radius:2px;'
        f'{"background:" + col + ";" if col else "border:1px solid #3b82f6;"}'
        f'display:inline-block;"></span>'
        f'<span style="color:#475569;font-size:10.5px;">{lbl}</span></span>'
        for lbl, col in legend_items
    )
    st.markdown(
        f'<div style="display:flex;flex-wrap:wrap;gap:4px;margin-top:10px;">'
        f'{legend_html}</div>',
        unsafe_allow_html=True,
    )

    # Summary cards
    n_domains = len(set(f'{d["db"]}:{d["acc"]}' for d in domains))
    n_rep     = len(rtx) + sum(r[2] for r in repeats)
    dis_len   = (
        sum(e - s + 1 for s, e, t in ss_segs if t == "D") if ss_segs
        else sum(1 for s in plddt if s < 50)
    )
    pct_high  = round(sum(1 for s in plddt if s > 70) / len(plddt) * 100, 1) if plddt else 0

    st.markdown(_cards(
        ("Domains found",       n_domains,       None),
        ("Repeat units",        n_rep,           None),
        ("Disordered residues", dis_len,         None),
        ("% high confidence",   f"{pct_high}%",  None),
    ), unsafe_allow_html=True)

    # Feature chips
    chip_parts = []
    if n_domains:
        db_counts: dict = {}
        for d in domains:
            db_counts[d["db"]] = db_counts.get(d["db"], 0) + 1
        chip_parts.append(_chip(f'{n_domains} domain annotation(s)', "#3b82f6"))
        chip_parts.append(_chip(f'Most hits: {max(db_counts, key=db_counts.get)}', "#0891b2"))
    if rtx:
        chip_parts.append(_chip(f'{len(rtx)} RTX motif(s)', "#a855f7"))
    if repeats:
        chip_parts.append(_chip(f'{len(repeats)} tandem repeat region(s)', "#f59e0b"))
    if dis_len > 10:
        chip_parts.append(_chip(f'{dis_len} disordered residues', "#ef4444"))
    if pct_high >= 70:
        chip_parts.append(_chip(f'{pct_high}% high-confidence structure', "#22c55e"))
    elif 0 < pct_high <= 30:
        chip_parts.append(_chip(f'Low overall confidence ({pct_high}%)', "#f97316"))

    if chip_parts:
        st.markdown(
            '<div style="display:flex;flex-wrap:wrap;gap:6px;margin-top:14px;">'
            + "".join(chip_parts) + "</div>",
            unsafe_allow_html=True,
        )


# ── Tool status sidebar list ───────────────────────────────────────────────────

_ALL_TOOLS = ["InterProScan", "BLASTp", "Phobius", "HMMER", "FoldSeek"]


def _status_list(results: dict, active: list) -> str:
    rows = ""
    for tool in _ALL_TOOLS:
        is_active = tool in active
        if not is_active:
            dot   = "background:#0f1929"
            ncls  = ""
            badge = '<span class="dp-badge" style="color:#1e2d4a;">n/a</span>'
        elif tool in results:
            ok    = results[tool]["ok"]
            dot   = f"background:{'#22c55e' if ok else '#ef4444'}"
            ncls  = "" if ok else "err"
            lbl   = "done" if ok else "failed"
            lcol  = "#22c55e" if ok else "#ef4444"
            badge = f'<span class="dp-badge" style="color:{lcol};">{lbl}</span>'
        else:
            dot   = "background:#475569"
            ncls  = "active"
            badge = '<span class="dp-badge" style="color:#475569;">ready</span>'

        rows += (
            f'<div class="dp-status-row">'
            f'<span class="dp-dot" style="{dot};"></span>'
            f'<span class="dp-tool-name {ncls}">{_esc(tool)}</span>'
            f"{badge}</div>"
        )
    return rows


# ── Guard: EMAIL required ──────────────────────────────────────────────────────

if not EMAIL:
    st.markdown("<div style='height:28px'></div>", unsafe_allow_html=True)
    st.error(
        "**EBI_EMAIL not configured.** "
        "Add `EBI_EMAIL=your@email.com` to `.env` and restart, "
        "or set it as a Streamlit Cloud secret."
    )
    st.stop()

# ── Page HTML constants ─────────────────────────────────────────────────────────

_LANDING_HTML = """
<style>
.dpp-land {
  display: flex;
  align-items: center;
  min-height: 75vh;
  position: relative;
  overflow: hidden;
  border-radius: 10px;
  background: radial-gradient(ellipse at 32% 50%, #0d1e3a 0%, #06090f 70%);
}

/* ── Left: protein visualization (55%) ── */
.dpp-lvis {
  flex: 0 0 55%;
  position: relative;
  min-height: 75vh;
  overflow: hidden;
}
.dpp-lvis-bg {
  position: absolute; inset: 0;
  background:
    radial-gradient(ellipse at 38% 42%, rgba(6,182,212,0.22) 0%, transparent 52%),
    radial-gradient(ellipse at 62% 22%, rgba(139,92,246,0.18) 0%, transparent 46%),
    radial-gradient(ellipse at 28% 74%, rgba(244,114,182,0.14) 0%, transparent 42%),
    radial-gradient(ellipse at 72% 72%, rgba(16,185,129,0.12) 0%, transparent 40%);
}

/* ─── Helix ribbons (elongated, angled, glowing) ─── */
.dpp-hx { position: absolute; border-radius: 50%; pointer-events: none; }
.dpp-hx1 {
  width: 400px; height: 58px;
  top: 8%; left: 1%;
  background: linear-gradient(135deg,
    transparent 0%, rgba(6,182,212,0.5) 18%,
    rgba(34,211,238,1.0) 50%, rgba(6,182,212,0.5) 82%, transparent 100%);
  transform: rotate(-22deg);
  box-shadow: 0 0 38px rgba(34,211,238,0.90), 0 0 76px rgba(6,182,212,0.50), 0 0 130px rgba(6,182,212,0.22);
  animation: dpp-hx1a 13s ease-in-out infinite;
}
.dpp-hx2 {
  width: 310px; height: 50px;
  top: 43%; left: 18%;
  background: linear-gradient(135deg,
    transparent 0%, rgba(244,114,182,0.5) 18%,
    rgba(251,113,133,1.0) 50%, rgba(244,114,182,0.5) 82%, transparent 100%);
  transform: rotate(14deg);
  box-shadow: 0 0 32px rgba(251,113,133,0.90), 0 0 65px rgba(244,114,182,0.50), 0 0 110px rgba(244,114,182,0.22);
  animation: dpp-hx2a 16s ease-in-out infinite 2s;
}
.dpp-hx3 {
  width: 240px; height: 44px;
  top: 20%; left: 48%;
  background: linear-gradient(135deg,
    transparent 0%, rgba(52,211,153,0.5) 18%,
    rgba(110,231,183,1.0) 50%, rgba(52,211,153,0.5) 82%, transparent 100%);
  transform: rotate(-36deg);
  box-shadow: 0 0 28px rgba(110,231,183,0.88), 0 0 56px rgba(52,211,153,0.46), 0 0 95px rgba(52,211,153,0.20);
  animation: dpp-hx3a 19s ease-in-out infinite 4.5s;
}
.dpp-hx4 {
  width: 270px; height: 46px;
  top: 68%; left: 38%;
  background: linear-gradient(135deg,
    transparent 0%, rgba(245,158,11,0.5) 18%,
    rgba(251,191,36,1.0) 50%, rgba(245,158,11,0.5) 82%, transparent 100%);
  transform: rotate(9deg);
  box-shadow: 0 0 30px rgba(251,191,36,0.88), 0 0 60px rgba(245,158,11,0.46), 0 0 100px rgba(245,158,11,0.20);
  animation: dpp-hx4a 14s ease-in-out infinite 3.5s;
}

/* ─── Beta-sheet arrows ─── */
.dpp-sh { position: absolute; pointer-events: none; }
.dpp-sh1 {
  width: 210px; height: 46px;
  top: 52%; left: 3%;
  background: linear-gradient(90deg,
    rgba(139,92,246,0.92) 0%, rgba(167,139,250,1.00) 62%, rgba(139,92,246,0.30) 100%);
  clip-path: polygon(0 24%, 78% 24%, 78% 0%, 100% 50%, 78% 100%, 78% 76%, 0 76%);
  filter: drop-shadow(0 0 16px rgba(167,139,250,0.90)) drop-shadow(0 0 34px rgba(139,92,246,0.48));
  animation: dpp-sh1a 17s ease-in-out infinite 1.2s;
}
.dpp-sh2 {
  width: 168px; height: 38px;
  top: 30%; left: 3%;
  background: linear-gradient(90deg,
    rgba(239,68,68,0.88) 0%, rgba(252,165,165,1.00) 62%, rgba(239,68,68,0.28) 100%);
  clip-path: polygon(0 24%, 78% 24%, 78% 0%, 100% 50%, 78% 100%, 78% 76%, 0 76%);
  filter: drop-shadow(0 0 13px rgba(252,165,165,0.88)) drop-shadow(0 0 28px rgba(239,68,68,0.44));
  animation: dpp-sh2a 21s ease-in-out infinite 6s;
}

/* ─── Loop-arc connectors ─── */
.dpp-lp { position: absolute; border-radius: 50%; pointer-events: none; }
.dpp-lp1 {
  width: 115px; height: 76px;
  top: 33%; left: 30%;
  border-top: 2.5px solid rgba(34,211,238,0.58);
  animation: dpp-lp1a 13s ease-in-out infinite 1.8s;
}
.dpp-lp2 {
  width: 82px; height: 54px;
  top: 60%; left: 27%;
  border-top: 2px solid rgba(167,139,250,0.52);
  animation: dpp-lp2a 18s ease-in-out infinite 4s;
}

/* ─── Large ambient depth orbs ─── */
.dpp-orb { position: absolute; border-radius: 50%; pointer-events: none; }
.dpp-orb1 {
  width: 480px; height: 480px; top: -15%; left: -12%;
  background: radial-gradient(circle, rgba(6,182,212,0.14) 0%, transparent 62%);
  filter: blur(44px);
}
.dpp-orb2 {
  width: 400px; height: 400px; top: 28%; left: 16%;
  background: radial-gradient(circle, rgba(139,92,246,0.12) 0%, transparent 62%);
  filter: blur(38px);
}
.dpp-orb3 {
  width: 360px; height: 360px; top: 52%; left: 46%;
  background: radial-gradient(circle, rgba(244,114,182,0.10) 0%, transparent 62%);
  filter: blur(34px);
}

/* Edge fades into page background */
.dpp-rfade {
  position: absolute; top: 0; right: 0; bottom: 0; width: 190px;
  background: linear-gradient(to right, transparent, #06090f);
  pointer-events: none;
}
.dpp-tbfade {
  position: absolute; inset: 0;
  background: linear-gradient(to bottom, #06090f 0%, transparent 7%, transparent 93%, #06090f 100%);
  pointer-events: none;
}

/* ── Keyframes ── */
@keyframes dpp-hx1a{0%,100%{transform:rotate(-22deg) translate(0,0) scale(1)}35%{transform:rotate(-19deg) translate(10px,-7px) scale(1.03)}70%{transform:rotate(-25deg) translate(-7px,9px) scale(0.97)}}
@keyframes dpp-hx2a{0%,100%{transform:rotate(14deg) translate(0,0) scale(1)}40%{transform:rotate(16deg) translate(-11px,7px) scale(1.04)}80%{transform:rotate(11deg) translate(9px,-9px) scale(0.96)}}
@keyframes dpp-hx3a{0%,100%{transform:rotate(-36deg) translate(0,0) scale(1)}50%{transform:rotate(-33deg) translate(7px,11px) scale(1.05)}}
@keyframes dpp-hx4a{0%,100%{transform:rotate(9deg) translate(0,0) scale(1)}45%{transform:rotate(12deg) translate(-9px,-7px) scale(1.04)}}
@keyframes dpp-sh1a{0%,100%{transform:translate(0,0);opacity:.92}55%{transform:translate(13px,-9px);opacity:1}}
@keyframes dpp-sh2a{0%,100%{transform:translate(0,0);opacity:.88}65%{transform:translate(-11px,11px);opacity:.97}}
@keyframes dpp-lp1a{0%,100%{opacity:.46;transform:translate(0,0)}50%{opacity:.74;transform:translate(9px,-6px)}}
@keyframes dpp-lp2a{0%,100%{opacity:.40;transform:translate(0,0)}55%{opacity:.66;transform:translate(-7px,8px)}}
@keyframes dpp-pulse{0%,100%{opacity:1}50%{opacity:0.3}}

/* ── Right: text panel (45%) ── */
.dpp-ltxt {
  flex: 0 0 45%;
  display: flex;
  flex-direction: column;
  justify-content: center;
  padding: 0 52px 0 12px;
  position: relative;
  z-index: 2;
}
.dpp-lbr {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  font-size: 9.5px;
  font-weight: 700;
  letter-spacing: 0.20em;
  text-transform: uppercase;
  color: #22d3ee;
  margin-bottom: 24px;
  padding: 5px 14px;
  border: 1px solid rgba(34,211,238,0.22);
  border-radius: 20px;
  width: fit-content;
  background: rgba(34,211,238,0.05);
}
.dpp-lbr-dot {
  width: 6px; height: 6px; border-radius: 50%;
  background: #22d3ee;
  animation: dpp-pulse 2s infinite;
}
.dpp-lh {
  font-size: 46px !important; font-weight: 800 !important; line-height: 1.08 !important;
  color: #f1f5f9 !important; margin: 0 0 2px !important; letter-spacing: -0.03em !important;
}
.dpp-lhg {
  font-size: 46px !important; font-weight: 800 !important; line-height: 1.08 !important;
  background: linear-gradient(120deg, #22d3ee 0%, #818cf8 50%, #f472b6 100%) !important;
  -webkit-background-clip: text !important; -webkit-text-fill-color: transparent !important;
  background-clip: text !important; margin: 0 0 26px !important; letter-spacing: -0.03em !important;
}
.dpp-ldesc {
  font-size: 13.5px !important; line-height: 1.80 !important; color: #64748b !important;
  margin: 0 0 30px !important; max-width: 340px !important;
}
.dpp-ltools { display: flex !important; flex-wrap: wrap !important; gap: 6px !important; }
.dpp-ltool {
  background: rgba(15,23,42,0.9) !important; border: 1px solid rgba(51,65,85,0.60) !important;
  border-radius: 6px !important; padding: 4px 11px !important; font-size: 10.5px !important;
  color: #94a3b8 !important; font-weight: 500 !important; letter-spacing: 0.01em !important;
}
</style>

<div class="dpp-land">
  <div class="dpp-lvis">
    <div class="dpp-lvis-bg"></div>
    <div class="dpp-orb dpp-orb1"></div>
    <div class="dpp-orb dpp-orb2"></div>
    <div class="dpp-orb dpp-orb3"></div>
    <div class="dpp-hx dpp-hx1"></div>
    <div class="dpp-hx dpp-hx2"></div>
    <div class="dpp-hx dpp-hx3"></div>
    <div class="dpp-hx dpp-hx4"></div>
    <div class="dpp-sh dpp-sh1"></div>
    <div class="dpp-sh dpp-sh2"></div>
    <div class="dpp-lp dpp-lp1"></div>
    <div class="dpp-lp dpp-lp2"></div>
    <div class="dpp-rfade"></div>
    <div class="dpp-tbfade"></div>
  </div>
  <div class="dpp-ltxt">
    <span class="dpp-lbr"><span class="dpp-lbr-dot"></span>Dark Proteome Pipeline</span>
    <p class="dpp-lh">Illuminate the</p>
    <p class="dpp-lhg">dark proteome</p>
    <p class="dpp-ldesc">Automated structural and functional annotation for
uncharacterised microbial proteins — five complementary EBI REST APIs
in one pipeline.</p>
    <div class="dpp-ltools">
      <span class="dpp-ltool">InterProScan</span>
      <span class="dpp-ltool">BLASTp</span>
      <span class="dpp-ltool">FoldSeek</span>
      <span class="dpp-ltool">Phobius</span>
      <span class="dpp-ltool">HMMER</span>
    </div>
  </div>
</div>
"""

# ── Landing hero: 3Dmol viewer + Anime.js text (template; __VASE_B64__ filled at runtime) ──

_LANDING_HERO_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
*{margin:0;padding:0;box-sizing:border-box;}
html,body{background:#0a0e1a;width:100%;height:100%;overflow:hidden;
  font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','Inter',sans-serif;}
.hero{display:flex;width:100%;height:100%;}

/* ── Left: 3Dmol protein viewer ── */
.hero-left{flex:0 0 52%;position:relative;background:#0a0e1a;overflow:hidden;}
#mol-viewer{width:100%;height:100%;position:relative;}
.vload{
  position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
  background:#0a0e1a;z-index:10;transition:opacity 1s ease;
}
.vload-ring{
  width:36px;height:36px;border-radius:50%;
  border:2px solid rgba(34,211,238,0.15);
  border-top-color:#22d3ee;
  animation:spin-ring 1s linear infinite;
}
@keyframes spin-ring{to{transform:rotate(360deg)}}

/* ── Right: futuristic grid panel ── */
.hero-right{
  flex:0 0 48%;
  display:flex;flex-direction:column;justify-content:center;
  padding:0 56px 0 52px;
  background-color:#020408;
  background-image:
    linear-gradient(rgba(34,211,238,.032) 1px,transparent 1px),
    linear-gradient(90deg,rgba(34,211,238,.032) 1px,transparent 1px);
  background-size:44px 44px;
  position:relative;overflow:hidden;
}
.hero-right::before{
  content:'';position:absolute;left:0;top:0;bottom:0;width:1px;
  background:linear-gradient(to bottom,transparent 5%,rgba(34,211,238,.28) 50%,transparent 95%);
}
/* faint corner radiance */
.hero-right::after{
  content:'';position:absolute;right:-80px;top:-80px;
  width:320px;height:320px;border-radius:50%;
  background:radial-gradient(circle,rgba(139,92,246,.06) 0%,transparent 65%);
  pointer-events:none;
}

.badge{
  display:inline-flex;align-items:center;gap:8px;
  font-size:9px;font-weight:700;letter-spacing:.22em;text-transform:uppercase;
  color:#22d3ee;margin-bottom:28px;padding:5px 14px;
  border:1px solid rgba(34,211,238,.18);border-radius:20px;
  background:rgba(34,211,238,.04);width:fit-content;
}
.badge-dot{width:5px;height:5px;border-radius:50%;background:#22d3ee;animation:pulse-dot 2s infinite;}
@keyframes pulse-dot{0%,100%{opacity:1}50%{opacity:.2}}

.title-1{font-size:44px;font-weight:800;line-height:1.08;color:#f1f5f9;margin:0 0 3px;letter-spacing:-.03em;}
.title-2{
  font-size:44px;font-weight:800;line-height:1.08;
  background:linear-gradient(120deg,#22d3ee 0%,#818cf8 52%,#f472b6 100%);
  -webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;
  margin:0 0 28px;letter-spacing:-.03em;
}
.desc{font-size:13px;line-height:1.85;color:#475569;margin:0 0 32px;max-width:360px;}
.chips{display:flex;flex-wrap:wrap;gap:6px;}
.chip{
  background:rgba(2,4,8,.95);border:1px solid rgba(34,211,238,.15);
  border-radius:6px;padding:4px 12px;font-size:10px;color:#64748b;
  font-weight:500;letter-spacing:.02em;
}
.dpp-char{display:inline-block;will-change:transform,opacity;}
</style>
</head>
<body>
<div class="hero">

  <!-- Left: 3Dmol auto-rotating rainbow ribbon -->
  <div class="hero-left">
    <div id="mol-viewer"></div>
    <div class="vload" id="vload"><div class="vload-ring"></div></div>
  </div>

  <!-- Right: grid text panel -->
  <div class="hero-right">
    <span class="badge" id="badge"><span class="badge-dot"></span>Dark Proteome Pipeline</span>
    <p class="title-1" id="t1">Illuminate the</p>
    <p class="title-2" id="t2">dark proteome</p>
    <p class="desc"    id="desc">Automated structural and functional annotation for uncharacterised microbial proteins — five complementary EBI REST APIs in one pipeline.</p>
    <div class="chips" id="chips">
      <span class="chip">InterProScan</span>
      <span class="chip">BLASTp</span>
      <span class="chip">FoldSeek</span>
      <span class="chip">Phobius</span>
      <span class="chip">HMMER</span>
    </div>
  </div>
</div>

<!-- Pre-hide text; fallback timer reveals everything if CDNs fail -->
<script>
var _fb=setTimeout(function(){
  ['#badge','#t1','#t2','#desc'].forEach(function(s){
    var e=document.querySelector(s);if(e){e.style.opacity='1';e.style.transform='none';}
  });
  document.querySelectorAll('.chip,.dpp-char').forEach(function(e){e.style.opacity='1';e.style.transform='none';});
  var vl=document.getElementById('vload');if(vl)vl.style.opacity='0';
},5000);
['#badge','#t2','#desc'].forEach(function(s){
  var e=document.querySelector(s);if(e)e.style.opacity='0';
});
document.querySelectorAll('.chip').forEach(function(e){e.style.opacity='0';e.style.transform='translateY(6px)';});
</script>

<!-- 3Dmol.js protein viewer (VASE.pdb embedded as base64) -->
<script>
(function(){
  var b64='__VASE_B64__';
  function initViewer(){
    if(typeof $3Dmol==='undefined'){setTimeout(initViewer,80);return;}
    var container=document.getElementById('mol-viewer');
    var viewer=$3Dmol.createViewer(container,{
      backgroundColor:'#0a0e1a',
      antialias:true,
      disableFog:true,
    });
    viewer.addModel(atob(b64),'pdb');
    viewer.setStyle({},{cartoon:{colorscheme:'rainbow'}});
    viewer.zoomTo();
    viewer.zoom(0.88);
    viewer.spin('y',0.8);
    viewer.render();
    var vl=document.getElementById('vload');
    if(vl){vl.style.opacity='0';setTimeout(function(){vl.style.display='none';},1100);}
  }
  var s=document.createElement('script');
  s.src='https://3dmol.org/build/3Dmol-min.js';
  s.onload=function(){setTimeout(initViewer,60);};
  s.onerror=function(){var vl=document.getElementById('vload');if(vl)vl.style.opacity='0';};
  document.head.appendChild(s);
})();
</script>

<!-- Anime.js text animations for right panel -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/animejs/3.2.1/anime.min.js"
        onload="dppAnim()"
        onerror="clearTimeout(_fb);['#badge','#t1','#t2','#desc'].forEach(function(s){var e=document.querySelector(s);if(e){e.style.opacity='1';e.style.transform='none';}});document.querySelectorAll('.chip').forEach(function(e){e.style.opacity='1';e.style.transform='none';})">
</script>
<script>
function dppAnim(){
  clearTimeout(_fb);
  var t1=document.getElementById('t1');
  t1.innerHTML=t1.textContent.split('').map(function(c){
    return '<span class="dpp-char" style="opacity:0;transform:translateY(16px)">'+(c===' '?'&nbsp;':c)+'</span>';
  }).join('');
  document.getElementById('t2').style.transform='translateY(18px)';
  document.getElementById('desc').style.transform='translateY(14px)';
  document.getElementById('badge').style.transform='translateY(8px)';
  var tl=anime.timeline({easing:'easeOutExpo'});
  tl.add({targets:'#badge',opacity:[0,1],translateY:[8,0],duration:600});
  tl.add({targets:'#t1 .dpp-char',opacity:[0,1],translateY:[16,0],duration:550,delay:anime.stagger(28)},'-=300');
  tl.add({targets:'#t2',opacity:[0,1],translateY:[18,0],duration:700},'-=150');
  tl.add({targets:'#desc',opacity:[0,1],translateY:[14,0],duration:550},'-=480');
  tl.add({targets:'#chips .chip',opacity:[0,1],translateY:[6,0],duration:360,delay:anime.stagger(65)},'-=260');
}
</script>
</body>
</html>"""


@st.cache_resource
def _build_landing_html():
    """Read VASE.pdb once, base64-encode it, and inject into the hero template."""
    try:
        b64 = base64.b64encode(
            pathlib.Path("data/input/VASE.pdb").read_bytes()
        ).decode()
    except (FileNotFoundError, OSError):
        b64 = ""
    return _LANDING_HERO_TEMPLATE.replace("__VASE_B64__", b64)


# ── Page routing ───────────────────────────────────────────────────────────────

_page = st.session_state.get("page", "landing")

if _page == "landing":
    components.html(_build_landing_html(), height=620, scrolling=False)
    st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)
    _, _cta, _ = st.columns([2, 3, 2])
    with _cta:
        if st.button("Begin Analysis  →", type="primary", use_container_width=True):
            st.session_state["page"] = "analysis"
            st.rerun()
    st.stop()

# ── Analysis page ──────────────────────────────────────────────────────────────

col_left, col_right = st.columns([3, 9], gap="small")

# ── Left panel — input + status ────────────────────────────────────────────────

with col_left:
    if st.button("← Home", key="back_home"):
        st.session_state["page"] = "landing"
        st.rerun()
    st.markdown(
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;'
        'letter-spacing:0.12em;text-transform:uppercase;margin:18px 0 16px;">Input</p>',
        unsafe_allow_html=True,
    )

    fasta_file = st.file_uploader(
        "FASTA FILE",
        type=["fasta", "fa", "txt"],
        help="InterProScan · BLASTp · Phobius · HMMER",
    )
    st.markdown("<div style='height:6px'></div>", unsafe_allow_html=True)
    pdb_file = st.file_uploader(
        "PDB FILE",
        type=["pdb"],
        help="FoldSeek structural search (optional)",
    )
    st.markdown("<div style='height:2px'></div>", unsafe_allow_html=True)

    run = st.button(
        "Run Analysis",
        type="primary",
        disabled=not (fasta_file or pdb_file),
    )

    if run:
        try:
            fasta_text = fasta_file.read().decode() if fasta_file else None
            pdb_text   = pdb_file.read().decode()   if pdb_file   else None

            # ── Cache check: same files → load instantly without re-running ──
            _cached = _load_cache(fasta_text, pdb_text)
            if _cached:
                st.success(
                    "✅ Results restored from cache — same files were run previously. "
                    "To force a fresh run, clear the upload and re-upload your files."
                )
                st.session_state["results"]      = _cached["results"]
                st.session_state["fasta_text"]   = _cached["fasta_text"]
                st.session_state["pdb_text"]     = _cached["pdb_text"]
                st.session_state["active_tools"] = _cached["active_tools"]
                st.session_state["protein_name"] = _extract_protein_name(fasta_text)
                st.rerun()

            tasks: dict = {}
            if fasta_text:
                tasks["InterProScan"] = (run_interproscan, fasta_text)
                tasks["BLASTp"]       = (run_blast,        fasta_text)
                tasks["Phobius"]      = (run_phobius,       fasta_text)
                tasks["HMMER"]        = (run_hmmer,         fasta_text)
            if pdb_text:
                tasks["FoldSeek"]     = (run_foldseek,      pdb_text)

            n_tasks = len(tasks)
            _protein_name = _extract_protein_name(fasta_text)
            st.session_state["active_tools"] = list(tasks.keys())
            st.session_state["protein_name"] = _protein_name
            # Seed session state immediately so partial results appear on refresh
            st.session_state["results"]    = {}
            st.session_state["fasta_text"] = fasta_text
            st.session_state["pdb_text"]   = pdb_text

            results: dict = {}
            _pipeline_start = time.time()
            _start_str = time.strftime("%H:%M:%S")
            print(f"[PIPELINE] Starting {n_tasks} tools sequentially: {list(tasks.keys())}", flush=True)

            # Run sequentially — avoids OOM and thread contention on free-tier cloud.
            # Each tool's status is streamed live to the UI via st.status().
            with st.status(
                f"Running {n_tasks} tool{'s' if n_tasks != 1 else ''} — this takes 8–15 minutes…",
                expanded=True,
            ) as _status_box:
                st.write(f"🕐 Pipeline started at **{_start_str}** — keep this tab open")
                for idx, (name, (fn, arg)) in enumerate(tasks.items(), 1):
                    _pipeline_elapsed = int(time.time() - _pipeline_start)
                    st.write(
                        f"⏳  [{idx}/{n_tasks}]  **{name}** — running…"
                        f"  _(+{_pipeline_elapsed}s  |  {time.strftime('%H:%M:%S')})_"
                    )
                    if name == "BLASTp":
                        st.info(
                            "BLASTp can take 10–15 minutes from cloud servers — this is normal."
                        )
                    print(f"[PIPELINE] Starting {name}", flush=True)
                    _t0 = time.time()
                    try:
                        data = fn(arg)
                        results[name] = {"ok": True, "data": data}
                        _elapsed = round(time.time() - _t0, 1)
                        _total_so_far = int(time.time() - _pipeline_start)
                        st.write(
                            f"✅  [{idx}/{n_tasks}]  **{name}** — done in {_elapsed}s"
                            f"  _(pipeline total: {_total_so_far}s)_"
                        )
                        print(f"[PIPELINE] {name} done in {_elapsed}s", flush=True)
                    except Exception as _exc:
                        results[name] = {"ok": False, "error": str(_exc)}
                        _elapsed = round(time.time() - _t0, 1)
                        st.write(f"❌  [{idx}/{n_tasks}]  **{name}** — failed: {_exc}")
                        print(f"[PIPELINE] {name} FAILED in {_elapsed}s — {_exc}", flush=True)
                    # Write partial results to session_state after every tool so that
                    # a browser refresh mid-pipeline shows whatever completed so far.
                    st.session_state["results"] = dict(results)

                _total   = round(time.time() - _pipeline_start, 1)
                _n_ok    = sum(1 for r in results.values() if r["ok"])
                _all_ok  = _n_ok == n_tasks
                _status_box.update(
                    label=(f"{'All' if _all_ok else f'{_n_ok}/{n_tasks}'} tools completed"
                           f" — {_total}s total (started {_start_str})"),
                    state="complete" if _all_ok else "error",
                )
                print(f"[PIPELINE] All done in {_total}s — "
                      f"ok={[k for k,v in results.items() if v['ok']]} "
                      f"failed={[k for k,v in results.items() if not v['ok']]}", flush=True)

            st.session_state["results"]    = results
            st.session_state["fasta_text"] = fasta_text
            st.session_state["pdb_text"]   = pdb_text

            # Save to disk — same files will reload instantly even after a server restart.
            _save_cache(fasta_text, pdb_text, results, list(tasks.keys()))

        except Exception as _pipeline_crash:
            import traceback as _tb
            st.error(f"Pipeline crashed unexpectedly: {_pipeline_crash}")
            st.code(_tb.format_exc())
            print(f"[PIPELINE] CRASH: {_pipeline_crash}", flush=True)

    # Status list — reads session_state so it updates after each run
    _res    = st.session_state.get("results",      {})
    _active = st.session_state.get("active_tools", [])

    st.markdown(
        '<hr style="border:none;border-top:1px solid #0d1424;margin:22px 0 16px;">'
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;'
        'letter-spacing:0.12em;text-transform:uppercase;margin-bottom:6px;">Tool Status</p>'
        + _status_list(_res, _active),
        unsafe_allow_html=True,
    )

# ── Analysis empty-state HTML ───────────────────────────────────────────────────

_EMPTY_STATE_HTML = """
<div style="display:flex;flex-direction:column;align-items:center;
            justify-content:center;min-height:440px;">
  <div style="opacity:0.18;margin-bottom:22px;">
    <svg width="70" height="70" viewBox="0 0 70 70" fill="none">
      <circle cx="35" cy="35" r="32" stroke="#94a3b8" stroke-width="1.5"/>
      <path d="M18 35 Q18 20 35 20 Q52 20 52 35" stroke="#22d3ee" stroke-width="2.2"
            fill="none" stroke-linecap="round"/>
      <path d="M52 35 Q52 50 35 50 Q18 50 18 35" stroke="#a855f7" stroke-width="2.2"
            fill="none" stroke-linecap="round"/>
      <circle cx="35" cy="35" r="3.5" fill="#475569"/>
    </svg>
  </div>
  <p style="color:#334155;font-size:13px;text-align:center;
            line-height:1.80;margin:0;max-width:220px;">
    Upload a FASTA file in the left panel<br>then click
    <span style="color:#3b82f6;font-weight:600;">Run Analysis</span>
  </p>
</div>
"""

# ── Right panel — results ──────────────────────────────────────────────────────

with col_right:
    _res    = st.session_state.get("results",      {})
    _active = st.session_state.get("active_tools", [])

    _pdb_ss   = st.session_state.get("pdb_text")
    _fasta_ss = st.session_state.get("fasta_text")

    if not _res and not _pdb_ss:
        _had_prior = bool(st.session_state.get("active_tools"))
        if _had_prior:
            st.info(
                "**Results may have expired** — your session was reset. "
                "Re-upload your files and click **Run Analysis** to restore results from cache "
                "(if the same files were run before, they load instantly without re-running)."
            )
        st.markdown(_EMPTY_STATE_HTML, unsafe_allow_html=True)
    else:
        if _res:
            n_ok    = sum(r["ok"] for r in _res.values())
            n_total = len(_res)
            all_ok  = n_ok == n_total
            bar_col = "#22c55e" if all_ok else "#f59e0b"
            bar_bg  = "rgba(34,197,94,0.06)"  if all_ok else "rgba(245,158,11,0.06)"
            bar_bdr = "rgba(34,197,94,0.15)"  if all_ok else "rgba(245,158,11,0.15)"

            _sum_col, _dl_col, _clr_col = st.columns([2.4, 0.9, 0.7])
            with _sum_col:
                st.markdown(
                    f'<div style="background:{bar_bg};border:1px solid {bar_bdr};border-radius:8px;'
                    f'padding:10px 16px;display:flex;align-items:center;gap:8px;">'
                    f'<span style="color:{bar_col};font-size:12.5px;font-weight:600;">'
                    f'&#x2714;&ensp;{n_ok}/{n_total} tools completed</span></div>',
                    unsafe_allow_html=True,
                )
            with _dl_col:
                _pname = st.session_state.get("protein_name", "protein")
                _fname = f"dpp_{_pname}_{time.strftime('%Y%m%d_%H%M%S')}.json"
                st.download_button(
                    label="⬇ JSON",
                    data=_results_to_json_bytes(_res, _pname),
                    file_name=_fname,
                    mime="application/json",
                    help="Save all results locally before the session expires",
                    use_container_width=True,
                )
            with _clr_col:
                if st.button("✕ Clear", help="Clear results and return to home",
                             use_container_width=True, key="clear_results"):
                    for _k in ["results", "fasta_text", "pdb_text",
                               "active_tools", "protein_name"]:
                        st.session_state.pop(_k, None)
                    st.session_state["page"] = "landing"
                    st.rerun()
            st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)

            # Inline error banners for failed tools — full error text, always visible
            for name in _active:
                r = _res.get(name)
                if r and not r["ok"]:
                    st.error(f"**{name} failed** — {r['error']}")

        # Build tab list: Structure → Domain Map → tool results
        finished        = [t for t in _active if _res.get(t, {}).get("ok")]
        _has_ipr        = _res.get("InterProScan", {}).get("ok")
        _has_domain_map = bool(_pdb_ss or _has_ipr)

        tab_names = (
            (["🧬 Structure"]   if _pdb_ss          else []) +
            (["🗺️ Domain Map"]  if _has_domain_map  else []) +
            finished
        )

        if tab_names:
            tabs   = st.tabs(tab_names)
            offset = 0
            if _pdb_ss:
                with tabs[offset]:
                    _phobius_data = (
                        _res["Phobius"]["data"]
                        if _res.get("Phobius", {}).get("ok") else None
                    )
                    show_structure(_pdb_ss, _fasta_ss, _phobius_data)
                offset += 1
            if _has_domain_map:
                with tabs[offset]:
                    _ipr_data = _res["InterProScan"]["data"] if _has_ipr else None
                    show_domain_map(_ipr_data, _pdb_ss, _fasta_ss)
                offset += 1
            for i, name in enumerate(finished):
                with tabs[offset + i]:
                    data = _res[name]["data"]
                    try:
                        if   name == "InterProScan": show_interproscan(data)
                        elif name == "BLASTp":       show_blast(data)
                        elif name == "Phobius":      show_phobius(data)
                        elif name == "HMMER":        show_hmmer(data)
                        elif name == "FoldSeek":     show_foldseek(data)
                    except Exception as _tab_exc:
                        st.error(f"**{name} display error:** {_tab_exc}")
                        with st.expander("Raw data (for debugging)"):
                            try:
                                st.code(json.dumps(data, indent=2)[:4000]
                                        if isinstance(data, (dict, list)) else str(data)[:4000])
                            except Exception:
                                st.code(str(data)[:4000])

# ── Debug panel (only visible when ?debug=1 is in the URL) ───────────────────

_qp = st.query_params
if _qp.get("debug") == "1":
    st.markdown(
        '<hr style="border:none;border-top:1px solid #1e2d4a;margin:28px 0 12px;">'
        '<p style="color:#f59e0b;font-size:9.5px;font-weight:700;letter-spacing:.1em;'
        'text-transform:uppercase;margin:0 0 12px;">🛠 Debug Panel</p>',
        unsafe_allow_html=True,
    )
    _dcol1, _dcol2 = st.columns(2)

    with _dcol1:
        if st.button("Test BLAST connection", key="dbg_blast"):
            _test_seq = ">test\nMKTLLLTLVV"
            with st.spinner("Submitting 10-aa test sequence to NCBI BLAST…"):
                try:
                    _r = requests.post(
                        "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                        data={
                            "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
                            "ENTREZ_QUERY": "bacteria[organism]",
                            "QUERY": _test_seq, "FORMAT_TYPE": "JSON2",
                            "HITLIST_SIZE": 5, "MATRIX_NAME": "BLOSUM62",
                            "EXPECT": "0.001",
                            "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
                        },
                        timeout=60,
                    )
                    _rid_dbg = None
                    for _ln in _r.text.splitlines():
                        if _ln.startswith("    RID = "):
                            _rid_dbg = _ln.split("=", 1)[1].strip()
                    st.markdown(
                        f'<div style="background:rgba(34,197,94,0.08);border:1px solid '
                        f'rgba(34,197,94,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#94a3b8;">'
                        f'<strong style="color:#22c55e;">HTTP {_r.status_code}</strong><br>'
                        f'RID: <code style="color:#f59e0b;">{_rid_dbg or "NOT FOUND in response"}</code><br>'
                        f'Response length: {len(_r.text)} chars</div>',
                        unsafe_allow_html=True,
                    )
                    if _rid_dbg:
                        st.session_state["_dbg_blast_rid"] = _rid_dbg
                except Exception as _e:
                    st.markdown(
                        f'<div style="background:rgba(239,68,68,0.08);border:1px solid '
                        f'rgba(239,68,68,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#ef4444;">Submission error:<br>'
                        f'<code>{_html.escape(str(_e))}</code></div>',
                        unsafe_allow_html=True,
                    )

        # If we have a RID from a previous test, show a poll-once button
        _dbg_rid = st.session_state.get("_dbg_blast_rid")
        if _dbg_rid:
            if st.button(f"Poll status of RID {_dbg_rid}", key="dbg_blast_poll"):
                try:
                    _rp = requests.get(
                        "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                        params={
                            "CMD": "Get", "RID": _dbg_rid,
                            "FORMAT_OBJECT": "SearchInfo",
                            "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
                        },
                        timeout=30,
                    )
                    _poll_status = "UNKNOWN"
                    _hits_line = ""
                    for _ln in _rp.text.splitlines():
                        if "Status=" in _ln:
                            _poll_status = _ln.strip().split("=", 1)[1].strip()
                        if "ThereAreHits=" in _ln:
                            _hits_line = _ln.strip()
                    st.markdown(
                        f'<div style="background:rgba(59,130,246,0.08);border:1px solid '
                        f'rgba(59,130,246,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#94a3b8;">'
                        f'Poll status: <strong style="color:#3b82f6;">{_poll_status}</strong><br>'
                        f'{_html.escape(_hits_line) if _hits_line else ""}</div>',
                        unsafe_allow_html=True,
                    )
                except Exception as _e:
                    st.error(f"Poll error: {_e}")

    with _dcol2:
        if st.button("Test HMMER connection", key="dbg_hmmer"):
            with st.spinner("Sending HEAD request to EBI HMMER endpoint…"):
                _hmmer_url = "https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan/run"
                try:
                    _rh = requests.head(_hmmer_url, timeout=20)
                    st.markdown(
                        f'<div style="background:rgba(34,197,94,0.08);border:1px solid '
                        f'rgba(34,197,94,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#94a3b8;">'
                        f'<strong style="color:#22c55e;">HTTP {_rh.status_code}</strong><br>'
                        f'Server reachable ✓</div>',
                        unsafe_allow_html=True,
                    )
                except Exception as _e:
                    st.markdown(
                        f'<div style="background:rgba(239,68,68,0.08);border:1px solid '
                        f'rgba(239,68,68,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#ef4444;">HEAD request error:<br>'
                        f'<code>{_html.escape(str(_e))}</code></div>',
                        unsafe_allow_html=True,
                    )

            # Also submit a real minimal HMMER job and report the job ID
            st.markdown('<p style="color:#64748b;font-size:11px;margin:8px 0 4px;">Submitting minimal HMMER job…</p>', unsafe_allow_html=True)
            _test_seq_h = ">test\nMKTLLLTLVVVTIACSFA"
            try:
                _rj = requests.post(
                    _hmmer_url,
                    data={"email": EMAIL, "sequence": _test_seq_h, "database": "pfam", "E": "1.0"},
                    timeout=60,
                )
                _jid_dbg = _rj.text.strip()
                st.markdown(
                    f'<div style="background:rgba(59,130,246,0.08);border:1px solid '
                    f'rgba(59,130,246,0.2);border-radius:6px;padding:10px 14px;'
                    f'font-size:12px;color:#94a3b8;">'
                    f'Submission HTTP {_rj.status_code}<br>'
                    f'Job ID: <code style="color:#f59e0b;">{_html.escape(_jid_dbg[:80])}</code></div>',
                    unsafe_allow_html=True,
                )
                if _jid_dbg:
                    st.session_state["_dbg_hmmer_jid"] = _jid_dbg
            except Exception as _e:
                st.markdown(
                    f'<div style="background:rgba(239,68,68,0.08);border:1px solid '
                    f'rgba(239,68,68,0.2);border-radius:6px;padding:10px 14px;'
                    f'font-size:12px;color:#ef4444;">HMMER submit error:<br>'
                    f'<code>{_html.escape(str(_e))}</code></div>',
                    unsafe_allow_html=True,
                )

        # Poll the HMMER job once if we have a stored job ID
        _dbg_jid = st.session_state.get("_dbg_hmmer_jid")
        if _dbg_jid:
            if st.button(f"Poll HMMER status", key="dbg_hmmer_poll"):
                try:
                    _rps = requests.get(
                        f"https://www.ebi.ac.uk/Tools/services/rest/hmmer3_hmmscan/status/{_dbg_jid}",
                        timeout=20,
                    )
                    st.markdown(
                        f'<div style="background:rgba(59,130,246,0.08);border:1px solid '
                        f'rgba(59,130,246,0.2);border-radius:6px;padding:10px 14px;'
                        f'font-size:12px;color:#94a3b8;">'
                        f'HTTP {_rps.status_code} — Status: '
                        f'<strong style="color:#3b82f6;">{_html.escape(_rps.text.strip())}</strong></div>',
                        unsafe_allow_html=True,
                    )
                except Exception as _e:
                    st.error(f"Poll error: {_e}")

    st.markdown(
        f'<p style="color:#475569;font-size:10.5px;margin-top:10px;">'
        f'EMAIL env: {"set ✓" if EMAIL else "NOT SET ✗"} &nbsp;·&nbsp; '
        f'Python sees outbound HTTPS: check buttons above</p>',
        unsafe_allow_html=True,
    )

# ── Status bar ─────────────────────────────────────────────────────────────────

_res_sb = st.session_state.get("results", {})
if _res_sb:
    _ok  = sum(r["ok"] for r in _res_sb.values())
    _tot = len(_res_sb)
    _sc  = "#22c55e" if _ok == _tot else "#f59e0b"
    _st  = f"&#x25CF;&ensp;{_ok}/{_tot} tools completed"
else:
    _sc = "#334155"
    _st = "&#x25CB;&ensp;Ready"

st.markdown(f"""
<div style="position:fixed;bottom:0;left:0;right:0;background:#080c17;
            border-top:1px solid #1e2d4a;padding:7px 28px;
            display:flex;justify-content:space-between;align-items:center;
            z-index:100;font-size:11px;font-family:-apple-system,sans-serif;">
  <span style="color:{_sc};font-weight:600;">{_st}</span>
  <span style="color:#1e2d4a;">darkproteome.streamlit.app</span>
</div>
""", unsafe_allow_html=True)
