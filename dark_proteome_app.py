"""
DarkProteome — Professional scientific annotation tool for the microbial dark proteome.
Run:  streamlit run dark_proteome_app.py
"""

import base64
import html as _html
import io
import json
import math
import os
import tarfile
import tempfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from dotenv import load_dotenv

load_dotenv()
EMAIL  = os.environ.get("EBI_EMAIL", "")
GITHUB = "https://github.com/Akjaer05/dark-proteome-pipeline"

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

/* ── Column layout ──────────────────────────────────────── */
[data-testid="stHorizontalBlock"] {
    gap: 0 !important;
    align-items: stretch !important;
}
[data-testid="stHorizontalBlock"] > [data-testid="column"]:first-child {
    background: #080c17 !important;
    border-right: 1px solid #1e2d4a !important;
    padding: 28px 20px 60px !important;
    min-height: 65vh;
}
[data-testid="stHorizontalBlock"] > [data-testid="column"]:last-child {
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
    box-shadow: 0 6px 22px rgba(59, 130, 246, 0.45) !important;
    transform: translateY(-1px) !important;
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
    <span style="background:rgba(59,130,246,0.1); color:#3b82f6; font-size:10px;
                 font-weight:700; padding:2px 8px; border-radius:20px;
                 border:1px solid rgba(59,130,246,0.2); letter-spacing:0.06em;">v1.0</span>
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

# ── Hero ───────────────────────────────────────────────────────────────────────

def _badge(text, coming_soon=False):
    if coming_soon:
        return (
            f'<span style="background:rgba(100,116,139,0.05);color:#334155;'
            f'border:1px solid #1e2d4a;padding:3px 11px;border-radius:20px;'
            f'font-size:11px;font-weight:500;">'
            f'{text}&thinsp;<span style="font-size:9px;opacity:0.65;'
            f'font-style:italic;">soon</span></span>'
        )
    return (
        f'<span style="background:rgba(34,197,94,0.07);color:#22c55e;'
        f'border:1px solid rgba(34,197,94,0.18);padding:3px 11px;'
        f'border-radius:20px;font-size:11px;font-weight:500;">{text}</span>'
    )

_tool_badges = " ".join([
    _badge("InterProScan"),
    _badge("BLASTp"),
    _badge("Phobius"),
    _badge("HMMER"),
    _badge("FoldSeek"),
    _badge("SignalP 6.0", coming_soon=True),
    _badge("HHpred",      coming_soon=True),
])

st.markdown(f"""
<div style="background:linear-gradient(180deg,#0c1527 0%,#0a0e1a 100%);
            padding:44px 32px 38px; border-bottom:1px solid #1e2d4a;">
  <div style="max-width:820px;">
    <div style="display:inline-flex; align-items:center; gap:5px;
                background:rgba(59,130,246,0.07); color:#3b82f6;
                border:1px solid rgba(59,130,246,0.16);
                padding:3px 11px; border-radius:20px;
                font-size:9.5px; font-weight:700; letter-spacing:0.12em;
                text-transform:uppercase; margin-bottom:18px;">
      &#x25CF;&ensp;Structural Bioinformatics
    </div>
    <h1 style="font-size:27px; font-weight:700; color:#f0f6ff; line-height:1.32;
               letter-spacing:-0.02em; margin:0 0 12px;">
      Annotate hypothetical proteins<br>from the microbial dark proteome
    </h1>
    <p style="font-size:13.5px; color:#475569; line-height:1.65;
              margin:0 0 22px; max-width:600px;">
      Multi-tool parallel annotation pipeline for uncharacterised bacterial proteins.
      Upload a FASTA sequence and optionally an AlphaFold2 PDB structure to run all
      tools simultaneously.
    </p>
    <div style="display:flex; flex-wrap:wrap; gap:7px; align-items:center;">
      <span style="font-size:9.5px; color:#334155; font-weight:700;
                   letter-spacing:0.1em; text-transform:uppercase; margin-right:2px;">
        Tools
      </span>
      {_tool_badges}
    </div>
  </div>
</div>
""", unsafe_allow_html=True)

# ── Tool runners (unchanged) ───────────────────────────────────────────────────

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


def run_blast(sequence):
    url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    try:
        r = requests.put(url, params={
            "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "nr",
            "QUERY": sequence, "FORMAT_TYPE": "JSON2",
            "EMAIL": EMAIL, "TOOL": "dark-proteome-pipeline",
        }, timeout=60)
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
            }, timeout=60)
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
        }, timeout=300)
        with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
            names = [n for n in zf.namelist() if n.endswith(".json")]
            return json.loads(zf.read(names[0]))
    except requests.exceptions.Timeout:
        raise RuntimeError(
            "BLASTp timed out — NCBI nr database queries can take 10+ minutes under load. "
            "Try again; the job may still be running on NCBI's servers."
        )


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
    try:
        r = requests.post(f"{url}/run", data={
            "email": EMAIL, "sequence": sequence,
            "database": "pfam", "E": "1.0",
        }, timeout=60)
        r.raise_for_status()
        jid = r.text.strip()
        if _ebi_poll(url, jid) != "FINISHED":
            raise RuntimeError("Job did not finish successfully")
        r = requests.get(f"{url}/result/{jid}/out", timeout=120)
        r.raise_for_status()
        return r.text
    except requests.exceptions.Timeout:
        raise RuntimeError(
            "HMMER timed out — EBI servers may be slow. Try again in a few minutes."
        )


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
        s = requests.get(f"{url}/ticket/{ticket}", timeout=60).json().get("status", "UNKNOWN")
        if s in ("COMPLETE", "ERROR", "FAILED", "UNKNOWN"):
            break
        time.sleep(10)
    if s != "COMPLETE":
        raise RuntimeError(f"FoldSeek status: {s}")
    r = requests.get(f"{url}/result/download/{ticket}", timeout=120)
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
    try:
        bo2  = data["BlastOutput2"]
        bo2  = bo2[0] if isinstance(bo2, list) else bo2
        hits = bo2["report"]["results"]["search"]["hits"]
    except (KeyError, IndexError, TypeError):
        hits = []
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
    st.markdown(_cards(("Top hits", len(df), "nr database")), unsafe_allow_html=True)
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
    for line in text.splitlines():
        if "Scores for complete sequence" in line:
            in_hits = True
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
    above = int((df["Significant"] == "Yes").sum()) if not df.empty else 0
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


def _calc_dssp(pdb_text: str):
    """Return {H, E, C: pct} via DSSP, or None if the binary is unavailable."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", io.StringIO(pdb_text))
    model = structure[0]
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        pdbio.save(f)
        tmppath = f.name
    try:
        dssp_obj = None
        for exe in ("mkdssp", "dssp"):
            try:
                dssp_obj = DSSP(model, tmppath, dssp=exe)
                break
            except Exception:
                continue
        if dssp_obj is None:
            return None
        counts = {"H": 0, "E": 0, "C": 0}
        for key in dssp_obj:
            ss = dssp_obj[key][2]
            if ss in ("H", "G", "I"):
                counts["H"] += 1
            elif ss in ("E", "B"):
                counts["E"] += 1
            else:
                counts["C"] += 1
        total = sum(counts.values())
        return {k: round(v / total * 100, 1) for k, v in counts.items()} if total else None
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


# ── Structure tab — main display function ─────────────────────────────────────

def show_structure(pdb_text: str, fasta_text, phobius_text=None) -> None:
    plddt  = _parse_pdb_plddts(pdb_text)
    mean_pl = round(sum(plddt) / len(plddt), 1) if plddt else 0
    seq    = _parse_fasta_seq(fasta_text) if fasta_text else ""

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

    ss = _calc_dssp(pdb_text)

    col_v, col_a = st.columns([5, 4])

    with col_v:
        components.html(_3dmol_html(pdb_text, height=440), height=494, scrolling=False)
        if plddt:
            st.markdown(_plddt_bar_html(plddt), unsafe_allow_html=True)

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

# ── Two-column layout ──────────────────────────────────────────────────────────

col_left, col_right = st.columns([3, 9], gap="small")

# ── Left panel — input + status ────────────────────────────────────────────────

with col_left:
    st.markdown(
        '<p style="color:#3b82f6;font-size:9.5px;font-weight:700;'
        'letter-spacing:0.12em;text-transform:uppercase;margin-bottom:16px;">Input</p>',
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
        "Run Pipeline",
        type="primary",
        disabled=not (fasta_file or pdb_file),
    )

    if run:
        fasta_text = fasta_file.read().decode() if fasta_file else None
        pdb_text   = pdb_file.read().decode()   if pdb_file   else None

        tasks: dict = {}
        if fasta_text:
            tasks["InterProScan"] = (run_interproscan, fasta_text)
            tasks["BLASTp"]       = (run_blast,        fasta_text)
            tasks["Phobius"]      = (run_phobius,       fasta_text)
            tasks["HMMER"]        = (run_hmmer,         fasta_text)
        if pdb_text:
            tasks["FoldSeek"]     = (run_foldseek,      pdb_text)

        st.session_state["active_tools"] = list(tasks.keys())

        results: dict = {}
        with st.spinner(f"Running {len(tasks)} tools in parallel — typically 8–12 min…"):
            with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
                futures = {pool.submit(fn, arg): name for name, (fn, arg) in tasks.items()}
                for future in as_completed(futures):
                    name = futures[future]
                    try:
                        results[name] = {"ok": True,  "data": future.result()}
                    except Exception as exc:
                        results[name] = {"ok": False, "error": str(exc)}

        st.session_state["results"]    = results
        st.session_state["fasta_text"] = fasta_text
        st.session_state["pdb_text"]   = pdb_text

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

# ── Right panel — results ──────────────────────────────────────────────────────

with col_right:
    _res    = st.session_state.get("results",      {})
    _active = st.session_state.get("active_tools", [])

    _pdb_ss   = st.session_state.get("pdb_text")
    _fasta_ss = st.session_state.get("fasta_text")

    if not _res and not _pdb_ss:
        # Empty-state placeholder
        st.markdown("""
        <div style="display:flex;flex-direction:column;align-items:center;
                    justify-content:center;min-height:360px;gap:14px;opacity:0.55;">
          <div style="width:52px;height:52px;background:#080c17;border:1px solid #1e2d4a;
                      border-radius:14px;display:flex;align-items:center;
                      justify-content:center;">
            <svg width="22" height="22" viewBox="0 0 24 24" fill="none"
                 stroke="#1e2d4a" stroke-width="1.8"
                 stroke-linecap="round" stroke-linejoin="round">
              <circle cx="12" cy="12" r="10"/>
              <line x1="12" y1="8" x2="12" y2="12"/>
              <line x1="12" y1="16" x2="12.01" y2="16"/>
            </svg>
          </div>
          <div style="text-align:center;">
            <p style="color:#334155;font-size:14px;font-weight:500;margin:0 0 4px;">
              No results yet
            </p>
            <p style="color:#1e2d4a;font-size:12px;margin:0;">
              Upload a FASTA file and click Run Pipeline to begin annotation
            </p>
          </div>
        </div>
        """, unsafe_allow_html=True)
    else:
        if _res:
            n_ok    = sum(r["ok"] for r in _res.values())
            n_total = len(_res)
            all_ok  = n_ok == n_total
            bar_col = "#22c55e" if all_ok else "#f59e0b"
            bar_bg  = "rgba(34,197,94,0.06)"  if all_ok else "rgba(245,158,11,0.06)"
            bar_bdr = "rgba(34,197,94,0.15)"  if all_ok else "rgba(245,158,11,0.15)"
            st.markdown(
                f'<div style="background:{bar_bg};border:1px solid {bar_bdr};border-radius:8px;'
                f'padding:10px 16px;margin-bottom:18px;display:flex;align-items:center;gap:8px;">'
                f'<span style="color:{bar_col};font-size:12.5px;font-weight:600;">'
                f'&#x2714;&ensp;{n_ok}/{n_total} tools completed</span></div>',
                unsafe_allow_html=True,
            )
            # Inline error banners for failed tools
            for name in _active:
                r = _res.get(name)
                if r and not r["ok"]:
                    st.markdown(
                        f'<div style="background:rgba(239,68,68,0.05);'
                        f'border:1px solid rgba(239,68,68,0.15);border-radius:6px;'
                        f'padding:10px 14px;margin-bottom:8px;color:#ef4444;font-size:12px;">'
                        f'<strong>{_esc(name)}</strong> — {_esc(r["error"])}</div>',
                        unsafe_allow_html=True,
                    )

        # Build tab list: Structure first (if PDB uploaded), then tool results
        finished  = [t for t in _active if _res.get(t, {}).get("ok")]
        tab_names = (["🧬 Structure"] if _pdb_ss else []) + finished

        if tab_names:
            tabs   = st.tabs(tab_names)
            offset = 0
            if _pdb_ss:
                with tabs[0]:
                    _phobius_data = (
                        _res["Phobius"]["data"]
                        if _res.get("Phobius", {}).get("ok") else None
                    )
                    show_structure(_pdb_ss, _fasta_ss, _phobius_data)
                offset = 1
            for i, name in enumerate(finished):
                with tabs[offset + i]:
                    data = _res[name]["data"]
                    if   name == "InterProScan": show_interproscan(data)
                    elif name == "BLASTp":       show_blast(data)
                    elif name == "Phobius":      show_phobius(data)
                    elif name == "HMMER":        show_hmmer(data)
                    elif name == "FoldSeek":     show_foldseek(data)

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
