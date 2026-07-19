#!/usr/bin/env python3
"""
Dark Proteome Pipeline — Fold Classifier Validation Harness
============================================================
Loads  validation/ground_truth.json, finds PDB files in data/input/,
runs the structural fold classifier on each available protein, and prints:

  • Per-domain table: predicted vs expected fold, structural score, result
  • Per-domain fold accuracy (excluding vocab-gap domains and skipped proteins)
  • Outcome-category accuracy

Run from the project root:
    python3 validation/run_validation.py

What this script does in plain English:
  - For each of the 15 dissertation proteins, it looks for a matching PDB file.
  - If no PDB is found, it prints "SKIPPED" and moves on (never crashes).
  - For each found protein, it uses DSSP geometry (beta-strand count, spacing,
    orientation, Cα angles) to score 15 possible fold types.
  - It uses ONLY structural evidence — no InterProScan or BLAST results —
    so the score reflects what the 3D shape alone can tell us.
  - If the highest structural score is below 40/100, it abstains ("No confident
    fold assignment") instead of guessing.
  - It compares the prediction to the expert label in ground_truth.json and
    marks each domain CORRECT, WRONG, ABSTAIN✓, FALSE_CALL, or VOCAB_GAP
    (VOCAB_GAP = fold type not yet in the classifier vocabulary).
"""

import io
import json
import math
import os
import pathlib
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, PDBIO

try:
    from Bio.PDB.DSSP import DSSP as _DSSP
except Exception:
    _DSSP = None

ROOT = pathlib.Path(__file__).parent.parent

# ─────────────────────────────────────────────────────────────────────────────
#  Structural-analysis helpers (pure Python / BioPython, no Streamlit)
# ─────────────────────────────────────────────────────────────────────────────

ABSTAIN_THRESHOLD = 40  # max structural score below this → no confident assignment

# Keyword sets used by scoring functions (kept in sync with dark_proteome_app.py)
_BARREL_ACCESSIONS    = {"PS51208", "SSF103515", "PF00291"}
_AT_ACCESSIONS        = {"PF20696", "PF27415", "PF22399", "PS51208", "PF00291"}
_RHS_ACCESSIONS       = {"PF03917", "PF05593", "PF19434", "PF13796"}
_FHA_ACCESSIONS       = {"PF02413", "PF13312", "PF15609"}
_FN3_IG_ACCESSIONS    = {"PF00041", "PF00047", "PF13927", "SSF48726"}
_TPR_ACCESSIONS       = {"PF00515", "PF07719", "PF07720", "PF07721", "PF13176", "PF13174"}
_TIM_ACCESSIONS       = {"SSF51351", "G3DSA:3.20.20.70"}
_ROSSMANN_ACCESSIONS  = {"SSF52374", "SSF51735", "G3DSA:3.40.50.720"}
_OB_ACCESSIONS        = {"SSF50249", "PF01336"}
_PROPELLER_ACCESSIONS = {"PF00400", "PF01344", "PF01966", "SSF50978"}
_LECTIN_ACCESSIONS    = {"SSF49899", "PF00139", "PF00652", "PF14200"}
_COIL_ACCESSIONS      = {"PF05765", "SSF47370"}
_HEAT_ACCESSIONS      = {"PF02985", "PF00514", "PF13190", "PF19012"}
_CALYCIN_ACCESSIONS   = {"PF00061", "SSF50814"}


def _ipr_hits_for(ipr_data, accession_set):
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


def _blast_top_title(blast_data):
    if not blast_data:
        return ""
    for result in blast_data.get("BlastOutput2", []):
        hits = result.get("report", {}).get("results", {}).get("search", {}).get("hits", [])
        if hits:
            return hits[0].get("description", [{}])[0].get("title", "").lower()
    return ""


def _find_barrel_cutoff(ipr_data):
    if not ipr_data:
        return 0, None
    cutoff, label = None, None
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
                        label  = f"{acc} ({lib}: {name})" if name else f"{acc} ({lib})"
    return (cutoff or 0), label


def _rtx_matches(seq):
    """Return list of (start, end) for RTX nonapeptide GGXGXDXUX motifs."""
    if not seq:
        return []
    out = []
    for i in range(len(seq) - 8):
        w = seq[i:i+9].upper()
        if (w[0] == 'G' and w[1] == 'G' and w[3] == 'G' and w[5] == 'D'
                and w[7] in 'LIVMFWACYP'):
            out.append((i + 1, i + 9))
    return out


def _run_dssp_raw(pdb_text):
    """Return (ss_list, method_name).  ss_list is ['H','E','C',...] or None."""
    # 1. pydssp (pure Python, best accuracy, optional dependency)
    try:
        import pydssp
        raw    = pydssp.read_pdbtext(pdb_text)
        coords = raw[0] if isinstance(raw, (list, tuple)) else raw
        ss_arr = pydssp.assign(coords, out_type="c3")
        return [str(s) for s in ss_arr], "pydssp"
    except Exception:
        pass

    # 2. Phi/ψ dihedral heuristic via BioPython (always available)
    try:
        from Bio.PDB.Polypeptide import PPBuilder
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("p", io.StringIO(pdb_text))
        ppb    = PPBuilder()
        ss     = []
        for pp in ppb.build_peptides(struct):
            for phi, psi in pp.get_phi_psi_list():
                if phi is None or psi is None:
                    ss.append("C")
                    continue
                ph, ps = math.degrees(phi), math.degrees(psi)
                if -90 <= ph <= -30 and -77 <= ps <= -17:
                    ss.append("H")
                elif -170 <= ph <= -50 and (100 <= ps <= 180 or -180 <= ps <= -150):
                    ss.append("E")
                else:
                    ss.append("C")
        if ss:
            return ss, "Ramachandran approximation (lower accuracy)"
    except Exception:
        pass

    # 3. mkdssp binary via BioPython
    parser = PDBParser(QUIET=True)
    struct  = parser.get_structure("p", io.StringIO(pdb_text))
    model   = struct[0]
    pdbio   = PDBIO(); pdbio.set_structure(struct)
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
        pdbio.save(f); tmp = f.name
    try:
        for exe in ("mkdssp", "dssp"):
            try:
                if _DSSP is None: continue
                d = _DSSP(model, tmp, dssp=exe)
                return [d[k][2] for k in d], "DSSP"
            except Exception:
                continue
        return None, "none"
    finally:
        os.unlink(tmp)


def _run_dssp_strands(pdb_text):
    """Return (strand_list, ss_pct_dict, method_name)."""
    ss, method = _run_dssp_raw(pdb_text)
    if ss is None:
        return [], None, "none"
    strands, in_s, s0 = [], False, 0
    for i, c in enumerate(ss):
        beta = c in ("E", "B")
        if beta and not in_s:   in_s, s0 = True, i + 1
        elif not beta and in_s: in_s = False; strands.append((s0, i))
    if in_s: strands.append((s0, len(ss)))
    counts = {"H": 0, "E": 0, "C": 0}
    for c in ss:
        if c in ("H", "G", "I"): counts["H"] += 1
        elif c in ("E", "B"):    counts["E"] += 1
        else:                    counts["C"] += 1
    tot = sum(counts.values())
    pct = {k: round(v/tot*100, 1) for k, v in counts.items()} if tot else None
    return strands, pct, method


def _parse_pdb_plddts(pdb_text):
    """Extract per-residue pLDDT (B-factor of Cα, first model)."""
    p = PDBParser(QUIET=True)
    s = p.get_structure("p", io.StringIO(pdb_text))
    scores = []
    for model in s:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    scores.append(res["CA"].get_bfactor())
        break
    return scores


def _is_alphafold_pdb(pdb_text, plddts=None):
    """Return True when this looks like an AlphaFold/ESMFold structure.
    Checks header keywords first; falls back to B-factor range heuristic."""
    for line in pdb_text.splitlines()[:30]:
        lo = line.lower()
        if any(m in lo for m in ("alphafold", "esm", "colabfold", "af2", "esmfold")):
            return True
    scores = plddts if plddts is not None else _parse_pdb_plddts(pdb_text)
    return bool(scores and all(0.0 <= s <= 100.0 for s in scores))


def _parse_pdb_ca_coords(pdb_text):
    p = PDBParser(QUIET=True)
    s = p.get_structure("p", io.StringIO(pdb_text))
    coords = []
    for model in s:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    v = res["CA"].get_vector()
                    coords.append((float(v[0]), float(v[1]), float(v[2])))
        break
    return coords


# ── Scoring functions (structural_score, sequence_score, evidence_list) ───────
# [STR] prefix = structural evidence (DSSP/geometry)
# [SEQ] prefix = sequence evidence (IPR/BLAST) — shown as corroborating only


def _score_at(ipr_data, blast_data, seq, gap_sd):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _AT_ACCESSIONS):
        seq_s = min(100, seq_s + 25)
        hits.append(f"[SEQ] {acc}" + (f" {name}" if name else ""))
    for kw in ["autotransporter", "type v secretion", "outer membrane protein"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 15); hits.append(f"[SEQ] BLAST:{kw}"); break
    if seq and seq.rstrip()[-1].upper() == "F":
        seq_s = min(100, seq_s + 10); hits.append("[SEQ] C-term Phe (BAM signal)")
    if 0 < gap_sd < 8:
        str_s = min(100, str_s + 35); hits.append(f"[STR] gap SD={gap_sd:.1f}<8 (regular β-helix)")
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_rhs(ipr_data, blast_data, seq, gap_sd, par_pct, antipar_pct, first_s):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _RHS_ACCESSIONS):
        seq_s = min(100, seq_s + 30); hits.append(f"[SEQ] {acc}")
        if seq_s >= 60: break
    if seq and len(seq) > 1500:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {len(seq)} aa >1500")
    if gap_sd > 20:
        str_s = min(100, str_s + 15); hits.append(f"[STR] irregular SD={gap_sd:.1f}")
    if par_pct < 60 and antipar_pct < 60:
        str_s = min(100, str_s + 10); hits.append("[STR] mixed orientation")
    if first_s > 200:
        str_s = min(100, str_s + 10); hits.append(f"[STR] first strand at {first_s}")
    for kw in ["rhs repeat", "type vi secretion"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 15); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_fha(ipr_data, blast_data, seq, gap_sd):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _FHA_ACCESSIONS):
        seq_s = min(100, seq_s + 25); hits.append(f"[SEQ] {acc}")
        if seq_s >= 75: break
    if seq and len(seq) > 2000:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {len(seq)} aa >2000")
    if 0 < gap_sd < 8:
        str_s = min(100, str_s + 25); hits.append(f"[STR] SD={gap_sd:.1f}<8")
    for kw in ["two-partner secretion", "filamentous hemagglutinin"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 15); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_fn3(ipr_data, blast_data, n_strands, antipar_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _FN3_IG_ACCESSIONS)[:4]:
        seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] {acc}")
    if n_strands >= 6:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {n_strands} strands")
    if antipar_pct > 50:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {antipar_pct}% antiparallel")
    for kw in ["fibronectin", "immunoglobulin"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 10); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_tpr(ipr_data, blast_data, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _TPR_ACCESSIONS):
        seq_s = min(100, seq_s + 25); hits.append(f"[SEQ] {acc}")
        if seq_s >= 75: break
    if ss_pct and ss_pct.get("H", 0) > 40:
        str_s = min(100, str_s + 25); hits.append(f"[STR] {ss_pct['H']}% helix")
    if ss_pct and ss_pct.get("E", 100) < 10:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {ss_pct.get('E',0)}% strand")
    for kw in ["tetratricopeptide", "tpr repeat"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 15); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_tim(ipr_data, blast_data, seq, n_strands, par_pct, mean_gap, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _TIM_ACCESSIONS):
        seq_s = min(100, seq_s + 35); hits.append(f"[SEQ] {acc}")
        if seq_s >= 70: break
    if 6 <= n_strands <= 10:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {n_strands} strands (TIM:8)")
    if par_pct > 60:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {par_pct}% parallel")
    if 10 <= mean_gap <= 45:
        str_s = min(100, str_s + 10); hits.append(f"[STR] gap {mean_gap:.1f} aa")
    if seq and 200 <= len(seq) <= 450:
        str_s = min(100, str_s + 10); hits.append(f"[STR] size {len(seq)} aa")
    for kw in ["tim barrel", "triosephosphate"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_rossmann(ipr_data, blast_data, ss_pct, par_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _ROSSMANN_ACCESSIONS):
        seq_s = min(100, seq_s + 30); hits.append(f"[SEQ] {acc}")
        if seq_s >= 60: break
    if ss_pct and ss_pct.get("H", 0) > 25 and ss_pct.get("E", 0) > 15:
        str_s = min(100, str_s + 15); hits.append(f"[STR] mixed α/β")
    if par_pct > 50:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {par_pct}% parallel")
    for kw in ["rossmann", "nad-binding", "nucleotide binding"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_ob(ipr_data, blast_data, seq, n_strands, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _OB_ACCESSIONS):
        seq_s = min(100, seq_s + 35); hits.append(f"[SEQ] {acc}")
    if seq and len(seq) < 200:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {len(seq)} aa (OB small)")
    if 4 <= n_strands <= 7:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {n_strands} strands (OB:5)")
    if ss_pct and ss_pct.get("E", 0) > 30:
        str_s = min(100, str_s + 10); hits.append(f"[STR] {ss_pct['E']}% strand")
    for kw in ["ob fold", "oligonucleotide"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_propeller(ipr_data, blast_data, n_strands, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _PROPELLER_ACCESSIONS):
        seq_s = min(100, seq_s + 25); hits.append(f"[SEQ] {acc}")
    if n_strands >= 16:
        str_s = min(100, str_s + 25); hits.append(f"[STR] {n_strands} strands (≥16)")
    if ss_pct and ss_pct.get("E", 0) > 25:
        str_s = min(100, str_s + 10); hits.append(f"[STR] {ss_pct['E']}% strand")
    if ss_pct and ss_pct.get("H", 0) < 20:
        str_s = min(100, str_s + 10); hits.append(f"[STR] low helix {ss_pct.get('H',0)}%")
    for kw in ["wd40", "kelch", "beta-propeller"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_lectin(ipr_data, blast_data, ss_pct, antipar_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _LECTIN_ACCESSIONS):
        seq_s = min(100, seq_s + 30); hits.append(f"[SEQ] {acc}")
    if ss_pct and ss_pct.get("E", 0) > 25:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {ss_pct['E']}% strand")
    if antipar_pct > 50:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {antipar_pct}% antiparallel")
    for kw in ["lectin", "concanavalin", "jelly roll"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_coil(ipr_data, blast_data, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _COIL_ACCESSIONS):
        seq_s = min(100, seq_s + 35); hits.append(f"[SEQ] {acc}")
    if ss_pct and ss_pct.get("H", 0) > 70:
        str_s = min(100, str_s + 35); hits.append(f"[STR] {ss_pct['H']}% helix (>70%)")
    if ss_pct and ss_pct.get("E", 100) < 5:
        str_s = min(100, str_s + 15); hits.append(f"[STR] {ss_pct.get('E',0)}% strand")
    for kw in ["coiled-coil", "leucine zipper"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 15); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_heat(ipr_data, blast_data, seq, ss_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _HEAT_ACCESSIONS):
        seq_s = min(100, seq_s + 25); hits.append(f"[SEQ] {acc}")
        if seq_s >= 50: break
    if ss_pct and ss_pct.get("H", 0) > 60:
        str_s = min(100, str_s + 25); hits.append(f"[STR] {ss_pct['H']}% helix")
    if ss_pct and ss_pct.get("E", 100) < 15:
        str_s = min(100, str_s + 10); hits.append(f"[STR] low strand {ss_pct.get('E',0)}%")
    if seq and len(seq) > 500:
        str_s = min(100, str_s + 10); hits.append(f"[STR] {len(seq)} aa")
    for kw in ["heat repeat", "armadillo"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _score_calycin(ipr_data, blast_data, seq, n_strands, antipar_pct):
    str_s, seq_s, hits = 0, 0, []
    for acc, lib, name in _ipr_hits_for(ipr_data, _CALYCIN_ACCESSIONS):
        seq_s = min(100, seq_s + 45); hits.append(f"[SEQ] {acc}")
        if seq_s >= 90: break
    if seq and len(seq) < 200:
        str_s = min(100, str_s + 10); hits.append(f"[STR] {len(seq)} aa")
    if 6 <= n_strands <= 10:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {n_strands} strands (calycin:8)")
    if antipar_pct > 60:
        str_s = min(100, str_s + 20); hits.append(f"[STR] {antipar_pct}% antiparallel")
    for kw in ["lipocalin", "fatty acid-binding"]:
        if kw in _blast_top_title(blast_data):
            seq_s = min(100, seq_s + 20); hits.append(f"[SEQ] BLAST:{kw}"); break
    return min(100, str_s + seq_s), str_s, seq_s, hits


def _classify_segment(strands, ca_coords, plddt, seq, ss_pct,
                      ipr_data=None, blast_data=None):
    """Classify one protein segment (whole protein, or a domain slice).

    Returns dict with:
      best_fold         — predicted fold name (or 'No confident fold assignment')
      structural_score  — highest structural-only score (the decision-driving number)
      sequence_score    — sequence component for the winning fold (corroborating)
      fold_scores_str   — {fold: structural_score} for all 15 folds
      n_strands, ss_pct
    """
    n_total   = len(strands)
    rtx_hits  = _rtx_matches(seq) if seq else []

    # Score helix-only folds (available even with 0 strands)
    _, tpr_str,  tpr_seq,  _ = _score_tpr(ipr_data, blast_data, ss_pct)
    _, coil_str, coil_seq, _ = _score_coil(ipr_data, blast_data, ss_pct)
    _, heat_str, heat_seq, _ = _score_heat(ipr_data, blast_data, seq, ss_pct)

    if n_total < 2:
        helix_str = {"TPR solenoid": tpr_str, "Coiled-coil": coil_str,
                     "HEAT/ARM solenoid": heat_str}
        best = max(helix_str, key=helix_str.get)
        best_str = helix_str[best]
        if best_str < ABSTAIN_THRESHOLD:
            return dict(best_fold="No confident fold assignment",
                        structural_score=best_str, sequence_score=0,
                        fold_scores_str={}, n_strands=n_total, ss_pct=ss_pct)
        seq_map = {"TPR solenoid": tpr_seq, "Coiled-coil": coil_seq,
                   "HEAT/ARM solenoid": heat_seq}
        return dict(best_fold=best, structural_score=best_str,
                    sequence_score=seq_map[best],
                    fold_scores_str=helix_str, n_strands=n_total, ss_pct=ss_pct)

    barrel_cutoff, _ = _find_barrel_cutoff(ipr_data)
    ws = [(s, e) for s, e in strands if s < barrel_cutoff] if barrel_cutoff > 0 else strands
    wc = ca_coords[:barrel_cutoff-1] if barrel_cutoff and ca_coords else ca_coords
    wp = plddt[:barrel_cutoff-1] if barrel_cutoff and plddt else plddt
    if len(ws) < 2: ws, wc, wp, barrel_cutoff = strands, ca_coords, plddt, 0

    n  = len(ws)
    fs = ws[0][0] if ws else 0
    lengths = [e - s + 1 for s, e in ws]
    gaps    = [ws[i+1][0] - ws[i][1] - 1 for i in range(n - 1)]
    mg      = sum(gaps) / len(gaps) if gaps else 0
    gsd     = math.sqrt(sum((g - mg)**2 for g in gaps) / len(gaps)) if len(gaps) > 1 else 0

    # Item 7: repeat regularity — require SD < 30% of mean gap OR < 4 aa absolute.
    # The absolute floor prevents false-penalising tight RTX/FHA repeats where
    # mg ≈ 3–5 aa and any small variance looks proportionally large.
    regular_repeat = (gsd < max(4.0, 0.30 * mg)) if mg > 0 else True

    def _norm(v):
        m = math.sqrt(sum(x*x for x in v))
        return tuple(x/m for x in v) if m else (0., 0., 0.)

    dirs = []
    for s, e in ws:
        s0, e0 = s-1, e-1
        if 0 <= s0 < len(wc) and 0 <= e0 < len(wc):
            v = tuple(wc[e0][j] - wc[s0][j] for j in range(3))
            dirs.append(_norm(v))
        else:
            dirs.append(None)

    angles = []
    for i in range(len(dirs)-1):
        if dirs[i] and dirs[i+1]:
            ca = max(-1., min(1., sum(a*b for a, b in zip(dirs[i], dirs[i+1]))))
            angles.append(math.degrees(math.acos(ca)))

    par_pct     = round(sum(1 for a in angles if a < 45)  / len(angles)*100) if angles else 0
    antipar_pct = round(sum(1 for a in angles if a > 135) / len(angles)*100) if angles else 0
    rtx_on      = any(not (re < s-3 or rs > e+3) for rs, re in rtx_hits for s, e in ws)

    spl = [wp[i] for s, e in ws for i in range(s-1, min(e, len(wp)))]
    spl_mean = round(sum(spl)/len(spl), 1) if spl else 0

    # pLDDT > 70 for repeat regions (Item 7)
    high_plddt_repeats = spl_mean > 70

    # Score all folds
    _, at_str,   at_seq,   _ = _score_at(ipr_data, blast_data, seq, gsd)
    _, rhs_str,  rhs_seq,  _ = _score_rhs(ipr_data, blast_data, seq, gsd,
                                           par_pct, antipar_pct, fs)
    _, fha_str,  fha_seq,  _ = _score_fha(ipr_data, blast_data, seq, gsd)
    _, fn3_str,  fn3_seq,  _ = _score_fn3(ipr_data, blast_data, n_total, antipar_pct)
    _, tim_str,  tim_seq,  _ = _score_tim(ipr_data, blast_data, seq, n_total,
                                           par_pct, mg, ss_pct)
    _, ros_str,  ros_seq,  _ = _score_rossmann(ipr_data, blast_data, ss_pct, par_pct)
    _, ob_str,   ob_seq,   _ = _score_ob(ipr_data, blast_data, seq, n_total, ss_pct)
    _, prop_str, prop_seq, _ = _score_propeller(ipr_data, blast_data, n_total, ss_pct)
    _, lec_str,  lec_seq,  _ = _score_lectin(ipr_data, blast_data, ss_pct, antipar_pct)
    _, cal_str,  cal_seq,  _ = _score_calycin(ipr_data, blast_data, seq, n_total, antipar_pct)

    # Beta-solenoid (purely geometric, Item 7: require regular spacing + high pLDDT)
    sol = min(30, n * 4) + max(0, int(30 - gsd * 6))
    sol += round(max(par_pct, antipar_pct) * 0.2) if angles else 0
    sol += 15 if rtx_on else 0
    sol += 5  if high_plddt_repeats else 0
    if not regular_repeat: sol = max(0, sol - 20)  # Item 7: penalise irregular spacing
    sol = min(100, max(0, sol))

    # Beta-barrel (purely geometric)
    brl = (30 if 8 <= n_total <= 24 else 0) + int(antipar_pct * 0.25)
    brl += 20 if 1 <= mg <= 5 else 0
    brl += max(0, int(20 - gsd * 3))
    brl += 15 if 5 <= sum(lengths)/n <= 12 else 0
    brl = min(100, max(0, brl))

    fold_scores_str = {
        "Beta-solenoid":           sol,
        "RHS solenoid":            rhs_str,
        "Autotransporter β-helix": at_str,
        "FHA β-helix":             fha_str,
        "FN3/Ig sandwich":         fn3_str,
        "Lectin/jelly-roll":       lec_str,
        "Beta-barrel":             brl,
        "Beta-propeller":          prop_str,
        "TIM barrel":              tim_str,
        "OB fold":                 ob_str,
        "Calycin/lipocalin":       cal_str,
        "Rossmann fold":           ros_str,
        "Coiled-coil":             coil_str,
        "HEAT/ARM solenoid":       heat_str,
        "TPR solenoid":            tpr_str,
    }
    seq_map = {
        "Beta-solenoid": 0, "RHS solenoid": rhs_seq,
        "Autotransporter β-helix": at_seq, "FHA β-helix": fha_seq,
        "FN3/Ig sandwich": fn3_seq, "Lectin/jelly-roll": lec_seq,
        "Beta-barrel": 0, "Beta-propeller": prop_seq,
        "TIM barrel": tim_seq, "OB fold": ob_seq,
        "Calycin/lipocalin": cal_seq, "Rossmann fold": ros_seq,
        "Coiled-coil": coil_seq, "HEAT/ARM solenoid": heat_seq,
        "TPR solenoid": tpr_seq,
    }

    best = max(fold_scores_str, key=fold_scores_str.get)
    best_str_score = fold_scores_str[best]

    if best_str_score < ABSTAIN_THRESHOLD:
        return dict(best_fold="No confident fold assignment",
                    structural_score=best_str_score, sequence_score=0,
                    fold_scores_str=fold_scores_str, n_strands=n_total, ss_pct=ss_pct)

    return dict(best_fold=best, structural_score=best_str_score,
                sequence_score=seq_map.get(best, 0),
                fold_scores_str=fold_scores_str, n_strands=n_total, ss_pct=ss_pct)


# ─────────────────────────────────────────────────────────────────────────────
#  pLDDT-based domain segmentation
# ─────────────────────────────────────────────────────────────────────────────

def _segment_by_plddt(plddt, min_domain=50, floor=70, gap_min=10):
    """Find boundaries where pLDDT drops below 'floor' for ≥ gap_min residues.
    Returns list of (start, end) in 1-based residue coordinates."""
    n = len(plddt)
    if n == 0:
        return []
    in_low, low_start = False, 0
    cuts = []
    for i, s in enumerate(plddt):
        if s < floor and not in_low:
            in_low, low_start = True, i
        elif s >= floor and in_low:
            in_low = False
            if (i - low_start) >= gap_min:
                cuts.append((low_start + i) // 2)
    splits = [0] + cuts + [n]
    segs = []
    for i in range(len(splits)-1):
        s, e = splits[i], splits[i+1]
        if (e - s) >= min_domain:
            segs.append((s+1, e))
    return segs or [(1, n)]


# ─────────────────────────────────────────────────────────────────────────────
#  Fold-name matching: ground-truth label → accepted classifier outputs
# ─────────────────────────────────────────────────────────────────────────────

# None = vocabulary gap (classifier has no class for this fold)
_GT_TO_PRED = {
    "GA3-like domain":                                      None,
    "beta-roll (RTX-type)":                                 ["Beta-solenoid"],
    "beta-helix solenoid (Type Va autotransporter)":        ["Autotransporter β-helix", "FHA β-helix"],
    "Ig-like beta-sandwich (tandem repeat)":                ["FN3/Ig sandwich"],
    "jelly-roll (CE2 esterase N-domain)":                   ["Lectin/jelly-roll"],
    "lectin beta-sandwich (x5 tandem, ConA-like)":          ["Lectin/jelly-roll", "FN3/Ig sandwich"],
    "beta-sandwich (Ig-like)":                              ["FN3/Ig sandwich"],
    "beta-helix repeat":                                    ["Beta-solenoid", "Autotransporter β-helix", "FHA β-helix"],
    "RTX beta-propeller (novel radial variant, 4-bladed)":  ["Beta-propeller", "Beta-solenoid"],
    "beta-helix (polysaccharide lyase-like, tandem x4-6)":  ["Beta-solenoid", "Autotransporter β-helix"],
    "zinc-metalloprotease fold (Peptidase_M26 catalytic)":  None,
    "G5 beta-sandwich":                                     ["FN3/Ig sandwich"],
    "LPXTG sortase anchor":                                 None,
    "beta-helix (FHA-type)":                                ["FHA β-helix", "Autotransporter β-helix"],
    "No confident fold assignment":                         ["No confident fold assignment"],
    "coiled-coil (WXG100/EsxA-like)":                      ["Coiled-coil"],
    "RHS beta-solenoid":                                    ["RHS solenoid"],
    "FN3/Ig beta-sandwich stalk":                           ["FN3/Ig sandwich"],
    "pyocin knob beta-propeller (C-term RBD)":              ["Beta-propeller"],
    "pentapeptide-repeat solenoid (left-handed beta-helix variant)": ["Beta-solenoid"],
    "TPR solenoid (alpha-alpha superhelix)":                ["TPR solenoid"],
    "calycin beta-barrel (lipocalin/FABP fold)":            ["Calycin/lipocalin"],
    "beta-helix solenoid passenger (pectin lyase-like, x3 tandem)": ["Autotransporter β-helix", "Beta-solenoid"],
    "beta-barrel translocator (Type Va)":                   ["Beta-barrel"],
}


def _check(predicted, expected_fold):
    """Return 'CORRECT', 'WRONG', or 'VOCAB_GAP'."""
    accepted = _GT_TO_PRED.get(expected_fold, [])
    if accepted is None:
        return "VOCAB_GAP"
    return "CORRECT" if predicted in accepted else "WRONG"


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    gt_path  = ROOT / "validation" / "ground_truth.json"
    pdb_dir  = ROOT / "data" / "input"

    with open(gt_path) as f:
        gt = json.load(f)

    proteins = gt["proteins"]

    W = [13, 4, 17, 34, 28, 5, 5, 11]
    hdr = (f"  {'Protein':<{W[0]}} {'Dom':>{W[1]}} {'Residues':<{W[2]}} "
           f"{'Expected fold':<{W[3]}} {'Predicted fold':<{W[4]}} "
           f"{'Str':>{W[5]}} {'Seq':>{W[6]}} {'Result':<{W[7]}}")
    rule = "  " + "─" * (sum(W) + 7 * 2)

    print()
    print(hdr)
    print(rule)

    n_class = n_ok = 0
    outcome_tbl: dict = {}

    for prot in proteins:
        name     = prot["name"]
        pdb_path = pdb_dir / f"{name.upper()}.pdb"

        if not pdb_path.exists():
            print(f"  {'SKIPPED':13} — no PDB found for {name}")
            continue

        pdb_text = pdb_path.read_text()
        plddts   = _parse_pdb_plddts(pdb_text)
        is_af    = _is_alphafold_pdb(pdb_text, plddts)
        if not is_af:
            print(f"\n  WARNING: {name}.pdb — B-factors do not look like pLDDT "
                  f"(possibly experimental; confidence track may be misleading)\n")

        ca_coords = _parse_pdb_ca_coords(pdb_text)
        strands_all, ss_pct, ss_method = _run_dssp_strands(pdb_text)
        print(f"\n  {name} ({len(plddts)} aa) | SS: {ss_method} | "
              f"{len(strands_all)} β-strands")
        print(rule)

        domains = prot["domains"]
        for d_idx, dom in enumerate(domains, 1):
            expected  = dom["expected_fold"]
            rr        = dom.get("residue_range")
            do_abstain = dom.get("expected_abstain_domain", False) or prot.get("expected_abstain", False)
            outcome   = prot["outcome"]

            if rr:
                r0, r1 = rr[0], rr[1]
                res_lbl = f"{r0}–{r1}"
                seg_str = [(s, e) for s, e in strands_all if s >= r0 and e <= r1]
                seg_ca  = ca_coords[r0-1:r1] if ca_coords else []
                seg_pl  = plddts[r0-1:r1]    if plddts else []
            elif prot["multi_domain"] and len(domains) > 1:
                segs = _segment_by_plddt(plddts)
                if d_idx - 1 < len(segs):
                    r0, r1 = segs[d_idx-1]
                    res_lbl = f"{r0}–{r1} *"
                    seg_str = [(s, e) for s, e in strands_all if s >= r0 and e <= r1]
                    seg_ca  = ca_coords[r0-1:r1] if ca_coords else []
                    seg_pl  = plddts[r0-1:r1] if plddts else []
                else:
                    r0, r1  = 1, len(plddts)
                    res_lbl = f"1–{len(plddts)}"
                    seg_str, seg_ca, seg_pl = strands_all, ca_coords, plddts
            else:
                r0, r1  = 1, len(plddts)
                res_lbl = f"1–{len(plddts)}"
                seg_str, seg_ca, seg_pl = strands_all, ca_coords, plddts

            res = _classify_segment(seg_str, seg_ca, seg_pl, None, ss_pct)
            pred     = res["best_fold"]
            str_sc   = res["structural_score"]
            seq_sc   = res["sequence_score"]
            verdict  = _check(pred, expected)

            if outcome not in outcome_tbl:
                outcome_tbl[outcome] = [0, 0]

            if verdict == "VOCAB_GAP":
                tag = "VOCAB_GAP"
            elif do_abstain:
                if pred == "No confident fold assignment":
                    tag = "ABSTAIN ✓"; n_ok += 1; n_class += 1
                    outcome_tbl[outcome][0] += 1; outcome_tbl[outcome][1] += 1
                else:
                    tag = "FALSE_CALL"; n_class += 1
                    outcome_tbl[outcome][1] += 1
            elif verdict == "CORRECT":
                tag = "CORRECT ✓"; n_ok += 1; n_class += 1
                outcome_tbl[outcome][0] += 1; outcome_tbl[outcome][1] += 1
            else:
                tag = "WRONG"; n_class += 1
                outcome_tbl[outcome][1] += 1

            exp_s  = (expected[:32]  + "…") if len(expected)  > 33 else expected
            pred_s = (pred[:26] + "…") if len(pred) > 27 else pred
            print(f"  {name:<{W[0]}} {d_idx:>{W[1]}} {res_lbl:<{W[2]}} "
                  f"{exp_s:<{W[3]}} {pred_s:<{W[4]}} "
                  f"{str_sc:>{W[5]}} {seq_sc:>{W[6]}} {tag:<{W[7]}}")

        print(rule)

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    print("═" * 72)
    print("  SUMMARY")
    print("═" * 72)
    if n_class:
        acc = n_ok / n_class * 100
        print(f"\n  Per-domain fold accuracy  (PDB found, vocab covered):  "
              f"{n_ok}/{n_class} = {acc:.1f}%")
    else:
        print("\n  No classifiable domains (all PDBs missing or all vocab gaps)")

    print(f"\n  Outcome-category breakdown:")
    for cat, (nc, nt) in sorted(outcome_tbl.items()):
        bar = "█" * nc + "░" * (nt - nc)
        print(f"    {cat:<26} {nc}/{nt}  [{bar}]")

    print()
    print("  Legend:")
    print("    Str = structural score (DSSP/geometry only, 0–100)")
    print("    Seq = sequence evidence component (IPR/BLAST, 0 here — structural-only run)")
    print("    *   = boundary inferred from pLDDT drops (no explicit range in ground truth)")
    print("    VOCAB_GAP = fold class not yet in the classifier vocabulary")
    print()


if __name__ == "__main__":
    main()
