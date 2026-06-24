---
title: Dark Proteome Pipeline
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
pinned: false
---

# Dark Proteome Pipeline

A bioinformatics pipeline for functional annotation of uncharacterised ("dark") proteins, built as part of a BSc dissertation project. The pipeline submits protein sequences and structures to multiple remote annotation services in parallel and saves all raw results to a named output folder per protein.

## Project structure

```
dark-proteome-pipeline/
├── scripts/
│   ├── run_pipeline.py        # Master script — runs all tools in parallel
│   ├── run_interproscan.py    # EBI InterProScan (domain/family annotation)
│   ├── run_blast.py           # NCBI BLASTp against nr
│   ├── run_phobius.py         # EBI Phobius (signal peptide + TM topology)
│   ├── run_hmmer.py           # EBI HMMER hmmscan against Pfam
│   ├── run_foldseek.py        # FoldSeek structural similarity search
│   ├── run_signalp.py         # DTU SignalP 6.0 (script ready; server intermittent)
│   ├── run_hhpred.py          # MPI HHpred structural homology (script ready; server intermittent)
│   └── parse_interproscan.py  # Parse InterProScan JSON into a flat TSV table
├── data/
│   ├── input/                 # Input FASTA and PDB files
│   └── output/                # Results, organised by protein name
├── logs/
├── .env                       # Local environment variables (not committed)
├── .gitignore
├── requirements.txt
└── README.md
```

## Setup

**1. Clone the repository**
```bash
git clone https://github.com/<your-username>/dark-proteome-pipeline.git
cd dark-proteome-pipeline
```

**2. Install dependencies**
```bash
pip install -r requirements.txt
```

**3. Configure your email**

EBI and NCBI fair-use policies require a valid email address with each job submission. Create a `.env` file in the project root:
```
EBI_EMAIL=your@email.com
```

## Usage

### Master pipeline (recommended)

Run all applicable tools for a protein in a single command. All tools submit and poll simultaneously, so total runtime equals the slowest tool rather than the sum of all.

```bash
# FASTA + PDB — runs InterProScan, BLASTp, Phobius, HMMER, and FoldSeek
python scripts/run_pipeline.py --fasta data/input/PROTEIN.fasta --pdb data/input/PROTEIN.pdb

# FASTA only — runs InterProScan, BLASTp, Phobius, HMMER
python scripts/run_pipeline.py --fasta data/input/PROTEIN.fasta

# PDB only — runs FoldSeek
python scripts/run_pipeline.py --pdb data/input/PROTEIN.pdb

# Custom output folder name
python scripts/run_pipeline.py --fasta data/input/PROTEIN.fasta --pdb data/input/PROTEIN.pdb --name MyProtein
```

Results are saved to `data/output/<protein_name>/` with one file per tool.

### Individual tools

Each tool can also be run standalone:

```bash
python scripts/run_interproscan.py data/input/PROTEIN.fasta
python scripts/run_blast.py        data/input/PROTEIN.fasta
python scripts/run_phobius.py      data/input/PROTEIN.fasta
python scripts/run_hmmer.py        data/input/PROTEIN.fasta
python scripts/run_foldseek.py     data/input/PROTEIN.pdb
```

### Parse InterProScan results

```bash
python scripts/parse_interproscan.py data/output/<result>.json
```

Writes a TSV file with one row per match location and the following columns:

| Column | Description |
|---|---|
| `protein_id` | Sequence identifier from the FASTA file |
| `database` | Source database (e.g. PFAM, PANTHER, GENE3D) |
| `signature_accession` | Database-specific accession |
| `signature_name` | Short name of the signature |
| `signature_description` | Description of the signature |
| `type` | Match type (DOMAIN, FAMILY, REPEAT, etc.) |
| `interpro_accession` | InterPro entry accession (if available) |
| `interpro_description` | InterPro entry description |
| `start` | Match start position (1-based) |
| `end` | Match end position (1-based) |
| `evalue` | E-value of the match (where available) |
| `go_terms` | GO term annotations (pipe-separated, where available) |

## Tools and output formats

| Tool | Input | Databases | Output format |
|---|---|---|---|
| InterProScan | FASTA | Pfam, PANTHER, GENE3D, SUPERFAMILY, and more | JSON |
| BLASTp | FASTA | NCBI nr | JSON |
| Phobius | FASTA | — | TSV (text) |
| HMMER hmmscan | FASTA | Pfam | Text (standard HMMER format) |
| FoldSeek | PDB | afdb50, afdb-swissprot, afdb-proteome, BFVD | tar.gz of M8 alignment files |
| SignalP 6.0 | FASTA | — | Text (DTU server intermittently unavailable) |
| HHpred | FASTA | PDB_mmCIF70 | .hhr text (MPI server intermittently unavailable) |

## Fair-use guidelines

| Service | Guidelines |
|---|---|
| EBI (InterProScan, Phobius, HMMER) | Max 30 simultaneous jobs; poll ≥ every 3 s; include email |
| NCBI BLAST | Max 1 request per 10 s; include email and tool name |
| FoldSeek | Poll ≥ every 10 s; server returns RATELIMIT if exceeded |
| DTU (SignalP) | Academic use only; do not automate excessively |

Full EBI guidelines: https://www.ebi.ac.uk/Tools/webservices/help/faq

## Requirements

- Python 3.9+
- `requests`
- `python-dotenv`
