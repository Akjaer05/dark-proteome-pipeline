# Dark Proteome Pipeline

A bioinformatics pipeline for functional annotation of uncharacterised ("dark") proteins, built as part of a BSc dissertation project. The pipeline submits protein sequences to the [EBI InterProScan REST API](https://www.ebi.ac.uk/Tools/services/rest/iprscan5) and parses the results into summary tables for downstream analysis.

## Project structure

```
dark-proteome-pipeline/
├── scripts/
│   ├── run_interproscan.py    # Submit a FASTA file to InterProScan, save raw JSON
│   └── parse_interproscan.py  # Parse JSON result into a flat TSV summary table
├── data/
│   ├── input/                 # Input FASTA files (one protein per file)
│   └── output/                # JSON results and TSV summary tables
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

EBI's fair-use policy requires a valid email address to be sent with each job submission. Create a `.env` file in the project root:
```
EBI_EMAIL=your@email.com
```

## Usage

### 1. Submit a FASTA file to InterProScan

```bash
python scripts/run_interproscan.py data/input/your_protein.fasta
```

- Submits the sequence to EBI InterProScan
- Polls for job completion (every 10 seconds, per EBI fair-use guidelines)
- Saves the raw JSON result to `data/output/`

### 2. Parse the JSON result into a summary table

```bash
python scripts/parse_interproscan.py data/output/<result>.json
```

Writes a TSV file to `data/output/` with one row per match location and the following columns:

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

## EBI fair-use guidelines

- Submit no more than 30 jobs simultaneously
- Poll no more than once every 3 seconds (this pipeline uses 10 s)
- Always include a valid email address
- Full guidelines: https://www.ebi.ac.uk/Tools/webservices/help/faq

## Requirements

- Python 3.9+
- `requests`
- `python-dotenv`
