# ViHoP: Virus-Host Pairing

**ViHoP** pairs giant viruses to protist hosts from environmental single-cell 10x RNA-seq data.

> Release **v1.0-Gossenkollesee** — Gossenkollesee giant-virus reference, Gossenkollesee 10x dataset reproduction bundle.

![Pipeline overview](sc_virus_host_detection_pipeline.png)

**Code availability:** https://github.com/amir-fromm/ViHoP/ (MIT License; tag `v1.0-Gossenkollesee`). Pipeline design and parameters are described in the manuscript **Methods**.

## Contents

```
ViHoP/
├── README.md
├── LICENSE
├── environment.yml          # lightweight conda env (demo + Python stages)
├── requirements.txt         # full conda lock from linux-64 run
├── run_pipeline.sh
├── pairing.config.example
├── bin/                     # pipeline scripts
├── demo/                    # small real-cell FASTQ demo (see Demo)
├── docs/
│   ├── sample_table_Gossenkollesee.tsv
│   ├── sample_table_fastq.tsv
│   └── replacement_taxa.tsv
├── config/
│   ├── pairing.config.Gossenkollesee          # manuscript Order/BLAST pairing
│   └── pairing.config.Gossenkollesee.full18s  # with 18S lastal
└── example_outputs/         # preview tables from the run
```

## 1. System requirements

### Operating system

- **Tested on:** Linux (Weizmann Wexac LSF cluster; RHEL/Rocky-class compute nodes)
- **Demo:** macOS or Linux desktop/laptop with conda + Python 3.10
- Full pipeline assumes an HPC scheduler (`bsub`); adapt `run_pipeline.sh` for SLURM/PBS if needed

### Software dependencies (versions used for Gossenkollesee)

**Python (conda):** Python 3.10 with numpy, pandas, scipy, scikit-learn, scprep, pyyaml, tqdm  
- Install from [`environment.yml`](environment.yml) (recommended for demo), or bit-for-bit Gossenkollesee lock via [`requirements.txt`](requirements.txt) (`conda create --name vihop --file requirements.txt`, linux-64).

**External binaries** (LSF modules in the original run; conda/bioconda or site modules elsewhere):

| Tool | Version used | Purpose |
|------|-------------|---------|
| CD-HIT | 4.6.6 | Viral gene deduplication |
| Cell Ranger | 10.0.0 | 10x count matrix generation |
| seqtk | — | FASTQ manipulation |
| TrimGalore | 0.6.5 | Read trimming |
| SPAdes | 3.15.3 | Per-cell RNA assembly |
| BBMap (bbduk) | — | Low-complexity filtering |
| BLAST+ | 2.11.0 | 18S rRNA search |
| LAST (lastal, lastdb) | 1574 | Viral protein & 18S alignment |
| barrnap | — | rRNA extraction from contigs |

Module load commands are set in `pairing.config`.

### Hardware

- **No non-standard hardware** (no GPU/FPGA required).
- **Demo:** normal laptop/desktop, &lt;4 GB RAM.
- **Full Gossenkollesee pipeline:** LSF cluster; Cell Ranger ~tens–hundreds of GB RAM per sample; intermediate storage can reach terabytes.

### Reference databases (full pipeline only)

| Database | Source | Used for |
|----------|--------|----------|
| Viral gene FASTA | Project-specific (e.g. Gossenkollesee metagenome genes) | Cell Ranger reference, infection detection |
| PR2 v5 | [pr2-database.org](https://pr2-database.org/) | Host 18S BLAST/lastal |
| metaPR2 v2 | [metapr2.org](https://app.metapr2.org/metapr2/) | Host 18S BLAST |
| Giant virus proteins | e.g. [GVDB](https://faylward.github.io/GVDB/) | Viral lastal |
| RefSeq + viruses (masked) | False-positive filtering | `lastal_search_viruses.sh` |
| Virus taxonomy CSV | Family-level classification | Host–virus pairing summary |

## 2. Installation guide

```bash
git clone https://github.com/amir-fromm/ViHoP.git
cd ViHoP
conda env create -f environment.yml
conda activate vihop
```

**Typical install time** on a normal desktop: about **5–15 minutes** (conda solve + package download; longer on first bioconda use or slow networks).

For a Gossenkollesee-identical linux-64 lockfile instead:

```bash
conda create --name vihop --file requirements.txt
conda activate vihop
```

## 3. Demo

A small dataset of **two real Gossenkollesee cells** (combined FASTQs + assembled transcripts from the manuscript analysis), with host annotations and viral lastal hits for the pairing stage. Sample IDs (`Alps5`, `Alps_Amir`) are the original 10x library folder names.

| Cell | Sample ID |
|------|-----------|
| `AATCGGCTCTGGATGC` | Alps5 |
| `ACATCCCTCTATGCCC` | Alps_Amir |

```bash
conda activate vihop
bash demo/run_demo.sh
```

- **What it runs:** `bin/blast_viruses_summary.py` (host–virus pairing summary; Gossenkollesee thresholds: identity ≥70%, length ≥30, `Order`, blast annotation).
- **Expected output:** tables matching `demo/expected/summary_blast_viruses.Order.tsv` and `demo/expected/virus_contigs.txt`; script prints `PASS`.
- **Expected runtime:** under **1 minute** on a normal desktop/laptop.
- FASTQs live under `demo/input/cells/` (~1.5 MB total). Details: [`demo/README.md`](demo/README.md).

This demo does **not** re-run Cell Ranger or SPAdes.

## 4. Instructions for use

### Setup for your data

```bash
cp pairing.config.example pairing.config
# Edit pairing.config: WD, database paths, FASTQ paths, queue names
# Update docs/sample_table_fastq.tsv: sample_name<TAB>fastq_path
```

For the Gossenkollesee manuscript analysis, start from `config/pairing.config.Gossenkollesee` (or `config/pairing.config.Gossenkollesee.full18s` if enabling 18S lastal).

### Run

```bash
conda activate vihop
sh run_pipeline.sh
```

Stages are toggled with `RUN_*` flags. Recommended: enable one stage at a time on first execution.

| Stage | Config flag | Key outputs |
|-------|-------------|-------------|
| 1. Database editing | `RUN_DB_EDIT` | Cell Ranger reference, `*.final.fasta`, `*.gtf` |
| 2. Cell Ranger count | `RUN_CELL_RANGER` | `{sample}/outs/filtered_feature_bc_matrix/` |
| 3. Matrix editing | `RUN_READ_MTX` | `Combined_RAW/data_raw.pickle.gz` |
| 4. Choose infected cells | `RUN_FIND_INFECTED` | `list_cells.tsv` |
| 5. Read extraction | `RUN_EXTRACT_READS` | `assemble_first_cells/all_cells/{cell}/` |
| 6. Assembly | `RUN_ASSEMBLE_CELLS` | `{cell}/assembly/transcripts.fasta` |
| 7–8. Read/transcript editing | `RUN_READ_EDITING` / `RUN_TRANSCRIPTS_EDITING` | combined FASTQs / final transcripts |
| 9–11. Host annotation | `RUN_BLAST` / `RUN_18S_LAST` / `RUN_CELL_ANNOTATION` | `blast_summary.{Level}.csv`, etc. |
| 12–13. Viral lastal + pairing | `RUN_VIRAL_LAST` / `RUN_VIRAL_SUMMARY` | `summary_blast_viruses.{Level}.tsv` |

### Reproducing manuscript Gossenkollesee results (optional)

Use `config/pairing.config.Gossenkollesee` with the three 10x libraries `Alps2`, `Alps5`, `Alps_Amir` (see `docs/sample_table_Gossenkollesee.tsv`):

- Host level `Order`; BLAST ≥98% identity, ≥300 bp; chimeric threshold 100%
- Viral lastal identity ≥70%, length ≥30; `ANNOTATION_METHOD=blast`
- Key result table: `summary_blast_viruses.Order.tsv` (preview in `example_outputs/`)

**Runtime:** full Gossenkollesee dataset (3 samples, ~10k cells) takes **days** on a busy LSF cluster. Memory: Cell Ranger was run at high RAM (order of ~288 GB per job in the original setup); SPAdes ~2 GB per cell.

`run_pipeline.sh` expects `RUN_READ_MTX` and `RUN_FIND_INFECTED` (provided in the Gossenkollesee configs above). Copy the chosen file to `pairing.config` and set `WD` before running.

## License

MIT — see [`LICENSE`](LICENSE).

## Citation

If you use ViHoP, please cite the preprint: Single-cell-resolved natural giant virus infections reveal conserved protist immunity (https://www.biorxiv.org/content/10.64898/2026.10.11.758280v1) and upstream tools (Cell Ranger, SPAdes, PR2, LAST, etc.).

## Contact

amir.fromm@weizmann.ac.il
