# Huge-ref Cell Ranger (Host_Virus_CR)

Run Cell Ranger count against the ~33 GB `Host_Virus_CR` transcriptome from
`MetaT_Eilat_261007` using the Snakemake ViHoP pipeline.

## Why special settings

- STAR mmaps ~30 GB; each ALIGN_AND_COUNT chunk RSS is ~42–46 GB.
- Without `--localmem`, Cell Ranger plans from host `MemTotal` and ignores the LSF cgroup.
- Martian’s default chunk cap (~37 GB) kills STAR; use `--overrides` with the **full** stage id
  (see `config/cellranger_align_mem_huge.json`). A short key like `"ALIGN_AND_COUNT"` is not enough.
- `ulimit -v unlimited` must run **inside** the job shell (the wrapper does this).
- LSF `rusage[mem=…]` is **per slot**. The Snakemake LSF plugin sets
  `rusage = resources.mem_mb / threads` when `mem_mb` is given as **total** job memory.
  MetaT config: `threads: 8`, `mem_mb: 160000` → 20 GB/slot × 8 = **160 GB** total,
  with `--localcores=8 --localmem=160 --localvmem=240`.

## Config

| File | Role |
|------|------|
| `config/config_metat_eilat.yaml` | Samples, Host_Virus_CR paths, long queue, 160 GB, CR 10.1.0 |
| `config/cellranger_align_mem_huge.json` | Full-stage Martian mem override |

Pilot sample is `Eilat_04Mai` only. After the ALIGN gate passes, add
`Eilat_10X_1`, `Eilat_10X_3`, `Eilat_10X_8` to `input.samples` and set
`force_clean: false` for resume-friendly restarts.

## Separate workdir (Alps vs Eilat)

Do **not** write MetaT/Eilat outputs into the repo tree (or an Alps project). Use `--workdir`:

```text
repo (code only):     /home/labs/vardi/amirz/snakemake_vihop/ViHoP
Eilat run directory:  /home/labs/vardi/amirz/pairing/MetaT_Eilat_261007/snakemake_vihop
```

`results/`, `logs/`, and `.snakemake/` land in the workdir. Cell Ranger outs
default to `<workdir>/results/cellranger/<sample>/outs/…`. Override with
`--bam-root /path/to/parent` (parent of `<sample>/outs`).

## Commands (no submit until ready)

On a Wexac login node, with the `vihop` conda env active:

```bash
REPO=/home/labs/vardi/amirz/snakemake_vihop/ViHoP
WD=/home/labs/vardi/amirz/pairing/MetaT_Eilat_261007/snakemake_vihop

# 1) Dry-run: expect one cellranger_count job, Host_Virus_CR paths, no mkref
$REPO/vihop --workdir "$WD" --configfile config/config_metat_eilat.yaml -n

# 2) Inspect resources (optional): Total RAM ≈ 160 GB, queue=long, -n 8
$REPO/vihop --workdir "$WD" --configfile config/config_metat_eilat.yaml -n --printshellcmds

# 3) Pilot ONLY (when you explicitly ask to submit)
$REPO/vihop --workdir "$WD" --configfile config/config_metat_eilat.yaml --jobs 1
```

Do **not** use profile default `jobs: 50` for huge-ref; always pass `--jobs 1` (pilot)
or at most `--jobs 2` after the first sample clears ALIGN.

## Validation gates

1. **Dry-run** shows `Eilat_04Mai` → `Host_Virus_CR` and does not schedule mkref/dedup.
2. **bsub line**: `-q long`, `-n 8`, `rusage[mem=20000]` (or equivalent totaling ~160 GB).
3. After ALIGN_AND_COUNT starts, check a chunk `_jobinfo` or monitor line:
   **`allowed` must be > 50G** (not 37G). If still 37G, stop — fix overrides / wipe pipestance
   with `force_clean: true` and re-run the pilot only.
4. RSS should stay under ~160 GB (shared STAR index + a few live chunks).
5. Only then add the other three samples (`--jobs 1` or `2`).

## Resume vs wipe

- `force_clean: true` (MetaT pilot default): deletes `results/cellranger/<sample>` before count
  so an old 37G cap stored in the pipestance cannot stick.
- After a successful override is confirmed, set `force_clean: false` so interrupted jobs
  reattach. The wrapper removes a stale `_lock` only if no `cellranger` process is running.

## choose_cells → candidate_cells.tsv

After Cell Ranger finishes, flip stages to:

```yaml
stages:
  cellranger: false
  choose_cells: true
  downstream: false
```

Default (`min_umi` / `min_genes: null`) writes **all** barcodes from each
`filtered_feature_bc_matrix` as `barcode<TAB>sample` into `results/candidate_cells.tsv`.
Set integer `min_umi` / `min_genes` to require minimum UMIs or genes per cell.

```bash
$REPO/vihop --workdir "$WD" --configfile config/config_metat_eilat.yaml --jobs 1
```

Then enable `downstream` (candidate file must already exist for extraction rule gen).

## Downstream

BAMs are under `results/cellranger/<sample>/outs/` by default (`bam_root: null`).
Huge BAMs need ample `disk_mb`.
