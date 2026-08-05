# Reproduce manuscript figures

Standalone Jupyter notebooks for reproducing Alps single-cell analysis figures and supplementary panels. Each notebook is self-contained (no `de_alps_utils.py`) and writes outputs under the cluster data tree at `/home/labs/vardi/amirz/pairing/Alps_All_MetaG`.

## Environment

Python 3.8+ with:

- `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `scipy`
- `statsmodels` (FDR correction in correlation / DE notebooks)
- `holoviews` / `bokeh` (Figure 2 Sankey only)

Recommended: `conda activate py38` (or equivalent env used on the Weizmann cluster).

## Data root

All paths are configured relative to:

```
/home/labs/vardi/amirz/pairing/Alps_All_MetaG
```

Shared external inputs include ViHoP Cell Ranger matrices, eggNOG annotations (`MetaG_eggnog.emapper.annotations`), and viral genome files under `genomes_of_interest/`.

## Recommended run order

Run preprocessing first; most figure notebooks depend on `Alps_preprocessed_Dinoflagellata.h5ad`.

| Step | Notebook | Main output |
|------|----------|-------------|
| 1 | `subset_dinoflagellate_data.ipynb` | `Alps_Dinoflagellata.h5ad` (422 cells) |
| 2 | `preprocess_data.ipynb` | `Alps_preprocessed_Dinoflagellata.h5ad` (372 cells) |
| 3 | `figure_s1_abundance.ipynb` | `gene_expression/figure_s1_abundance/*.png` |
| 4 | `figure_2_sankey.ipynb` | Sankey plot (ViHoP homology) |
| 5 | `figure_3_viral_genomes.ipynb` | `gene_expression/figure_3_viral_genomes/*.png` |
| 6 | `figure_4_dino_de.ipynb` | DE volcano + `differentially_expressed.Dinoflagellata.csv` |
| 7 | `figure_4_scatter.ipynb` | UMAP scatter PNGs (viral UMIs, viperin) |
| 8 | `figure_s5_correlation.ipynb` | Viperin correlation TSV + lollipop plot |

Steps 3–8 can be run in any order after step 2, except where noted below.

## Notebooks

### Data preparation

- **`subset_dinoflagellate_data.ipynb`** — Merge three Alps scRNA-seq samples and retain ViHoP-assigned Dinophyceae cells.
- **`preprocess_data.ipynb`** — QC, normalization, scaling, PCA, UMAP, and t-SNE on the dinoflagellate subset.

### Main figures

- **`figure_2_sankey.ipynb`** — Sankey diagram linking host cells and viruses from homology search results (ViHoP directory).
- **`figure_3_viral_genomes.ipynb`** — Circular pseudogenome expression maps for AIL3.217, AIL3.53, and AIL3.4.
- **`figure_4_dino_de.ipynb`** — Differential expression (infected vs uninfected) with broken-axis volcano plot.
- **`figure_4_scatter.ipynb`** — UMAP maps of % viral UMIs and normalized viperin (`140101_m_AIL3`) expression.

### Supplementary figures

- **`figure_s1_abundance.ipynb`** — Eukaryotic class abundance bar charts (18S, scRNA, scDNA).
- **`figure_s5_correlation.ipynb`** — Host–virus Spearman correlations for viperin; exports `host_virus_correlations.significant.viperin.tsv` and lollipop plot. Uses precomputed MetaG-to-genome mapping in `gene_expression/AIL3.53_METAG.tsv`.

## Usage

Open and run notebooks from this directory in Jupyter or VS Code. Each notebook has a configuration cell at the top with input/output paths; adjust `base_path` if your data live elsewhere.

```bash
cd /home/labs/vardi/amirz/github_scripts/submission/reproduce_figures
jupyter notebook
```
