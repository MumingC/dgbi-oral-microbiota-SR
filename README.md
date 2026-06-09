# Oral microbiota in DGBI: analysis code

Analysis code and summary result tables for a systematic review and re-analysis
of the oral microbiota in disorders of gut-brain interaction (DGBI), covering
gastro-oesophageal reflux disease (GERD), laryngopharyngeal reflux disease
(LPRD), and diarrhoea-predominant irritable bowel syndrome (IBS-D).

Raw 16S rRNA amplicon data from seven published studies were retrieved from
public repositories, reprocessed through a uniform DADA2/QIIME 2 pipeline, and
re-analysed in R for alpha/beta diversity, differential abundance (MaAsLin2),
random-effects meta-analysis, and predicted functional potential (PICRUSt2).

## Studies re-analysed

| Study | Disease | 16S region | n samples | Layout |
|-------|---------|-----------|-----------|--------|
| Ziganshina 2020 | GERD | V3–V4 | 26 | PE |
| Qian 2023 | GERD | V3–V4 | 60 | PE |
| Hao 2022 | GERD | V3–V4 | 24 | SE |
| Kawar 2021 | GERD | V1–V3 | 128 | SE |
| Zheng 2024 | LPRD | V3–V4 | 182 | PE |
| Li 2025 | IBS-D | V3–V4 | 35 | PE |
| Tang 2023 | IBS-D | V4–V5 | 62 | PE |

## Pipeline overview

Raw FASTQ are retrieved with `download_sra.sh` (NCBI SRA / DDBJ), then denoised
in QIIME 2 / DADA2 to ASV feature tables and taxonomy (run upstream; see the
Reproducibility section). The R code in this repository takes those feature
tables and runs the re-analysis. Paths are anchored with `here::here()`.

Core scripts (`R/scripts/`):

```
01_build_phyloseq.R          Assemble per-study phyloseq objects
02_alpha_diversity.R         Alpha diversity + Wilcoxon
03_beta_diversity.R          Beta diversity / PERMANOVA (adonis2) / PERMADISP
04_differential_abundance.R  Per-study MaAsLin2 (genus; SILVA & eHOMD)
05_meta_analysis.R           Random-effects meta-analysis per disease
06_picrust2_differential.R   MaAsLin2 + meta-analysis on PICRUSt2 KO / pathway
08_heatmap.R                 Genus effect-size heatmap
09_heatmap_picrust2.R        Functional (PICRUSt2 pathway) heatmap
11_fig2_tally_combined.R     Combined Figure 2 tally plot (manuscript)
```

Supplementary, sensitivity, and QC scripts (`R/scripts/supplementary/`):

```
02b_alpha_diversity_by_disease.R   Alpha figure grouped by disease
03b_beta_v3v4_gerd_pooled.R        Pooled V3–V4 GERD beta diversity
05b_sensitivity_amplicon.R         Meta-analysis restricted to V3–V4 studies
07_tally.R                         Supplementary forests + genus tally data
08b_study_correlation.R            Per-study coefficient correlation check
10_ko_menaquinol.R                 Menaquinol (vitamin K2) KO-level heatmap
qc_silva_ehomd_crossref.R          SILVA vs eHOMD taxonomy cross-reference
```

Helper functions live in `R/lib/`: `studies.R` (study registry),
`build_phyloseq.R`, `filter_phyloseq.R`, and `heatmap_helpers.R`.

A fixed seed (`set.seed(123)`) is set in the permutation-based scripts (`03`,
`03b`) so PERMANOVA / PERMADISP p-values are reproducible.

## How to run

```r
# Open the R/ directory as your working directory, then:
source("00_setup.R")     # installs + loads CRAN/Bioconductor packages
source("scripts/01_build_phyloseq.R")
# ... then run the core scripts 02 to 11 in order.
# Supplementary / sensitivity analyses: source("scripts/supplementary/<name>.R")
```

All paths are resolved with `here::here()`, so the code runs from any machine
once the input feature tables are placed under the study directories listed in
`R/lib/studies.R`. No paths need editing.

### Requirements
- R >= 4.4
- CRAN: tidyverse, vegan, ape, patchwork, here, broom, ggpubr, metafor
- Bioconductor: phyloseq, Maaslin2, ANCOMBC, microbiome

`00_setup.R` installs anything missing on first run.

## What's in this repository

- `R/00_setup.R`: package install/load and project paths.
- `R/lib/`: shared helper functions and the study registry.
- `R/scripts/`: core analysis pipeline (01 to 11).
- `R/scripts/supplementary/`: sensitivity, supplementary-figure, and QC scripts.
- `R/results/`: summary result tables only (top-level diversity, meta-analysis,
  and PERMANOVA TSVs, plus per-study MaAsLin2 `significant_results.tsv` and
  `all_results.tsv` for genus and PICRUSt2 KO/pathway, SILVA & eHOMD).
- `R/write_session_info.R`: helper to capture exact package versions (see below).
- `download_sra.sh`: raw-data retrieval script (edit `OUTPUT_DIR` before running).

Not tracked (see `.gitignore`): raw FASTQ, the cached phyloseq object, MaAsLin2
`features/` and `fits/` intermediates, run logs, and figures, all of which are
regenerable from the scripts.

## Reproducibility

Upstream denoising (raw FASTQ to ASV feature tables) was performed in QIIME 2 /
DADA2 outside this repository; per-study truncation parameters followed each
run's quality profiles. This repository covers the R re-analysis layer that
takes the resulting feature tables and produces the diversity, differential-
abundance, meta-analysis, and functional results reported in the manuscript.
Raw data are publicly available via the accessions handled in `download_sra.sh`.

Permutation-based tests use a fixed seed (`set.seed(123)`). To record exact
package versions, after `source("00_setup.R")` run
`source("write_session_info.R")`, which writes `reproducibility/sessionInfo.txt`.
For stricter pinning, initialise [`renv`](https://rstudio.github.io/renv/) in the
`R/` directory (`renv::init(); renv::snapshot()`) and commit the `renv.lock`.

## License & citation

Code is released under the MIT License (see `LICENSE`). If you use it, please
cite the associated article and this repository (see `CITATION.cff`).
