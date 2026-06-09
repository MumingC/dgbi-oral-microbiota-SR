# Oral microbiota in disorders of gut-brain interaction (DGBI) — 16S re-analysis pipeline

Reproducible re-analysis of seven publicly-available 16S rRNA studies of the
oral microbiota in three DGBI conditions (GERD, LPRD, IBS-D). This
repository covers both halves of the analysis:

1. **Upstream (this directory)** — raw FASTQ → ASV / genus / functional
   tables. Run on a Linux server (EC2 in our case).
2. **Downstream (`R/` subtree)** — phyloseq objects, alpha / beta diversity
   (incl. weighted UniFrac), MaAsLin2 differential abundance, metafor
   random-effects meta-analysis, figures. Run on a local workstation.

One repository, one URL — the manuscript can cite a single DOI for
end-to-end reproducibility.

## Included studies

| # | Study | Disease | BioProject | Layout | Region | Samples used |
|---|-------|---------|------------|--------|--------|--------------|
| 1 | Ziganshina 2020 | GERD  | PRJNA598080  | PE | V3-V4 | 26  (12 ctrl / 14 case) |
| 2 | Qian 2023       | GERD  | PRJNA894717  | PE | V3-V4 | 60  (30 / 30) |
| 3 | Hao 2022        | GERD  | PRJNA824804  | SE | V3-V4 | 24 oral+saliva (8 / 16) |
| 4 | Kawar 2021      | GERD  | PRJNA674379  | SE | V1-V3 | 128 (92 / 36) |
| 5 | Zheng 2024      | LPRD  | PRJNA996485  | PE | V3-V4 | 182 (102 / 80; LPRD_8W excluded) |
| 6 | Li 2025         | IBS-D | PRJNA1118066 | PE | V3-V4 | 35  (12 / 23) |
| 7 | Tang 2023       | IBS-D | PRJNA873889  | PE | V4-V5 | 62 tongue (10 / 52) |

Kim 2023 (FD) was deferred — Ion Torrent multiplexed FASTQ needs a barcode
file the authors have not yet released.

## Environment

| Tool | Version | Notes |
|------|---------|-------|
| QIIME2 amplicon distribution | `2025.10` | conda env `qiime2-amplicon-2025.10` — used for import, cutadapt, DADA2, classify-sklearn, MAFFT/FastTree, exports |
| q2-picrust2 plugin           | `2024.5`  | conda env `q2-picrust2-2024.5` — used only for the PICRUSt2 stage. 2024.5 cannot read 2025.10 archive v7.0, so the pipeline re-exports the table/rep-seqs to biom+fasta and re-imports them under 2024.5. |
| SRA Toolkit                  | `3.4.1`   | `prefetch`, `fasterq-dump` |
| EDirect                      | latest    | `esearch`, `efetch` for BioProject → run-list resolution |
| Python                       | `3.10+`   | only standard library used |

Reference databases (paths configurable via env vars):

| DB | Version | Region classifiers expected at |
|----|---------|--------------------------------|
| SILVA  | v138 nr99 (99% identity)         | `$SILVA_CLF_DIR/silva-138-99-{V3-V4,V1-V3,V4-V5}-classifier.qza` |
| eHOMD  | v16.02                           | `$EHOMD_CLF_DIR/ehomd-16.02-{V3-V4,V1-V3,V4-V5,V1-V2}-classifier.qza` |

Both classifiers are run against the same rep-seqs for every study (dual
taxonomy). SILVA is the main analysis DB; eHOMD is supplementary.

The downstream `R/` layer additionally requires R >= 4.4 with the packages
listed under [Downstream R analysis](#downstream-r-analysis).

## Repository layout

```
.
├── README.md
├── download_sra.sh         # SRA pull (prefetch + fasterq-dump per BioProject)
├── qiime2_pipeline.sh      # full QIIME2 + PICRUSt2 pipeline per study
├── study_params.tsv        # per-study primers / trunc / batching / RC flags
│                           #   — also Supplementary Table S1
└── R/                      # downstream phyloseq + MaAsLin2 + metafor
    ├── 00_setup.R                  # install/load CRAN + Bioconductor packages, project paths
    ├── write_session_info.R        # capture exact package versions
    ├── lib/                        # shared helpers + study registry
    │   ├── studies.R               #   study registry (paths, groups, regions)
    │   ├── build_phyloseq.R
    │   ├── filter_phyloseq.R
    │   ├── heatmap_helpers.R
    │   ├── menaquinol_ko_names.tsv
    │   └── metacyc_pathway_names.tsv
    ├── scripts/                    # core analysis pipeline
    │   ├── 01_build_phyloseq.R          # assemble per-study phyloseq objects
    │   ├── 02_alpha_diversity.R         # alpha diversity + Wilcoxon
    │   ├── 03_beta_diversity.R          # beta diversity / PERMANOVA (adonis2) / PERMADISP
    │   ├── 04_differential_abundance.R  # per-study MaAsLin2 (genus; SILVA & eHOMD)
    │   ├── 05_meta_analysis.R           # random-effects meta-analysis per disease
    │   ├── 06_picrust2_differential.R   # MaAsLin2 + meta-analysis on PICRUSt2 KO / pathway
    │   ├── 08_heatmap.R                 # genus effect-size heatmap
    │   ├── 09_heatmap_picrust2.R        # functional (PICRUSt2 pathway) heatmap
    │   ├── 11_fig2_tally_combined.R     # combined Figure 2 tally plot (manuscript)
    │   └── supplementary/          # sensitivity, supplementary-figure, QC scripts
    │       ├── 02b_alpha_diversity_by_disease.R   # alpha figure grouped by disease
    │       ├── 03b_beta_v3v4_gerd_pooled.R        # pooled V3-V4 GERD beta diversity
    │       ├── 05b_sensitivity_amplicon.R         # meta-analysis restricted to V3-V4 studies
    │       ├── 07_tally.R                          # supplementary forests + genus tally data
    │       ├── 08b_study_correlation.R            # per-study coefficient correlation check
    │       ├── 10_ko_menaquinol.R                 # menaquinol (vitamin K2) KO-level heatmap
    │       └── qc_silva_ehomd_crossref.R          # SILVA vs eHOMD taxonomy cross-reference
    └── results/                    # summary result tables only (see note below)
```

`R/results/` holds summary result tables only: top-level diversity,
meta-analysis, and PERMANOVA TSVs, plus per-study MaAsLin2
`significant_results.tsv` and `all_results.tsv` for genus and PICRUSt2
KO/pathway, under both SILVA and eHOMD.

Not tracked (see `.gitignore`): raw FASTQ, the cached phyloseq object, MaAsLin2
`features/`/`fits/` intermediates, run logs, and figures — all regenerable from
the scripts.

## Inputs / outputs

```
$RAW_ROOT/                                  # default: ~/sra_downloads/
├── Ziganshina_2020_GERD/  SRR*_{1,2}.fastq.gz
├── Qian_2023_GERD/        SRR*_{1,2}.fastq.gz
├── Li_2025_IBS/           SRR*_{1,2}.fastq.gz
├── Zheng_2024_LPRD/       SRR*_{1,2}.fastq.gz
├── Tang_2023_IBS/         SRR*_{1,2}.fastq.gz
├── Hao_2022_GERD/         SRR*.fastq         (uncompressed; SE)
└── Kawar_2021_GERD/       SRR*.fastq         (uncompressed; SE)

$OUT_ROOT/{DISEASE}/{Author_Year}/           # default: ~/oral_microbiota_DGBI/
├── manifest.tsv
├── {paired-end,single-end}-demux.qza
├── trimmed-demux.qza  trimmed-demux.qzv
├── table.qza  rep-seqs.qza  denoising-stats.qza
├── (split-batch only) table-long.qza  table-short.qza
├── taxonomy-silva.qza   taxonomy-ehomd.qza
├── aligned-rep-seqs.qza  masked-aligned-rep-seqs.qza
├── unrooted-tree.qza   rooted-tree.qza   rooted-tree.nwk
├── genus-table-silva.qza  genus-table-ehomd.qza
├── q2-picrust2_output/
├── exported-silva/      asv-table.tsv  genus-table.tsv  taxonomy.tsv
├── exported-ehomd/      asv-table.tsv  genus-table.tsv  taxonomy.tsv
├── exported-picrust2/   pathway_abundance.tsv  ko_metagenome.tsv
└── metadata.tsv                              # sample-id, group, ...
```

These TSV bundles plus `rooted-tree.nwk` and `metadata.tsv` are the only
files consumed by the R analysis (phyloseq + MaAsLin2 + metafor).

## Diversity is computed in R, not here

`qiime2_pipeline.sh` deliberately does NOT call `core-metrics-phylogenetic`.
All alpha- and beta-diversity numbers reported in the manuscript — Shannon,
observed features, Chao1, Bray-Curtis, weighted UniFrac, the PERMANOVA
tests — are produced by the R scripts (`R/scripts/02_alpha_diversity.R`,
`R/scripts/03_beta_diversity.R`) from the exported feature table, rooted
tree, and metadata. Keeping the QIIME2 diversity step out of the pipeline
avoids ambiguity about which numbers a reader is looking at.

## Run order

```bash
# 1) raw data — once, ~30 GB
./download_sra.sh                                  # all 7 BioProjects

# 2) classifiers — once. Train SILVA v138 nr99 region classifiers
#    (V3-V4, V1-V3, V4-V5) into $SILVA_CLF_DIR, and the eHOMD v16.02
#    classifiers into $EHOMD_CLF_DIR. (Classifier training is not in this
#    repo; standard q2-feature-classifier `fit-classifier-naive-bayes`
#    procedure on the SILVA/eHOMD reference seqs trimmed to each region.)

# 3) upstream pipeline — per study, ~hours
conda activate qiime2-amplicon-2025.10
./qiime2_pipeline.sh                               # all 7 studies
./qiime2_pipeline.sh Li_2025 Tang_2023             # subset

# 4) downstream analysis — local workstation
cd R/
Rscript 00_setup.R
Rscript scripts/01_build_phyloseq.R
# ...through scripts/11_fig2_tally_combined.R, then scripts/supplementary/<name>.R
```

The upstream stages are idempotent: each step skips if its output artefact
already exists, so partial runs can be resumed without re-doing finished
work.

### Split-batch studies (Zheng 2024, Tang 2023)

Each of these two studies contains samples from two sequencing batches with
different post-trimming read-length distributions. DADA2 requires a single
`trunc-len` per run, so each batch is denoised separately and the resulting
feature tables are merged with `qiime feature-table merge` /
`qiime feature-table merge-seqs`.

The pipeline expects `trimmed-demux-long.qza` and `trimmed-demux-short.qza`
(filtered subsets of `trimmed-demux.qza`) before invoking the DADA2 stage.
The split decision itself is informed by inspecting `trimmed-demux.qzv`'s
per-sample read-length quantiles — a manual step that varies between
datasets.

### Reverse-complement of rep-seqs (Hao 2022, Kawar 2021)

For these two studies the DADA2 rep-seqs ended up in reverse-complement
orientation. **Reverse-complementing is applied *only* before the PICRUSt2
stage**, not before taxonomy: `classify-sklearn` (Naive Bayes) classifies
both strands, so SILVA and eHOMD taxonomy is unaffected. PICRUSt2's
`hmmalign` step, however, requires forward orientation, so the pipeline
RC's the rep-seqs FASTA in-place inside `_picrust_reimport/` immediately
before re-importing for PICRUSt2 (controlled by the per-study
`rc_rep_seqs = 1` flag in `study_params.tsv`).

## Downstream R analysis

The `R/` subtree takes the exported feature tables, rooted tree, and
metadata and produces every diversity, differential-abundance,
meta-analysis, and functional result in the manuscript. All paths are
resolved with `here::here()`, so the code runs from any machine once the
input feature tables are placed under the study directories listed in
`R/lib/studies.R` — no paths need editing.

```r
# Open the R/ directory as your working directory, then:
source("00_setup.R")                       # installs + loads CRAN/Bioconductor packages
source("scripts/01_build_phyloseq.R")
# ... then run the core scripts 02 to 11 in order.
# Supplementary / sensitivity analyses: source("scripts/supplementary/<name>.R")
```

A fixed seed (`set.seed(123)`) is set in the permutation-based scripts
(`03_beta_diversity.R`, `supplementary/03b_beta_v3v4_gerd_pooled.R`) so
PERMANOVA / PERMADISP p-values are reproducible.

### Requirements

- R >= 4.4
- CRAN: tidyverse, vegan, ape, patchwork, here, broom, ggpubr, metafor
- Bioconductor: phyloseq, Maaslin2, ANCOMBC, microbiome

`00_setup.R` installs anything missing on first run.

### Downstream parameters (as reported in Methods)

- **MaAsLin2**: TSS + LOG normalisation, `min_prevalence = 0.1`,
  `min_abundance = 0.0001`, fixed effect = `group`, no covariates.
- **Meta-analysis**: `metafor::rma`, DerSimonian-Laird (DL) random-effects.
- **Multiple-testing correction**: Benjamini-Hochberg within disease group
  (not pooled across diseases).
- **LPRD** has only Zheng (k = 1) → not meta-analysed; reported descriptively.
- **Heatmap concordance**: per-study q < 0.25 marked with a star; bold box on
  meta-analysis padj < 0.05.
- **Sensitivity analysis**: V3-V4-only re-run (drops Kawar V1-V3 and Tang
  V4-V5).

## Notes on individual studies

- **Li 2025 trunc**: the initial trunc 272/271 dropped 10/35 samples; the
  successful re-run used 270/270, which is what the pipeline uses.
- **Hao 2022 platform**: Roche LS454 GS FLX Titanium, single-end (confirmed
  via ENA record for SRR18688579; Hao 2022 Int J Cancer
  [doi:10.1002/ijc.34191](https://doi.org/10.1002/ijc.34191)). The
  ~400–450 bp pyrosequencing reads have primers already trimmed at
  submission, which is why cutadapt finds 0% primer matches; the cutadapt
  step is effectively a no-op for this study.
- **Kawar 2021 read direction**: the deposited single-end FASTQs are the
  *reverse* read of the V1-V3 amplicon — the original paper's Methods
  state "reverse sequences from the FASTQ files were analyzed using
  QIIME2" (Kawar 2021 *Sci Rep*). The 27F primer therefore sits at the 3'
  end of the read, not the 5'.
- **Kawar 2021 5' primer (empirical check)**: across 1500 reads from 3
  accessions:
  - 519R `ATTACCGCGGCTGCTGG` (and degenerate variants
    `GWATTACCGCGGCKGCTG`, `GTATTACCGCGGCTGCTG`): 0/1500 anywhere.
  - 27F `AGAGTTTGATCMTGGCTCAG`: 0/1500.
  - The two dominant 5' 14-mers (`CACGTAGTTAGCCG` and `CACGGAATTAGCCG`,
    together ~75% of reads) reverse-complement to `CGGCTAACTACGTG` /
    `CGGCTAATTCCGTG`, i.e. ≈ E. coli 16S position 537 — immediately
    downstream of where 519R binds.
  - Conclusion: 519R was already trimmed at submission. The pipeline
    therefore sets `fwd_primer=""` for Kawar and skips cutadapt, then
    applies `rc_rep_seqs=1` so PICRUSt2 sees forward-oriented sequences.
- **Kawar 2021 truncation length**: our pipeline uses 220 nt, deliberately
  shorter than the paper's 258 nt (after Q25 trimming). This is a
  re-analysis choice driven by our own quality-curve inspection, not a
  transcription slip from the paper.

All four notes are also called out inline in `qiime2_pipeline.sh` near the
relevant SPEC entries and in `study_params.tsv`.

## Provenance — what's verbatim vs reconstructed

| Item | Source | Confidence |
|------|--------|-----------|
| BioProject IDs + the prefetch / fasterq-dump loop | original `~/.bash_history` line 654-655 | recovered verbatim |
| Ziganshina primers | original `~/.bash_history` line 623 | recovered verbatim |
| Final DADA2 trunc values for all 7 studies | EC2 `pipeline_log.md` (the "DONE" entries) plus the 2026-04-14 batch-comparison entries for Zheng/Tang | recovered |
| Other six studies' primers | from original-paper Methods sections (compiled into the per-study reference table on EC2) — not present in bash_history | reconstructed from sources |
| PICRUSt2 invocation + re-import workflow | EC2 `run_picrust2.sh` (preserved as-is) | recovered verbatim |
| RC step for Hao / Kawar | EC2 `pipeline_log.md` "PICRUSt2 RC issue" note + `handoff.md` Gotchas | recovered |
| Split-batch merge for Zheng / Tang | EC2 `pipeline_log.md` 2026-04-14, 2026-04-15 entries | recovered |
| Hao 2022 platform = LS454 | ENA run record for SRR18688579 | verified post-hoc |
| Kawar 2021 read direction = reverse, 519R absent | empirical scan of 1500 reads (see note above) | verified post-hoc |

## Reproducibility

Permutation-based tests use a fixed seed (`set.seed(123)`). To record exact
package versions, after `source("00_setup.R")` run
`source("write_session_info.R")`, which writes `reproducibility/sessionInfo.txt`.
For stricter pinning, initialise [`renv`](https://rstudio.github.io/renv/) in
the `R/` directory (`renv::init(); renv::snapshot()`) and commit the
`renv.lock`.

## Release & DOI

The repository is versioned with git tags. Each manuscript-relevant version
is archived to Zenodo for a citable DOI:

```bash
# tag the version submitted to the journal
git tag -a v1.0.0 -m "Submission to <journal>: <date>"
git push origin v1.0.0

# Zenodo (via the GitHub integration) picks up the tag, builds an archive,
# and mints a DOI. Replace [repository URL] placeholders in the manuscript
# and Supplementary Information with that DOI (form: 10.5281/zenodo.NNNNN).
```

If the GitHub ↔ Zenodo integration isn't enabled on the account, do it
manually at https://zenodo.org/account/settings/github/ before pushing the
tag. The DOI is minted on the first release and a new sub-DOI is minted on
each subsequent tag; cite the version-specific DOI in the manuscript so
reviewers always land on exactly the code that produced the figures.

## Citing / licensing

Code is released under the MIT License (see `LICENSE`). If you use it, please
cite the associated article and this repository (see `CITATION.cff`).

The processed outputs in this repository (ASV / genus / functional tables,
rooted trees) are derived from public SRA submissions; reuse is bound by
the licence terms of the original studies. SILVA and eHOMD reference
databases retain their own licences. PICRUSt2 predictions should be
treated as hypothesis-generating, especially in oral microbiome contexts
where reference coverage is sparser than gut.
