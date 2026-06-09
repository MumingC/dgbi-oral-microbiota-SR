#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# qiime2_pipeline.sh — clean re-implementation of the 7-study QIIME2 pipeline
#
# Stages (per study):
#   1. Build manifest from $RAW_ROOT/<STUDY_DIR>/
#   2. Import as PE or SE artifact
#   3. Cutadapt primer trimming
#   4. DADA2 denoise (per-study trunc lengths; batched for Zheng & Tang)
#   5. Merge batches where applicable (Zheng, Tang)
#   6. Optional reverse-complement of rep-seqs (Hao, Kawar)
#   7. classify-sklearn against SILVA v138 (99%) AND eHOMD v16.02
#   8. MAFFT + FastTree2 rooted tree
#   9. core-metrics-phylogenetic at per-study sampling depth
#  10. q2-picrust2 full-pipeline (NSTI > 2.0 excluded)
#  11. Genus-collapse + TSV exports (MaAsLin2 / phyloseq input)
#
# Env requirements:
#   - QIIME2 amplicon distribution 2025.10 (`qiime2-amplicon-2025.10` env)
#     used for steps 1-9 and 11
#   - q2-picrust2 v2024.5 (`q2-picrust2-2024.5` env) used for step 10
#     (2024.5 cannot read 2025.10 archive v7.0 — re-export + re-import)
#   - SILVA v138 99% region classifiers under $SILVA_CLF_DIR
#   - eHOMD v16.02 region classifiers under $EHOMD_CLF_DIR
#
# Usage:
#   conda activate qiime2-amplicon-2025.10
#   ./qiime2_pipeline.sh                   # all 7 studies
#   ./qiime2_pipeline.sh Li_2025 Tang_2023 # subset (by directory basename)
#
# Skip flags (env vars):
#   SKIP_PICRUST2=1   skip q2-picrust2 stage
#
# Note: alpha/beta diversity (Shannon, Chao1, Bray-Curtis, weighted UniFrac)
# is NOT computed here. All reported diversity in the manuscript is produced
# by the R scripts (R/scripts/02_alpha_diversity.R,
# R/scripts/03_beta_diversity.R) from the exported feature table, rooted
# tree, and metadata. Keeping the QIIME2 core-metrics step out of this
# script avoids any confusion about which numbers are reported.
# ---------------------------------------------------------------------------

set -euo pipefail

# --- paths -----------------------------------------------------------------
RAW_ROOT="${RAW_ROOT:-$HOME/sra_downloads}"
OUT_ROOT="${OUT_ROOT:-$HOME/oral_microbiota_DGBI}"
SILVA_CLF_DIR="${SILVA_CLF_DIR:-$OUT_ROOT/classifiers}"
EHOMD_CLF_DIR="${EHOMD_CLF_DIR:-$HOME/databases/classifiers}"
THREADS="${THREADS:-4}"
QIIME_ENV="${QIIME_ENV:-qiime2-amplicon-2025.10}"
PICRUST_ENV="${PICRUST_ENV:-q2-picrust2-2024.5}"

# --- per-study registry ----------------------------------------------------
# This bash associative-array registry is kept in sync with the standalone
# table `study_params.tsv` (same directory). The TSV is the human-readable /
# manuscript-supplementary version; this array is what the script actually
# reads. If you edit one, edit the other.
#
# Each row is a study spec; arrays are queried by study name.
#
# Fields used:
#   DISEASE         GERD | LPRD | IBS-D
#   RAW_DIR         subdir of $RAW_ROOT
#   OUT_SUBDIR      subdir of $OUT_ROOT (e.g. GERD/Ziganshina_2020)
#   LAYOUT          PE | SE
#   REGION          V3-V4 | V1-V3 | V4-V5  (selects classifier file)
#   FWD_PRIMER      cutadapt -p-front-f / -p-front
#   REV_PRIMER      cutadapt -p-front-r (PE only; ignored for SE)
#   TRUNC_F         DADA2 --p-trunc-len-f (PE) or --p-trunc-len (SE)
#   TRUNC_R         DADA2 --p-trunc-len-r (PE only)
#   BATCH_MODE      none | split (split: a second short batch is also denoised
#                                   then merged with the main one)
#   BATCH2_TRUNC_F  short-batch trunc-f (split only)
#   BATCH2_TRUNC_R  short-batch trunc-r (split only)
#   RC_REP_SEQS     0 | 1   reverse-complement rep-seqs before PICRUSt2
#                          (sklearn taxonomy is orientation-agnostic; PICRUSt2
#                          hmmalign requires forward orientation)
#
# Primer references: Ziganshina recovered from bash_history line 623; all
# others sourced from the per-paper CLAUDE.md table. Trunc values verified
# against pipeline_log.md final-run entries (not the initial failed values
# that the older CLAUDE.md table still shows for Li 2025).

declare -A SPEC=(
  # -------- Ziganshina 2020 (GERD, PE V3-V4) ------------------------------
  [Ziganshina_2020.DISEASE]="GERD"
  [Ziganshina_2020.RAW_DIR]="Ziganshina_2020_GERD"
  [Ziganshina_2020.OUT_SUBDIR]="GERD/Ziganshina_2020"
  [Ziganshina_2020.LAYOUT]="PE"
  [Ziganshina_2020.REGION]="V3-V4"
  [Ziganshina_2020.FWD_PRIMER]="CCTACGGGNGGCWGCAG"
  [Ziganshina_2020.REV_PRIMER]="GACTACHVGGGTATCTAATCC"
  [Ziganshina_2020.TRUNC_F]="230"
  [Ziganshina_2020.TRUNC_R]="220"
  [Ziganshina_2020.BATCH_MODE]="none"
  [Ziganshina_2020.RC_REP_SEQS]="0"

  # -------- Qian 2023 (GERD, PE V3-V4) ------------------------------------
  [Qian_2023.DISEASE]="GERD"
  [Qian_2023.RAW_DIR]="Qian_2023_GERD"
  [Qian_2023.OUT_SUBDIR]="GERD/Qian_2023"
  [Qian_2023.LAYOUT]="PE"
  [Qian_2023.REGION]="V3-V4"
  [Qian_2023.FWD_PRIMER]="CCTAYGGGRBGCASCAG"
  [Qian_2023.REV_PRIMER]="GGACTACHVGGGTWTCTAAT"
  [Qian_2023.TRUNC_F]="227"
  [Qian_2023.TRUNC_R]="223"
  [Qian_2023.BATCH_MODE]="none"
  [Qian_2023.RC_REP_SEQS]="0"

  # -------- Li 2025 (IBS-D, PE V3-V4) -------------------------------------
  # NOTE: CLAUDE.md still shows the *failed* initial trunc 272/271. The final
  # successful run was 270/270 (pipeline_log.md 2026-04-14).
  [Li_2025.DISEASE]="IBS-D"
  [Li_2025.RAW_DIR]="Li_2025_IBS"
  [Li_2025.OUT_SUBDIR]="IBS/Li_2025"
  [Li_2025.LAYOUT]="PE"
  [Li_2025.REGION]="V3-V4"
  [Li_2025.FWD_PRIMER]="ACTCCTACGGGAGGCAGCAG"
  [Li_2025.REV_PRIMER]="GGACTACHVGGGTWTCTAAT"
  [Li_2025.TRUNC_F]="270"
  [Li_2025.TRUNC_R]="270"
  [Li_2025.BATCH_MODE]="none"
  [Li_2025.RC_REP_SEQS]="0"

  # -------- Zheng 2024 (LPRD, PE V3-V4, two seq batches) ------------------
  [Zheng_2024.DISEASE]="LPRD"
  [Zheng_2024.RAW_DIR]="Zheng_2024_LPRD"
  [Zheng_2024.OUT_SUBDIR]="GERD/Zheng_2024"
  [Zheng_2024.LAYOUT]="PE"
  [Zheng_2024.REGION]="V3-V4"
  [Zheng_2024.FWD_PRIMER]="CCTACGGRRBGCASCAGKVRVGAAT"
  [Zheng_2024.REV_PRIMER]="GGACTACNVGGGTWTCTAATCC"
  [Zheng_2024.TRUNC_F]="219"      # long batch
  [Zheng_2024.TRUNC_R]="224"
  [Zheng_2024.BATCH_MODE]="split"
  [Zheng_2024.BATCH2_TRUNC_F]="212"
  [Zheng_2024.BATCH2_TRUNC_R]="220"
  [Zheng_2024.RC_REP_SEQS]="0"

  # -------- Hao 2022 (GERD, SE V3-V4) -------------------------------------
  # Platform: Roche LS454 GS FLX Titanium (confirmed via ENA run record for
  # SRR18688579; ref: Hao 2022 Int J Cancer doi:10.1002/ijc.34191). Reads
  # are long single-end pyrosequencing (~400-450 bp) with primers already
  # trimmed at submission, which is why cutadapt finds 0% primer matches.
  # The FWD_PRIMER below is the expected V3-V4 forward primer; in this run
  # it acts as a no-op against the deposited reads, which is expected.
  [Hao_2022.DISEASE]="GERD"
  [Hao_2022.RAW_DIR]="Hao_2022_GERD"
  [Hao_2022.OUT_SUBDIR]="GERD/Hao_2022"
  [Hao_2022.LAYOUT]="SE"
  [Hao_2022.REGION]="V3-V4"
  [Hao_2022.FWD_PRIMER]="GGAGGCAGCAGTRRGGAAT"
  [Hao_2022.REV_PRIMER]=""
  [Hao_2022.TRUNC_F]="400"
  [Hao_2022.TRUNC_R]=""
  [Hao_2022.BATCH_MODE]="none"
  [Hao_2022.RC_REP_SEQS]="1"      # forward-orient before PICRUSt2 hmmalign

  # -------- Kawar 2021 (GERD, SE V1-V3) -----------------------------------
  # Platform: Illumina MiSeq, single-end (ENA SRR13021985; Kawar 2021 Sci Rep).
  # The deposited reads are the REVERSE direction of the V1-V3 amplicon —
  # the paper's Methods states "reverse sequences from the FASTQ files
  # were analyzed using QIIME2". The 27F primer therefore sits at the
  # 3' end of the read; the 519R primer, if present at all, would be at
  # the 5'.
  #
  # Empirical primer check (1500 reads, 3 accessions):
  #   - 519R ATTACCGCGGCTGCTGG (and degenerate variants GWATTACCGCGGCKGCTG,
  #     GTATTACCGCGGCTGCTG):  0/1500 anywhere in the reads.
  #   - 27F  AGAGTTTGATCMTGGCTCAG:                                  0/1500.
  #   - Two dominant 5' 14-mers cover ~75% of reads:
  #       CACGTAGTTAGCCG  -> RC: CGGCTAACTACGTG  (E. coli 16S pos ~537)
  #       CACGGAATTAGCCG  -> RC: CGGCTAATTCCGTG  (same region, alt phylum)
  #     i.e. reads start immediately AFTER where 519R would have bound;
  #     the primer was trimmed at submission, not deposited in the FASTQ.
  #
  # FWD_PRIMER is therefore empty — cutadapt is skipped because there is
  # no 5' primer to remove. RC_REP_SEQS=1 puts rep-seqs in forward
  # orientation before PICRUSt2 hmmalign (sklearn taxonomy is orientation-
  # agnostic, so this only matters for PICRUSt2).
  #
  # Truncation: 220 nt is our quality-derived choice, deliberately shorter
  # than the paper's 258 nt (after Q25 trimming) — a re-analysis decision,
  # not a transcription error.
  [Kawar_2021.DISEASE]="GERD"
  [Kawar_2021.RAW_DIR]="Kawar_2021_GERD"
  [Kawar_2021.OUT_SUBDIR]="GERD/Kawar_2021"
  [Kawar_2021.LAYOUT]="SE"
  [Kawar_2021.REGION]="V1-V3"
  [Kawar_2021.FWD_PRIMER]=""
  [Kawar_2021.REV_PRIMER]=""
  [Kawar_2021.TRUNC_F]="220"
  [Kawar_2021.TRUNC_R]=""
  [Kawar_2021.BATCH_MODE]="none"
  [Kawar_2021.RC_REP_SEQS]="1"

  # -------- Tang 2023 (IBS-D, PE V4-V5, two seq batches) ------------------
  [Tang_2023.DISEASE]="IBS-D"
  [Tang_2023.RAW_DIR]="Tang_2023_IBS"
  [Tang_2023.OUT_SUBDIR]="IBS/Tang_2023"
  [Tang_2023.LAYOUT]="PE"
  [Tang_2023.REGION]="V4-V5"
  [Tang_2023.FWD_PRIMER]="GTGCCAGCMGCCGCGGTAA"
  [Tang_2023.REV_PRIMER]="CCGTCAATTCMTTTGAGTTT"
  [Tang_2023.TRUNC_F]="274"
  [Tang_2023.TRUNC_R]="245"
  [Tang_2023.BATCH_MODE]="split"
  [Tang_2023.BATCH2_TRUNC_F]="224"
  [Tang_2023.BATCH2_TRUNC_R]="222"
  [Tang_2023.RC_REP_SEQS]="0"
)

DEFAULT_STUDIES=(Ziganshina_2020 Qian_2023 Li_2025 Zheng_2024 Hao_2022 Kawar_2021 Tang_2023)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
spec() { echo "${SPEC[$1.$2]:-}"; }

silva_classifier() {
  echo "$SILVA_CLF_DIR/silva-138-99-${1}-classifier.qza"
}
ehomd_classifier() {
  echo "$EHOMD_CLF_DIR/ehomd-16.02-${1}-classifier.qza"
}

write_pe_manifest() {
  local raw_dir="$1" out_file="$2"
  python3 - "$raw_dir" "$out_file" <<'PY'
import glob, os, sys
raw, out = sys.argv[1], sys.argv[2]
fwd = sorted(glob.glob(os.path.join(raw, "*_1.fastq.gz")))
with open(out, "w") as f:
    f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for fp in fwd:
        sid = os.path.basename(fp).replace("_1.fastq.gz", "")
        rp = fp.replace("_1.fastq.gz", "_2.fastq.gz")
        if os.path.exists(rp):
            f.write(f"{sid}\t{fp}\t{rp}\n")
print(f"PE manifest: {len(fwd)} samples -> {out}")
PY
}

write_se_manifest() {
  local raw_dir="$1" out_file="$2"
  python3 - "$raw_dir" "$out_file" <<'PY'
import glob, os, sys
raw, out = sys.argv[1], sys.argv[2]
# Try both .fastq.gz and uncompressed .fastq (Hao/Kawar are uncompressed).
files = sorted(glob.glob(os.path.join(raw, "*.fastq.gz"))
               + glob.glob(os.path.join(raw, "*.fastq")))
# Drop paired suffixes if any sneak in
files = [f for f in files if not f.endswith(("_1.fastq", "_1.fastq.gz",
                                              "_2.fastq", "_2.fastq.gz"))]
with open(out, "w") as f:
    f.write("sample-id\tabsolute-filepath\n")
    for fp in files:
        sid = os.path.basename(fp).replace(".fastq.gz", "").replace(".fastq", "")
        f.write(f"{sid}\t{fp}\n")
print(f"SE manifest: {len(files)} samples -> {out}")
PY
}

reverse_complement_fasta() {
  # Reverse-complement every sequence in $1 -> $2 (preserves header order).
  python3 - "$1" "$2" <<'PY'
import sys
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]
with open(sys.argv[1]) as fin, open(sys.argv[2], "w") as fout:
    hdr, seq = None, []
    def flush():
        if hdr is not None:
            fout.write(hdr + "\n" + rc("".join(seq)) + "\n")
    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            flush(); hdr, seq = line, []
        else:
            seq.append(line)
    flush()
PY
}

min_sample_depth() {
  # Floor of the minimum per-sample total in a QIIME2 FeatureTable[Frequency].
  local table_qza="$1"
  local tmp; tmp=$(mktemp -d)
  qiime tools export --input-path "$table_qza" --output-path "$tmp" >/dev/null
  biom summarize-table -i "$tmp/feature-table.biom" -o "$tmp/summary.txt" >/dev/null
  awk -F': ' '/Min:/{gsub(",", "", $2); printf "%d\n", $2; exit}' "$tmp/summary.txt"
  rm -rf "$tmp"
}

# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

stage_import_and_trim() {
  local study="$1" outdir="$2" raw_dir="$3"
  local layout fwd rev
  layout=$(spec "$study" LAYOUT)
  fwd=$(spec "$study" FWD_PRIMER)
  rev=$(spec "$study" REV_PRIMER)

  cd "$outdir"

  if [[ "$layout" == "PE" ]]; then
    if [[ ! -f manifest.tsv ]]; then
      write_pe_manifest "$raw_dir" manifest.tsv
    fi
    if [[ ! -f paired-end-demux.qza ]]; then
      qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-format PairedEndFastqManifestPhred33V2 \
        --input-path manifest.tsv \
        --output-path paired-end-demux.qza
    fi
    if [[ ! -f trimmed-demux.qza ]]; then
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences paired-end-demux.qza \
        --p-front-f "$fwd" \
        --p-front-r "$rev" \
        --p-error-rate 0.1 \
        --p-cores "$THREADS" \
        --o-trimmed-sequences trimmed-demux.qza \
        --verbose 2>&1 | tail -n 40
    fi
  else
    if [[ ! -f manifest.tsv ]]; then
      write_se_manifest "$raw_dir" manifest.tsv
    fi
    if [[ ! -f single-end-demux.qza ]]; then
      qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-format SingleEndFastqManifestPhred33V2 \
        --input-path manifest.tsv \
        --output-path single-end-demux.qza
    fi
    if [[ ! -f trimmed-demux.qza ]]; then
      if [[ -z "$fwd" ]]; then
        # No 5' primer in the deposited reads (e.g. Kawar 2021 reverse-
        # direction reads): skip cutadapt and pass the demux straight to
        # DADA2 under the trimmed-demux.qza name.
        cp single-end-demux.qza trimmed-demux.qza
      else
        qiime cutadapt trim-single \
          --i-demultiplexed-sequences single-end-demux.qza \
          --p-front "$fwd" \
          --p-error-rate 0.1 \
          --p-cores "$THREADS" \
          --o-trimmed-sequences trimmed-demux.qza \
          --verbose 2>&1 | tail -n 40
      fi
    fi
  fi

  qiime demux summarize \
    --i-data trimmed-demux.qza \
    --o-visualization trimmed-demux.qzv
}

dada2_one_batch() {
  # $1 study, $2 input qza, $3 trunc_f, $4 trunc_r (empty for SE),
  # $5 out_table, $6 out_repseqs, $7 out_stats, $8 layout
  local study="$1" demux="$2" tf="$3" tr="$4"
  local ot="$5" os="$6" oss="$7" layout="$8"
  if [[ "$layout" == "PE" ]]; then
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs "$demux" \
      --p-trunc-len-f "$tf" \
      --p-trunc-len-r "$tr" \
      --p-trim-left-f 0 --p-trim-left-r 0 \
      --p-n-threads "$THREADS" \
      --o-table "$ot" \
      --o-representative-sequences "$os" \
      --o-denoising-stats "$oss" \
      --o-base-transition-stats "base-transition-stats.qza"
  else
    qiime dada2 denoise-single \
      --i-demultiplexed-seqs "$demux" \
      --p-trunc-len "$tf" \
      --p-trim-left 0 \
      --p-n-threads "$THREADS" \
      --o-table "$ot" \
      --o-representative-sequences "$os" \
      --o-denoising-stats "$oss" \
      --o-base-transition-stats "base-transition-stats.qza"
  fi
}

stage_dada2() {
  local study="$1" outdir="$2"
  local layout batch_mode tf tr tf2 tr2
  layout=$(spec "$study" LAYOUT)
  batch_mode=$(spec "$study" BATCH_MODE)
  tf=$(spec "$study" TRUNC_F)
  tr=$(spec "$study" TRUNC_R)

  cd "$outdir"

  if [[ "$batch_mode" == "none" ]]; then
    if [[ ! -f table.qza ]]; then
      dada2_one_batch "$study" trimmed-demux.qza "$tf" "$tr" \
        table.qza rep-seqs.qza denoising-stats.qza "$layout"
    fi
  else
    # Split mode: the user is expected to have produced two demux artifacts,
    # `trimmed-demux-long.qza` and `trimmed-demux-short.qza`, before this
    # stage runs. (Sample lists per batch are determined from the per-read
    # length histogram of `trimmed-demux.qzv` — manual step.)
    tf2=$(spec "$study" BATCH2_TRUNC_F)
    tr2=$(spec "$study" BATCH2_TRUNC_R)
    for required in trimmed-demux-long.qza trimmed-demux-short.qza; do
      if [[ ! -f "$required" ]]; then
        echo "TODO: $study is split-batch; please produce $required by" >&2
        echo "      filtering trimmed-demux.qza on the long/short sample IDs" >&2
        echo "      (see pipeline_log.md 2026-04-14 entries for the split)." >&2
        return 1
      fi
    done
    [[ -f table-long.qza ]] || \
      dada2_one_batch "$study" trimmed-demux-long.qza "$tf" "$tr" \
        table-long.qza rep-seqs-long.qza denoising-stats-long.qza "$layout"
    [[ -f table-short.qza ]] || \
      dada2_one_batch "$study" trimmed-demux-short.qza "$tf2" "$tr2" \
        table-short.qza rep-seqs-short.qza denoising-stats-short.qza "$layout"
    if [[ ! -f table.qza ]]; then
      qiime feature-table merge \
        --i-tables table-long.qza --i-tables table-short.qza \
        --o-merged-table table.qza
      qiime feature-table merge-seqs \
        --i-data rep-seqs-long.qza --i-data rep-seqs-short.qza \
        --o-merged-data rep-seqs.qza
    fi
  fi

  qiime metadata tabulate --m-input-file denoising-stats.qza \
    --o-visualization denoising-stats.qzv 2>/dev/null || true
  qiime feature-table summarize --i-table table.qza --o-visualization table.qzv
  qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
}

stage_taxonomy() {
  local study="$1" outdir="$2"
  local region silva ehomd
  region=$(spec "$study" REGION)
  silva=$(silva_classifier "$region")
  ehomd=$(ehomd_classifier "$region")

  cd "$outdir"

  for clf_pair in "silva:$silva" "ehomd:$ehomd"; do
    local tag="${clf_pair%%:*}"
    local clf="${clf_pair##*:}"
    if [[ ! -f "$clf" ]]; then
      echo "MISSING classifier: $clf" >&2; return 1
    fi
    if [[ ! -f "taxonomy-${tag}.qza" ]]; then
      qiime feature-classifier classify-sklearn \
        --i-classifier "$clf" \
        --i-reads rep-seqs.qza \
        --p-n-jobs "$THREADS" \
        --o-classification "taxonomy-${tag}.qza"
    fi
    qiime metadata tabulate --m-input-file "taxonomy-${tag}.qza" \
      --o-visualization "taxonomy-${tag}.qzv"
    if [[ -f metadata.tsv ]]; then
      qiime taxa barplot \
        --i-table table.qza --i-taxonomy "taxonomy-${tag}.qza" \
        --m-metadata-file metadata.tsv \
        --o-visualization "taxa-bar-plots-${tag}.qzv"
    fi
  done
}

stage_phylogeny() {
  local outdir="$1"
  cd "$outdir"
  if [[ ! -f rooted-tree.qza ]]; then
    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --p-n-threads "$THREADS" \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza
  fi
  mkdir -p exported-tree
  qiime tools export --input-path rooted-tree.qza --output-path exported-tree/
  cp -f exported-tree/tree.nwk rooted-tree.nwk
}

stage_genus_and_export() {
  local outdir="$1"
  cd "$outdir"
  for tag in silva ehomd; do
    if [[ ! -f "genus-table-${tag}.qza" ]]; then
      qiime taxa collapse \
        --i-table table.qza --i-taxonomy "taxonomy-${tag}.qza" \
        --p-level 6 --o-collapsed-table "genus-table-${tag}.qza"
    fi
    mkdir -p "exported-${tag}"
    qiime tools export --input-path "genus-table-${tag}.qza" \
      --output-path "exported-${tag}/"
    biom convert -i "exported-${tag}/feature-table.biom" \
      -o "exported-${tag}/genus-table.tsv" --to-tsv
    qiime tools export --input-path "taxonomy-${tag}.qza" --output-path "exported-${tag}/"
    qiime tools export --input-path table.qza --output-path "exported-${tag}/"
    biom convert -i "exported-${tag}/feature-table.biom" \
      -o "exported-${tag}/asv-table.tsv" --to-tsv
  done
}

# NOTE: stage_core_metrics() removed deliberately. Diversity (alpha + beta,
# incl. weighted UniFrac) is computed downstream by the R analysis from
# exported-{silva,ehomd}/asv-table.tsv + rooted-tree.nwk + metadata.tsv;
# running core-metrics-phylogenetic here would be wasted work and could
# mislead a reader into thinking those numbers come from QIIME2.
#
# If you ever want a QIIME2 diversity snapshot for sanity-checking, run:
#   qiime diversity core-metrics-phylogenetic \
#     --i-phylogeny rooted-tree.qza --i-table table.qza \
#     --p-sampling-depth "$(min_sample_depth table.qza)" \
#     --m-metadata-file metadata.tsv \
#     --output-dir core-metrics-results/
# but do not include those numbers in the manuscript.

stage_picrust2() {
  [[ "${SKIP_PICRUST2:-0}" == "1" ]] && return 0
  local study="$1" outdir="$2"
  local rc
  rc=$(spec "$study" RC_REP_SEQS)

  cd "$outdir"

  # Export 2025.10 artefacts -> re-import under 2024.5 env.
  mkdir -p _picrust_reimport
  qiime tools export --input-path table.qza    --output-path _picrust_reimport/ >/dev/null
  qiime tools export --input-path rep-seqs.qza --output-path _picrust_reimport/ >/dev/null

  if [[ "$rc" == "1" ]]; then
    reverse_complement_fasta _picrust_reimport/dna-sequences.fasta \
      _picrust_reimport/dna-sequences.rc.fasta
    mv _picrust_reimport/dna-sequences.rc.fasta _picrust_reimport/dna-sequences.fasta
  fi

  # Switch envs for PICRUSt2.
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$PICRUST_ENV"

  qiime tools import \
    --type 'FeatureTable[Frequency]' \
    --input-path _picrust_reimport/feature-table.biom \
    --output-path _picrust_reimport/table_24.qza

  qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path _picrust_reimport/dna-sequences.fasta \
    --output-path _picrust_reimport/rep-seqs_24.qza

  rm -rf q2-picrust2_output
  qiime picrust2 full-pipeline \
    --i-table _picrust_reimport/table_24.qza \
    --i-seq  _picrust_reimport/rep-seqs_24.qza \
    --output-dir q2-picrust2_output \
    --p-placement-tool epa-ng \
    --p-threads "$THREADS" \
    --p-hsp-method mp \
    --p-max-nsti 2 \
    --verbose 2>&1 | tee picrust2_run.log

  mkdir -p exported-picrust2
  for t in pathway_abundance ko_metagenome ec_metagenome; do
    [[ -f q2-picrust2_output/${t}.qza ]] || continue
    qiime tools export --input-path "q2-picrust2_output/${t}.qza" \
      --output-path "exported-picrust2/${t}/"
    biom convert -i "exported-picrust2/${t}/feature-table.biom" \
      -o "exported-picrust2/${t}.tsv" --to-tsv
  done

  conda activate "$QIIME_ENV"
}

stage_maaslin2_inputs() {
  # MaAsLin2 / phyloseq inputs are already produced by stage_genus_and_export
  # (exported-silva/, exported-ehomd/) and stage_picrust2 (exported-picrust2/).
  # We just ensure the per-study sample-metadata.tsv is alongside them.
  local outdir="$1"
  cd "$outdir"
  if [[ -f metadata.tsv ]]; then
    cp -f metadata.tsv exported-silva/sample-metadata.tsv  2>/dev/null || true
    cp -f metadata.tsv exported-ehomd/sample-metadata.tsv  2>/dev/null || true
    cp -f metadata.tsv exported-picrust2/sample-metadata.tsv 2>/dev/null || true
  else
    echo "TODO: $outdir is missing metadata.tsv (group / sample-id mapping)" >&2
  fi
}

# ---------------------------------------------------------------------------
# Per-study driver
# ---------------------------------------------------------------------------
process_study() {
  local study="$1"
  local raw_subdir out_subdir raw_dir outdir
  raw_subdir=$(spec "$study" RAW_DIR)
  out_subdir=$(spec "$study" OUT_SUBDIR)
  if [[ -z "$raw_subdir" || -z "$out_subdir" ]]; then
    echo "Unknown study: $study" >&2; return 1
  fi
  raw_dir="$RAW_ROOT/$raw_subdir"
  outdir="$OUT_ROOT/$out_subdir"
  mkdir -p "$outdir"

  echo
  echo "======================================================================"
  echo "  $study  ($(spec "$study" DISEASE), $(spec "$study" REGION), $(spec "$study" LAYOUT))"
  echo "======================================================================"

  stage_import_and_trim   "$study" "$outdir" "$raw_dir"
  stage_dada2             "$study" "$outdir"
  stage_taxonomy          "$study" "$outdir"
  stage_phylogeny                  "$outdir"
  stage_genus_and_export           "$outdir"
  stage_picrust2          "$study" "$outdir"
  stage_maaslin2_inputs            "$outdir"

  echo "  $study done."
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
main() {
  if [[ -z "${CONDA_PREFIX:-}" || "$(basename "$CONDA_PREFIX")" != "$QIIME_ENV" ]]; then
    echo "ERROR: activate the QIIME2 env first:  conda activate $QIIME_ENV" >&2
    exit 1
  fi
  local studies=( "$@" )
  [[ ${#studies[@]} -eq 0 ]] && studies=( "${DEFAULT_STUDIES[@]}" )
  for s in "${studies[@]}"; do
    process_study "$s"
  done
  echo
  echo "All requested studies finished."
}

main "$@"
