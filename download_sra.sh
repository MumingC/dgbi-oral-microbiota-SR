#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# download.sh — fetch raw 16S FASTQs for the 7 included studies
#
# Reproduces the SRA-pull loop actually used to populate ~/sra_downloads/.
# Source of the loop body: bash_history line 654-655 (single inline run).
#
# Each study is downloaded into its own subdirectory of $OUT_ROOT. For every
# run accession we:
#   1. prefetch the SRA blob (max-size 50G to allow large lanes)
#   2. fasterq-dump --split-3 (writes _1/_2 for PE, .fastq for SE)
#   3. gzip the FASTQs
#   4. remove the SRA blob to save disk
#
# Requires SRA Toolkit (>= 3.x; we used 3.4.1 — installed from
# https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz)
# and EDirect (esearch, efetch) for run-list resolution.
# ---------------------------------------------------------------------------

set -euo pipefail

OUT_ROOT="${OUT_ROOT:-$HOME/sra_downloads}"
THREADS="${THREADS:-4}"
MAX_SRA_SIZE="${MAX_SRA_SIZE:-50G}"

# BioProject : local-directory-name. Order matches the original pull.
STUDIES=(
  "PRJNA598080:Ziganshina_2020_GERD"   # PE V3-V4, 26 samples
  "PRJNA824804:Hao_2022_GERD"          # SE V3-V4, 48 samples (uncompressed .fastq)
  "PRJNA674379:Kawar_2021_GERD"        # SE V1-V3, 128 samples (uncompressed .fastq)
  "PRJNA894717:Qian_2023_GERD"         # PE V3-V4, 60 samples
  "PRJNA996485:Zheng_2024_LPRD"        # PE V3-V4, 252 samples (two seq batches)
  "PRJNA873889:Tang_2023_IBS"          # PE V4-V5, 124 samples (two seq batches)
  "PRJNA1118066:Li_2025_IBS"           # PE V3-V4, 35 samples
)

# --- preflight ---------------------------------------------------------------
for tool in prefetch fasterq-dump esearch efetch gzip; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required tool not on PATH: $tool" >&2
    exit 1
  fi
done

mkdir -p "$OUT_ROOT"

download_run() {
  local srr="$1"
  echo "  $(date +%F_%T)  $srr"
  prefetch "$srr" --max-size "$MAX_SRA_SIZE"
  fasterq-dump "$srr" --split-3 --threads "$THREADS"
  rm -rf "$srr"
  # gzip whatever fasterq-dump produced (paired -> _1/_2, single -> .fastq)
  if [[ -f "${srr}_1.fastq" ]]; then
    gzip "${srr}_1.fastq" "${srr}_2.fastq"
  elif [[ -f "${srr}.fastq" ]]; then
    # NOTE: Hao 2022 and Kawar 2021 were kept uncompressed in the original run
    # (QIIME2 accepts both). Comment out the next line if you want to match
    # the original layout exactly.
    gzip "${srr}.fastq"
  else
    echo "  WARNING: no FASTQ output found for $srr" >&2
  fi
}

# --- main loop ---------------------------------------------------------------
for entry in "${STUDIES[@]}"; do
  prj="${entry%%:*}"
  name="${entry##*:}"
  dir="$OUT_ROOT/$name"
  mkdir -p "$dir"
  cd "$dir"

  echo "=== $(date +%F_%T) START $name ($prj) ==="

  # Resolve run accessions from SRA for this BioProject.
  esearch -db sra -query "${prj}[BioProject]" \
    | efetch -format runinfo \
    | cut -d',' -f1 \
    | grep -v "^Run" \
    | grep -v "^$" \
    | while read -r srr; do
        download_run "$srr"
      done

  echo "=== $(date +%F_%T) END   $name ==="
done

echo "All downloads complete. Raw FASTQ under: $OUT_ROOT"
