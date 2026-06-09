#!/bin/bash
# Download raw 16S FASTQ for the oral microbiota / DGBI systematic review.
# Pipeline per run: prefetch (.sra) -> fasterq-dump (FASTQ) -> gzip.
#
# Requires SRA Toolkit >= 3.0 (prefetch, fasterq-dump) and Entrez Direct
# (esearch/efetch) to resolve SRR accessions from BioProject / DRA IDs.
# pigz is used for compression if available, otherwise gzip.
#
# Usage: set OUTPUT_DIR below, then  chmod +x download_sra.sh && ./download_sra.sh

set -euo pipefail

# ======================== USER SETTINGS ========================
OUTPUT_DIR="./sra_downloads"       # <-- Change this to your desired path
THREADS=6                          # Number of threads for fasterq-dump
TEMP_DIR="./sra_tmp"               # Temporary directory for fasterq-dump
# ===============================================================

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/download_log_$(date +%Y%m%d_%H%M%S).txt"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "SRA Download Script - Oral Microbiota SR"
echo "Started: $(date)"
echo "Output:  $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "========================================"

###############################################################################
# STEP 1: Define BioProject accessions
###############################################################################

# NCBI SRA BioProjects
declare -A NCBI_PROJECTS
NCBI_PROJECTS=(
    ["Hao_2022_GERD"]="PRJNA824804"
    # NOTE: PRJNA46307 (Hao 2022) is the HMP project with thousands of runs.
    # The paper likely reused public HMP data as controls. You may want to
    # skip this or download only specific samples. Uncomment if needed:
    # ["Hao_2022_HMP_controls"]="PRJNA46307"
    ["Kawar_2021_GERD"]="PRJNA674379"
    ["Ziganshina_2020_GERD"]="PRJNA598080"
    ["Qian_2023_GERD"]="PRJNA894717"
    ["Zheng_2024_LPRD"]="PRJNA996485"
    ["Tang_2023_IBS"]="PRJNA873889"
    ["Li_2025_IBS"]="PRJNA1118066"
    ["Kwiatkowska_2024_FC"]="PRJNA925675"
)

# DDBJ DRA accessions (synced to NCBI, so SRA Toolkit can access them)
declare -A DDBJ_PROJECTS
DDBJ_PROJECTS=(
    ["Tanaka_2022_IBS"]="DRA013075"
    ["Kim_2023_FD"]="DRA015691"
)

###############################################################################
# STEP 2: Function to get SRR accessions from BioProject
###############################################################################

get_srr_from_bioproject() {
    local project_id="$1"
    local study_name="$2"
    local srr_file="${OUTPUT_DIR}/${study_name}_SRR_list.txt"

    echo ""
    echo "--- Fetching SRR accessions for ${study_name} (${project_id}) ---"

    # Use EDirect to get run accessions
    esearch -db sra -query "${project_id}[BioProject]" | \
        efetch -format runinfo | \
        grep -v "^Run" | \
        cut -d',' -f1 | \
        grep -v "^$" > "$srr_file"

    local count=$(wc -l < "$srr_file")
    echo "  Found ${count} runs, saved to ${srr_file}"
    echo "$srr_file"
}

###############################################################################
# STEP 3: Function to get SRR accessions from DDBJ DRA
###############################################################################

get_srr_from_dra() {
    local dra_id="$1"
    local study_name="$2"
    local srr_file="${OUTPUT_DIR}/${study_name}_SRR_list.txt"

    echo ""
    echo "--- Fetching SRR accessions for ${study_name} (${dra_id}) ---"

    # DDBJ DRA data is mirrored to NCBI, so we can search by DRA accession
    esearch -db sra -query "${dra_id}" | \
        efetch -format runinfo | \
        grep -v "^Run" | \
        cut -d',' -f1 | \
        grep -v "^$" > "$srr_file"

    # If no results via DRA ID, try alternative search
    if [ ! -s "$srr_file" ]; then
        echo "  Direct DRA search returned 0 results. Trying alternative..."
        esearch -db sra -query "${dra_id}[Accession]" | \
            efetch -format runinfo | \
            grep -v "^Run" | \
            cut -d',' -f1 | \
            grep -v "^$" > "$srr_file"
    fi

    local count=$(wc -l < "$srr_file")
    echo "  Found ${count} runs, saved to ${srr_file}"
    echo "$srr_file"
}

###############################################################################
# STEP 4: Function to download and convert one SRR accession
###############################################################################

download_and_convert() {
    local srr_id="$1"
    local study_dir="$2"

    echo "  [prefetch]     ${srr_id} ..."
    prefetch "$srr_id" --output-directory "$study_dir" --max-size 50G

    echo "  [fasterq-dump] ${srr_id} to FASTQ ..."
    fasterq-dump "$srr_id" \
        --split-3 \
        --outdir "$study_dir" \
        --temp "$TEMP_DIR" \
        --threads "$THREADS"

    echo "  [gzip]         Compressing FASTQ files ..."
    # Use pigz if available (parallel gzip), otherwise gzip
    if command -v pigz &> /dev/null; then
        pigz -p "$THREADS" "${study_dir}/${srr_id}"*.fastq 2>/dev/null || true
    else
        gzip "${study_dir}/${srr_id}"*.fastq 2>/dev/null || true
    fi

    # Clean up .sra file to save space
    rm -rf "${study_dir}/${srr_id}" 2>/dev/null || true

    echo "  [done]         ${srr_id}"
}

###############################################################################
# STEP 5: Main download loop — NCBI projects
###############################################################################

echo ""
echo "========================================"
echo "Processing NCBI SRA BioProjects..."
echo "========================================"

for study_name in "${!NCBI_PROJECTS[@]}"; do
    project_id="${NCBI_PROJECTS[$study_name]}"
    study_dir="${OUTPUT_DIR}/${study_name}"
    mkdir -p "$study_dir"

    srr_file=$(get_srr_from_bioproject "$project_id" "$study_name")

    if [ ! -s "$srr_file" ]; then
        echo "  WARNING: No SRR accessions found for ${study_name}. Skipping."
        continue
    fi

    echo "  Downloading runs for ${study_name}..."
    while IFS= read -r srr_id; do
        [ -z "$srr_id" ] && continue
        # Skip if already downloaded
        if ls "${study_dir}/${srr_id}"*.fastq.gz 1>/dev/null 2>&1; then
            echo "  [skip] ${srr_id} already exists"
            continue
        fi
        download_and_convert "$srr_id" "$study_dir"
    done < "$srr_file"

    echo "=== ${study_name} complete ==="
done

###############################################################################
# STEP 6: Main download loop — DDBJ projects
###############################################################################

echo ""
echo "========================================"
echo "Processing DDBJ DRA projects..."
echo "========================================"

for study_name in "${!DDBJ_PROJECTS[@]}"; do
    dra_id="${DDBJ_PROJECTS[$study_name]}"
    study_dir="${OUTPUT_DIR}/${study_name}"
    mkdir -p "$study_dir"

    srr_file=$(get_srr_from_dra "$dra_id" "$study_name")

    if [ ! -s "$srr_file" ]; then
        echo "  WARNING: No runs found for ${study_name} (${dra_id})."
        echo "  You may need to download manually from DDBJ:"
        echo "    https://ddbj.nig.ac.jp/resource/sra-submission/${dra_id}"
        continue
    fi

    echo "  Downloading runs for ${study_name}..."
    while IFS= read -r srr_id; do
        [ -z "$srr_id" ] && continue
        if ls "${study_dir}/${srr_id}"*.fastq.gz 1>/dev/null 2>&1; then
            echo "  [skip] ${srr_id} already exists"
            continue
        fi
        download_and_convert "$srr_id" "$study_dir"
    done < "$srr_file"

    echo "=== ${study_name} complete ==="
done

###############################################################################
# STEP 7: Summary
###############################################################################

echo ""
echo "========================================"
echo "Download Summary"
echo "========================================"
echo "Completed: $(date)"
echo ""
echo "Files per study:"
for dir in "${OUTPUT_DIR}"/*/; do
    if [ -d "$dir" ]; then
        study=$(basename "$dir")
        fq_count=$(ls "$dir"/*.fastq.gz 2>/dev/null | wc -l || echo 0)
        echo "  ${study}: ${fq_count} FASTQ files"
    fi
done
echo ""
echo "Log saved to: ${LOG_FILE}"
echo "========================================"
