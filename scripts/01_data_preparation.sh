#!/bin/bash
#
# REGENIE GWAS - Step 1: Data Preparation and Quality Control
# This script performs quality control on genotype data before REGENIE analysis
#
# Usage: bash 01_data_preparation.sh
#

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# Configuration
# =============================================================================

# Input files
GENOTYPE_PREFIX="data/genotypes"  # PLINK format (.bed/.bim/.fam)
PHENOTYPE_FILE="data/phenotypes.txt"
COVARIATE_FILE="data/covariates.txt"

# Output directory
OUTPUT_DIR="output/qc"
mkdir -p ${OUTPUT_DIR}

# QC parameters
MAF_THRESHOLD=0.01          # Minor allele frequency threshold
GENO_THRESHOLD=0.1          # Genotype missingness threshold
MIND_THRESHOLD=0.1          # Individual missingness threshold
HWE_THRESHOLD=1e-15         # Hardy-Weinberg equilibrium p-value threshold

# =============================================================================
# Functions
# =============================================================================

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

check_file() {
    if [ ! -f "$1" ]; then
        log "ERROR: File not found: $1"
        exit 1
    fi
}

# =============================================================================
# Quality Control Pipeline
# =============================================================================

log "Starting data preparation and QC..."

# Check input files exist
log "Checking input files..."
check_file "${GENOTYPE_PREFIX}.bed"
check_file "${GENOTYPE_PREFIX}.bim"
check_file "${GENOTYPE_PREFIX}.fam"
check_file "${PHENOTYPE_FILE}"
check_file "${COVARIATE_FILE}"

# Step 1: Calculate missingness
log "Step 1: Calculating missingness rates..."
plink \
    --bfile ${GENOTYPE_PREFIX} \
    --missing \
    --out ${OUTPUT_DIR}/missingness

# Step 2: Remove variants and individuals with high missingness
log "Step 2: Filtering by missingness (variants < ${GENO_THRESHOLD}, individuals < ${MIND_THRESHOLD})..."
plink \
    --bfile ${GENOTYPE_PREFIX} \
    --geno ${GENO_THRESHOLD} \
    --mind ${MIND_THRESHOLD} \
    --make-bed \
    --out ${OUTPUT_DIR}/step1_missingness

# Step 3: Filter by MAF
log "Step 3: Filtering by MAF > ${MAF_THRESHOLD}..."
plink \
    --bfile ${OUTPUT_DIR}/step1_missingness \
    --maf ${MAF_THRESHOLD} \
    --make-bed \
    --out ${OUTPUT_DIR}/step2_maf

# Step 4: Filter by HWE
log "Step 4: Filtering by Hardy-Weinberg equilibrium (p > ${HWE_THRESHOLD})..."
plink \
    --bfile ${OUTPUT_DIR}/step2_maf \
    --hwe ${HWE_THRESHOLD} \
    --make-bed \
    --out ${OUTPUT_DIR}/step3_hwe

# Step 5: Calculate allele frequencies
log "Step 5: Calculating allele frequencies..."
plink \
    --bfile ${OUTPUT_DIR}/step3_hwe \
    --freq \
    --out ${OUTPUT_DIR}/allele_frequencies

# Step 6: Generate QC report
log "Step 6: Generating QC summary report..."
cat > ${OUTPUT_DIR}/qc_report.txt <<EOF
REGENIE GWAS - Quality Control Report
Generated: $(date)
========================================

Input Data:
-----------
Genotype file: ${GENOTYPE_PREFIX}
Phenotype file: ${PHENOTYPE_FILE}
Covariate file: ${COVARIATE_FILE}

QC Parameters:
--------------
MAF threshold: ${MAF_THRESHOLD}
Genotype missingness: ${GENO_THRESHOLD}
Individual missingness: ${MIND_THRESHOLD}
HWE p-value threshold: ${HWE_THRESHOLD}

Results Summary:
----------------
Initial variants: $(wc -l < ${GENOTYPE_PREFIX}.bim)
Initial individuals: $(wc -l < ${GENOTYPE_PREFIX}.fam)

After QC:
Variants passed QC: $(wc -l < ${OUTPUT_DIR}/step3_hwe.bim)
Individuals passed QC: $(wc -l < ${OUTPUT_DIR}/step3_hwe.fam)

Output Files:
-------------
QC-filtered genotypes: ${OUTPUT_DIR}/step3_hwe.*
Missingness report: ${OUTPUT_DIR}/missingness.*
Allele frequencies: ${OUTPUT_DIR}/allele_frequencies.frq

Next Step:
----------
Run: bash scripts/02_regenie_step1.sh
EOF

cat ${OUTPUT_DIR}/qc_report.txt

log "Data preparation and QC complete!"
log "QC-filtered data: ${OUTPUT_DIR}/step3_hwe.*"
log "QC report: ${OUTPUT_DIR}/qc_report.txt"
