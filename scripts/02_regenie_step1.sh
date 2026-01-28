#!/bin/bash
#
# REGENIE GWAS - Step 2: REGENIE Step 1 (Whole Genome Regression)
# This script runs REGENIE step 1 to fit whole genome regression model
#
# Usage: bash 02_regenie_step1.sh
#

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# Configuration
# =============================================================================

# Input files from QC step
QC_GENOTYPE="output/qc/step3_hwe"
PHENOTYPE_FILE="data/phenotypes.txt"
COVARIATE_FILE="data/covariates.txt"

# Output directory
OUTPUT_DIR="output/regenie_step1"
mkdir -p ${OUTPUT_DIR}

# REGENIE Step 1 parameters
THREADS=8                    # Number of threads to use
BSIZE=1000                   # Block size (number of SNPs per block)
MIN_MAC=20                   # Minimum minor allele count

# Create SNP list for Step 1 (use pruned SNPs for better efficiency)
STEP1_SNPS="${OUTPUT_DIR}/step1_snps"

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
# REGENIE Step 1 Pipeline
# =============================================================================

log "Starting REGENIE Step 1 (Whole Genome Regression)..."

# Check input files
log "Checking input files..."
check_file "${QC_GENOTYPE}.bed"
check_file "${QC_GENOTYPE}.bim"
check_file "${QC_GENOTYPE}.fam"
check_file "${PHENOTYPE_FILE}"
check_file "${COVARIATE_FILE}"

# Step 1: LD pruning to select independent SNPs for Step 1
log "Step 1: Performing LD pruning for variant selection..."
plink \
    --bfile ${QC_GENOTYPE} \
    --indep-pairwise 1000 100 0.9 \
    --out ${STEP1_SNPS}

log "Selected $(wc -l < ${STEP1_SNPS}.prune.in) independent SNPs for REGENIE Step 1"

# Extract pruned SNPs
plink \
    --bfile ${QC_GENOTYPE} \
    --extract ${STEP1_SNPS}.prune.in \
    --make-bed \
    --out ${STEP1_SNPS}

# Step 2: Run REGENIE Step 1
log "Step 2: Running REGENIE Step 1..."
log "This may take some time depending on dataset size..."

regenie \
    --step 1 \
    --bed ${STEP1_SNPS} \
    --phenoFile ${PHENOTYPE_FILE} \
    --covarFile ${COVARIATE_FILE} \
    --bsize ${BSIZE} \
    --minMAC ${MIN_MAC} \
    --lowmem \
    --lowmem-prefix ${OUTPUT_DIR}/tmp_rg \
    --threads ${THREADS} \
    --out ${OUTPUT_DIR}/fit_out

# Step 3: Verify output files
log "Step 3: Verifying REGENIE Step 1 output..."

if [ -f "${OUTPUT_DIR}/fit_out_pred.list" ]; then
    log "SUCCESS: REGENIE Step 1 predictions file created"
else
    log "ERROR: REGENIE Step 1 failed - predictions file not found"
    exit 1
fi

# Generate summary report
log "Step 4: Generating summary report..."
cat > ${OUTPUT_DIR}/step1_report.txt <<EOF
REGENIE Step 1 - Whole Genome Regression Report
Generated: $(date)
================================================

Input Data:
-----------
Genotype file: ${QC_GENOTYPE}
Phenotype file: ${PHENOTYPE_FILE}
Covariate file: ${COVARIATE_FILE}

Step 1 Parameters:
------------------
Block size: ${BSIZE}
Minimum MAC: ${MIN_MAC}
Threads: ${THREADS}
Pruned SNPs used: $(wc -l < ${STEP1_SNPS}.bim)

Output Files:
-------------
Predictions: ${OUTPUT_DIR}/fit_out_pred.list
LOCO predictions: ${OUTPUT_DIR}/fit_out_*.loco
Log file: ${OUTPUT_DIR}/fit_out.log

Next Step:
----------
Run: bash scripts/03_regenie_step2.sh
EOF

cat ${OUTPUT_DIR}/step1_report.txt

log "REGENIE Step 1 complete!"
log "Predictions file: ${OUTPUT_DIR}/fit_out_pred.list"
log "Report: ${OUTPUT_DIR}/step1_report.txt"
