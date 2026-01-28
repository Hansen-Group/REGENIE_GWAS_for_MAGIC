#!/bin/bash
#
# REGENIE GWAS - Step 3: REGENIE Step 2 (Association Testing)
# This script runs REGENIE step 2 to perform association testing
#
# Usage: bash 03_regenie_step2.sh
#

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# Configuration
# =============================================================================

# Input files
QC_GENOTYPE="output/qc/step3_hwe"
PHENOTYPE_FILE="data/phenotypes.txt"
COVARIATE_FILE="data/covariates.txt"
STEP1_PRED="output/regenie_step1/fit_out_pred.list"

# Output directory
OUTPUT_DIR="output/regenie_step2"
mkdir -p ${OUTPUT_DIR}

# REGENIE Step 2 parameters
THREADS=8                    # Number of threads to use
BSIZE=400                    # Block size for step 2
MIN_MAC=5                    # Minimum minor allele count
MIN_INFO=0.6                 # Minimum imputation INFO score (if using imputed data)

# Association test parameters
P_THRESHOLD=5e-8             # Genome-wide significance threshold

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
# REGENIE Step 2 Pipeline
# =============================================================================

log "Starting REGENIE Step 2 (Association Testing)..."

# Check input files
log "Checking input files..."
check_file "${QC_GENOTYPE}.bed"
check_file "${QC_GENOTYPE}.bim"
check_file "${QC_GENOTYPE}.fam"
check_file "${PHENOTYPE_FILE}"
check_file "${COVARIATE_FILE}"
check_file "${STEP1_PRED}"

# Run REGENIE Step 2
log "Step 1: Running REGENIE Step 2 association testing..."
log "This may take considerable time for large datasets..."

regenie \
    --step 2 \
    --bed ${QC_GENOTYPE} \
    --phenoFile ${PHENOTYPE_FILE} \
    --covarFile ${COVARIATE_FILE} \
    --pred ${STEP1_PRED} \
    --bsize ${BSIZE} \
    --minMAC ${MIN_MAC} \
    --threads ${THREADS} \
    --out ${OUTPUT_DIR}/assoc_results

# Step 2: Count significant associations
log "Step 2: Analyzing results..."

# Count total variants tested
TOTAL_VARIANTS=$(grep -v "^#" ${OUTPUT_DIR}/assoc_results_*.regenie | grep -v "CHROM" | wc -l || echo "0")

# Count genome-wide significant variants
# Note: REGENIE outputs LOG10P in column 12 (standard format)
# For p-value threshold, we need -log10(p) comparison
LOG10P_THRESHOLD=$(awk -v p=${P_THRESHOLD} 'BEGIN {print -log(p)/log(10)}')
SIGNIFICANT_VARIANTS=$(awk -v thresh=${LOG10P_THRESHOLD} 'NR>1 && $12>thresh' ${OUTPUT_DIR}/assoc_results_*.regenie | wc -l || echo "0")

log "Total variants tested: ${TOTAL_VARIANTS}"
log "Genome-wide significant variants (p < ${P_THRESHOLD}): ${SIGNIFICANT_VARIANTS}"

# Step 3: Extract top hits
log "Step 3: Extracting top associations..."

# Extract variants with p < 1e-5 for review (LOG10P > 5)
awk 'NR==1 || $12>5' ${OUTPUT_DIR}/assoc_results_*.regenie \
    > ${OUTPUT_DIR}/top_associations_p1e-5.txt

# Extract genome-wide significant variants (LOG10P > 7.3 for p < 5e-8)
LOG10P_THRESHOLD=$(awk -v p=${P_THRESHOLD} 'BEGIN {print -log(p)/log(10)}')
awk -v thresh="${LOG10P_THRESHOLD}" 'NR==1 || $12>thresh' ${OUTPUT_DIR}/assoc_results_*.regenie \
    > ${OUTPUT_DIR}/significant_associations.txt

# Step 4: Generate summary report
log "Step 4: Generating summary report..."

cat > ${OUTPUT_DIR}/step2_report.txt <<EOF
REGENIE Step 2 - Association Testing Report
Generated: $(date)
============================================

Input Data:
-----------
Genotype file: ${QC_GENOTYPE}
Phenotype file: ${PHENOTYPE_FILE}
Covariate file: ${COVARIATE_FILE}
Step 1 predictions: ${STEP1_PRED}

Step 2 Parameters:
------------------
Block size: ${BSIZE}
Minimum MAC: ${MIN_MAC}
Threads: ${THREADS}

Results Summary:
----------------
Total variants tested: ${TOTAL_VARIANTS}
Suggestive associations (p < 1e-5): $(wc -l < ${OUTPUT_DIR}/top_associations_p1e-5.txt | awk '{print $1-1}')
Genome-wide significant (p < ${P_THRESHOLD}): ${SIGNIFICANT_VARIANTS}

Output Files:
-------------
Full results: ${OUTPUT_DIR}/assoc_results_*.regenie
Top associations (p < 1e-5): ${OUTPUT_DIR}/top_associations_p1e-5.txt
Significant hits (p < ${P_THRESHOLD}): ${OUTPUT_DIR}/significant_associations.txt
Log file: ${OUTPUT_DIR}/assoc_results.log

Top 10 Associations:
--------------------
EOF

# Add top 10 associations to report
awk 'NR==1 || NR<=11' ${OUTPUT_DIR}/top_associations_p1e-5.txt >> ${OUTPUT_DIR}/step2_report.txt

cat >> ${OUTPUT_DIR}/step2_report.txt <<EOF

Next Step:
----------
Run: bash scripts/04_process_results.sh
EOF

cat ${OUTPUT_DIR}/step2_report.txt

log "REGENIE Step 2 complete!"
log "Results: ${OUTPUT_DIR}/assoc_results_*.regenie"
log "Report: ${OUTPUT_DIR}/step2_report.txt"
