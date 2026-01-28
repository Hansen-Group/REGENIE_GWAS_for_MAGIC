#!/bin/bash
#
# REGENIE GWAS - Step 4: Process and Visualize Results
# This script processes REGENIE results and generates summary reports
#
# Usage: bash 04_process_results.sh
#

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# Configuration
# =============================================================================

# Input files
RESULTS_DIR="output/regenie_step2"
RESULTS_FILE="${RESULTS_DIR}/assoc_results_*.regenie"

# Output directory
OUTPUT_DIR="output/processed_results"
mkdir -p ${OUTPUT_DIR}

# Processing parameters
P_THRESHOLD=5e-8             # Genome-wide significance threshold
SUGGESTIVE_P=1e-5            # Suggestive threshold

# =============================================================================
# Functions
# =============================================================================

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# =============================================================================
# Results Processing Pipeline
# =============================================================================

log "Starting results processing and visualization..."

# Step 1: Combine results if multiple files
log "Step 1: Combining results files..."

# Find all result files
RESULT_FILES=$(ls ${RESULTS_FILE} 2>/dev/null || echo "")

if [ -z "${RESULT_FILES}" ]; then
    log "ERROR: No results files found matching ${RESULTS_FILE}"
    exit 1
fi

# Combine all results into one file
first_file=true
for file in ${RESULT_FILES}; do
    if [ "${first_file}" = true ]; then
        cat ${file} > ${OUTPUT_DIR}/combined_results.txt
        first_file=false
    else
        tail -n +2 ${file} >> ${OUTPUT_DIR}/combined_results.txt
    fi
done

# Step 2: Generate summary statistics
log "Step 2: Generating summary statistics..."

python3 - <<'PYTHON_SCRIPT'
import sys
import pandas as pd
import numpy as np

# Read results
try:
    df = pd.read_csv('output/processed_results/combined_results.txt', sep='\s+')
except Exception as e:
    print(f"Error reading results: {e}")
    sys.exit(1)

# Calculate genomic inflation factor (lambda)
if 'LOG10P' in df.columns:
    # Convert LOG10P to p-values
    df['P'] = 10 ** (-df['LOG10P'])
elif 'P' in df.columns:
    pass
else:
    print("Error: No p-value column found")
    sys.exit(1)

# Calculate lambda
chisq = -2 * np.log(df['P'])
median_chisq = np.median(chisq)
lambda_gc = median_chisq / 0.4549364  # 0.4549364 is chi-square(1) median

# Summary statistics
summary = {
    'total_variants': len(df),
    'lambda_gc': lambda_gc,
    'min_p': df['P'].min(),
    'significant_5e8': (df['P'] < 5e-8).sum(),
    'suggestive_1e5': (df['P'] < 1e-5).sum()
}

# Save summary
with open('output/processed_results/summary_stats.txt', 'w') as f:
    f.write("GWAS Summary Statistics\n")
    f.write("=" * 50 + "\n")
    f.write(f"Total variants tested: {summary['total_variants']:,}\n")
    f.write(f"Genomic inflation factor (λ): {summary['lambda_gc']:.3f}\n")
    f.write(f"Minimum p-value: {summary['min_p']:.2e}\n")
    f.write(f"Genome-wide significant (p < 5e-8): {summary['significant_5e8']}\n")
    f.write(f"Suggestive associations (p < 1e-5): {summary['suggestive_1e5']}\n")

print(f"Summary statistics calculated: λ = {summary['lambda_gc']:.3f}")
PYTHON_SCRIPT

cat ${OUTPUT_DIR}/summary_stats.txt

# Step 3: Extract significant loci
log "Step 3: Extracting significant loci..."

# Calculate LOG10P threshold from p-value threshold
LOG10P_THRESHOLD=$(awk -v p=${P_THRESHOLD} 'BEGIN {print -log(p)/log(10)}')

# Extract genome-wide significant variants (REGENIE outputs LOG10P in column 12)
awk -v thresh="${LOG10P_THRESHOLD}" 'NR==1 || $12>thresh {print}' \
    ${OUTPUT_DIR}/combined_results.txt \
    > ${OUTPUT_DIR}/genome_wide_significant.txt

# Step 4: Create region-based summary (clump significant SNPs)
log "Step 4: Identifying independent loci..."

# Extract SNPs for clumping
awk 'NR>1 {print $3, $12}' ${OUTPUT_DIR}/genome_wide_significant.txt \
    > ${OUTPUT_DIR}/sig_snps_for_clumping.txt

# Note: Actual clumping would require PLINK with LD information
# This is a placeholder for the clumping step
cat > ${OUTPUT_DIR}/clumping_instructions.txt <<EOF
To identify independent loci, run PLINK clumping:

plink \\
    --bfile output/qc/step3_hwe \\
    --clump ${OUTPUT_DIR}/sig_snps_for_clumping.txt \\
    --clump-p1 ${P_THRESHOLD} \\
    --clump-r2 0.1 \\
    --clump-kb 250 \\
    --out ${OUTPUT_DIR}/clumped_results

This will identify independent genome-wide significant loci.
EOF

cat ${OUTPUT_DIR}/clumping_instructions.txt

# Step 5: Call visualization script
log "Step 5: Generating plots..."

if [ -f "utils/plot_manhattan.R" ]; then
    log "Generating Manhattan plot..."
    Rscript utils/plot_manhattan.R \
        ${OUTPUT_DIR}/combined_results.txt \
        ${OUTPUT_DIR}/manhattan_plot.png 2>/dev/null || log "Note: R plotting script may need dependencies"
else
    log "Note: Manhattan plot script not found at utils/plot_manhattan.R"
fi

# Step 6: Generate final report
log "Step 6: Generating final report..."

cat > ${OUTPUT_DIR}/final_report.txt <<EOF
REGENIE GWAS - Final Analysis Report
Generated: $(date)
=====================================

Analysis Files:
---------------
Combined results: ${OUTPUT_DIR}/combined_results.txt
Genome-wide significant: ${OUTPUT_DIR}/genome_wide_significant.txt
Summary statistics: ${OUTPUT_DIR}/summary_stats.txt

Key Findings:
-------------
EOF

cat ${OUTPUT_DIR}/summary_stats.txt >> ${OUTPUT_DIR}/final_report.txt

cat >> ${OUTPUT_DIR}/final_report.txt <<EOF

Output Files:
-------------
- Combined results: ${OUTPUT_DIR}/combined_results.txt
- Significant hits: ${OUTPUT_DIR}/genome_wide_significant.txt
- Summary stats: ${OUTPUT_DIR}/summary_stats.txt
- Manhattan plot: ${OUTPUT_DIR}/manhattan_plot.png (if R is available)

Recommendations:
----------------
1. Review genomic inflation factor (λ) - should be close to 1.0
2. Examine Manhattan and Q-Q plots for systematic issues
3. Perform LD clumping to identify independent loci
4. Annotate significant variants with gene information
5. Validate top hits in independent cohorts

For publication-ready results:
- Export significant loci to standard format
- Generate regional association plots for top loci
- Perform functional annotation
- Check for known GWAS catalog hits
EOF

cat ${OUTPUT_DIR}/final_report.txt

log "Results processing complete!"
log "Final report: ${OUTPUT_DIR}/final_report.txt"
log "All output files are in: ${OUTPUT_DIR}/"
