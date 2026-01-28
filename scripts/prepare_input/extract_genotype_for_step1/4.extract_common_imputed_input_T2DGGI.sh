#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=log/T2DGGI%j.out
#SBATCH --error=log/T2DGGI%j.err
### Number of nodes
#SBATCH --cpus-per-task=10
### Memory
#SBATCH --mem-per-cpu=15gb
### Requesting time - format is <days>-<hours>:<minutes>:<seconds> (here, 48 hours)
#SBATCH --time=04:00:00

module load perl
module load gsl/2.5
module load bcftools
module load plink/1.9.0
module load --auto htslib

# COHORT="T2DGGI"  
COHORT=("Inter99" "Danfund" "Health06" "Health08" "Health10" "steno")
COMMON_DIR="target_path/common_imputed_genotype"
COMMON_SNPS="${COMMON_DIR}/common_snps.txt"
mkdir -p "$COMMON_DIR"

# Extract and sort SNPs from each cohort
for cohort in "${COHORT[@]}"; do
    input_bim=$(find "/projects/magic-AUDIT/data/${cohort}/genotype_imputed/michigan.hrc.posids.plink19/" -name "*.dose.bim" | head -n 1)
    awk '{print $2}' "$input_bim" | sort > "${COMMON_DIR}/${cohort}_sorted_snps.txt"
    #awk '{print $2}' "/projects/magic-AUDIT/data/${cohort}/genotype_imputed/michigan.hrc.posids.plink19/*dose.bim" | sort > "${COMMON_DIR}/${cohort}_sorted_snps.txt"
done

# Initialize common SNP list with the first cohort
cp "${COMMON_DIR}/${COHORT[0]}_sorted_snps.txt" "$COMMON_SNPS"

# Iteratively find common SNPs across all cohorts
for cohort in "${COHORT[@]:1}"; do
    comm -12 "$COMMON_SNPS" "${COMMON_DIR}/${cohort}_sorted_snps.txt" > "${COMMON_SNPS}_tmp"
    mv "${COMMON_SNPS}_tmp" "$COMMON_SNPS"
done

