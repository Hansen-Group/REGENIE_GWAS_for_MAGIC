#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=T2DGGI%j.out
#SBATCH --error=T2DGGI%j.err
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
COMMON_DIR="target_path/common_genotype"
COMMON_SNPS="${COMMON_DIR}/common_snps.txt"

# Extract and sort SNPs from each cohort
for cohort in "${COHORT[@]}"; do
    awk '{print $2}' "${cohort}/genotype.bim" | sort > "${COMMON_DIR}/${cohort}_sorted_snps.txt"
done

# Initialize common SNP list with the first cohort
cp "${COMMON_DIR}/${COHORT[0]}_sorted_snps.txt" "$COMMON_SNPS"

# Iteratively find common SNPs across all cohorts
for cohort in "${COHORT[@]:1}"; do
    comm -12 "$COMMON_SNPS" "${COMMON_DIR}/${cohort}_sorted_snps.txt" > "${COMMON_SNPS}_tmp"
    mv "${COMMON_SNPS}_tmp" "$COMMON_SNPS"
done


plink --bfile T2DGGI/genotype_all \
      --extract ${COMMON_SNPS} \
      --make-bed \
      --out T2DGGI/genotype_common

echo "SNP extraction completed."