#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=log/2.generate_genotype_input_%j.out
#SBATCH --error=log/%j.err
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
module load libdeflate/1.18 htslib

COHORT="cohort_name"  
output_dir="target_output_path/${COHORT}"
mkdir -p "$output_dir"
genotype_filter='FILTER="PASS;GENOTYPED"'  # or INFO/TYPED=1 or 'FILTER="GENOTYPED"' depending on cohort

# 1. extract chr 23

echo "extract genotyped variants in chr23 for cohort: $COHORT"
input_vcf="/projects/magic-AUDIT/data/${COHORT}/genotype_imputed/michigan.hrc.posids.rsids.vcf/chr23.dose.vcf.gz"
filtered_vcf="${output_dir}/genotyped_only.chr23.vcf"
bcftools filter -i "$genotype_filter" -o "$filtered_vcf" --threads 10 "$input_vcf"

# 2. concatenate chr 1-23

echo "concatenating chr 1-23 into a single vcf.gz"
bcftools concat ${output_dir}/genotyped_only.chr{1..23}.vcf -Oz -o ${output_dir}/genotyped.vcf.gz --threads 10

echo "converting vcf to plink format"

plink \
--vcf "${output_dir}/genotyped.vcf.gz" \
--double-id --keep-allele-order --make-bed \
--out "${output_dir}/genotype" --memory 60000 --threads 4