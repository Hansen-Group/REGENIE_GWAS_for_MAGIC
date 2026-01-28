#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=log/1.extract_genotype_%j.out
#SBATCH --error=log/1.extract_genotype_%j.err
### Number of nodes
#SBATCH --cpus-per-task=10
### Memory
#SBATCH --mem-per-cpu=15gb
### Requesting time - format is <days>-<hours>:<minutes>:<seconds> (here, 48 hours)
#SBATCH --time=04:00:00
#SBATCH --array 1-22 # for chr 1-22

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module load perl
module load gsl/2.5
module load bcftools
module load plink/1.9.0
module load libdeflate/1.18 htslib

COHORT="cohort_name" 
output_dir="target_output_path/${COHORT}"
genotype_filter='FILTER="PASS;GENOTYPED"'  # or INFO/TYPED=1 or 'FILTER="GENOTYPED"' depending on cohort
mkdir -p "$output_dir"
echo "extract genotyped variants in chr1-22 for cohort: $COHORT"

input_vcf="/projects/magic-AUDIT/data/${COHORT}/genotype_imputed/michigan.hrc.posids.rsids.vcf/chr"${SLURM_ARRAY_TASK_ID}".dose.vcf.gz"
filtered_vcf="${output_dir}/genotyped_only.chr"${SLURM_ARRAY_TASK_ID}".vcf"

bcftools filter -i "$genotype_filter" "$input_vcf" -o "$filtered_vcf" --threads 10
