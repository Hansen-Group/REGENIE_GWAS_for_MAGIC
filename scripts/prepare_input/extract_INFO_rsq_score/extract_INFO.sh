#!/bin/bash
#SBATCH --output=extract_rsq_score%j.out
#SBATCH --error=extract_rsq_score%j.err
### Number of nodes
#SBATCH --cpus-per-task=2
### Memory
#SBATCH --mem-per-cpu=15gb
### Requesting time - format is <days>-<hours>:<minutes>:<seconds> (here, 48 hours)
#SBATCH --time=04:00:00


COHORT="cohort_name"  
output_dir="target_path/${COHORT}"
mkdir -p "$output_dir"

BASE_DIR="/projects/magic-AUDIT/data/${COHORT}/genotype_imputed/michigan.hrc.vcf" 

OUT_F="${output_dir}/extracted_INFO_Female.tsv"
OUT_M="${output_dir}/extracted_INFO_Male.tsv"

echo "Extract Rsq from INFO files for $COHORT"

# autosomes
for chr in {1..22}; do
  f="$BASE_DIR/chr${chr}.info.gz"
  [[ -f "$f" ]] || { echo "Missing $f, skipping..."; continue; }
  echo "Processing $f"
  zcat "$f" | awk 'NR>1{print $1"\t"$7}' >> "$OUT_F"
  zcat "$f" | awk 'NR>1{print $1"\t"$7}' >> "$OUT_M"
done

# chrX female
female_info="$BASE_DIR/chrX.no.auto_female.info.gz"
if [[ -f "$female_info" ]]; then
  echo "Processing $female_info"
  zcat "$female_info" | awk 'NR>1{sub(/^X:/,"23:",$1); print $1"\t"$7}' >> "$OUT_F"
fi

# chrX male
male_info="$BASE_DIR/chrX.no.auto_male.info.gz"
if [[ -f "$male_info" ]]; then
  echo "Processing $male_info"
  zcat "$male_info" | awk 'NR>1{sub(/^X:/,"23:",$1); print $1"\t"$7}' >> "$OUT_M"
fi