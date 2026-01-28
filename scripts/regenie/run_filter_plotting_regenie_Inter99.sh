#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=qc_plot_%j.out
#SBATCH --error=qc_plot_%j.err
### Number of nodes
#SBATCH --cpus-per-task=4
### Memory
#SBATCH --mem-per-cpu=15gb
### Requesting time - format is <days>-<hours>:<minutes>:<seconds> (here, 48 hours)
#SBATCH --time=08:00:00

module load perl
module load gsl/2.5
COHORT="Inter99"  

# Loop through all .regenie files in the current directory
for sex in Male Female; do

  echo "Filtering maf with $sex in $COHORT"
  info_file="/projects/magic-AUDIT/data/regenie/input/INFO/${COHORT}/extracted_INFO_${sex}.txt"
  workdir="/projects/magic-AUDIT/data/regenie/output_perTrait/${COHORT}/${sex}/"
  output_dir="${workdir}filter_maf"
  mkdir -p "$output_dir"

  for gwas_file in ${workdir}/*.regenie; 
  do
  # Extract the basename of the GWAS file and define the output file name
  base_name=$(basename "$gwas_file")
  output_file="${output_dir}/${base_name}_info03_maf001"
  # Create a temporary file 
  temp_file="${output_dir}/temp_${base_name}_with_maf.txt"

  # Step 1: Create the SNP column and calculate MAF without adding the MAF column to the final output
  awk 'BEGIN {OFS=" "} {
    if (NR==1) { 
      print "SNP", $0   # Print header line with new SNP column
    } else { 
      MAF = ($6 > 0.5) ? 1 - $6 : $6;  # Calculate MAF from A1FREQ
      if (MAF >= 0.01) {               # Filter MAF > 0.01
        print $1 ":" $2, $0            # Printll H   the line without adding the MAF column
      }
    }
  }' "$gwas_file" > "$temp_file"

  # Step 2: Filter rows based on Rsq score and print final output
  awk 'NR==FNR {rsq[$1]=$2; next} FNR==1 || rsq[$1] >= 0.3' "$info_file" "$temp_file" > "$output_file"

  # Step 3: Clean up temporary file
  rm "$temp_file"
  done
  # # Step 4: Plot Manhattan and QQ plots
  echo "Plotting Manhattan and QQ plots with $sex in $COHORT"
  module load --auto R/4.3.3
  Rscript magic_QC.R $COHORT $sex
  
done
