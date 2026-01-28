#!/bin/bash
### Output files (comment out the next 2 lines to get the job name used instead)
#SBATCH --output=log/regenie_%j.out
#SBATCH --error=log/regenie_%j.err
### Number of nodes
#SBATCH --cpus-per-task=15
### Memory
#SBATCH --mem-per-cpu=15gb
### Requesting time - format is <days>-<hours>:<minutes>:<seconds> (here, 48 hours)
#SBATCH --time=08:00:00
#SBATCH --array=1-16 # depends on number of traits to analyze

module load perl
module load gsl/2.5
module load bcftools
module load plink/1.9.0
module load libdeflate/1.18 htslib
module load regenie/3.2.5
 
## Define cohort and analyses file
## The array number in line 11 (array) needs to match number of analyses in analyses_file
COHORT="Inter99"  
cohort="inter99"
info_file_path="${COHORT}_traits.tsv"
output_dir="target_path/output_perTrait/"
## Read in trait_name, phenotype, and covariates from analysis_file ##
trait_name=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1; exit}' ${info_file_path})
phenotype=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2; exit}' ${info_file_path})
covariates=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3; exit}' ${info_file_path})

## Run regenie code twice, once for male once for female
for sex in Male Female
do
## only add covariate flags if there are covariants to include
if [[ $covariates == no_covariates ]]
then
   covariate_command=""
else
   covariate_command="--covarFile /projects/magic-AUDIT/data/regenie/input/phenotype/${COHORT}/${COHORT}_MAGIC_${sex}_covariates.txt --covarCol $covariates"
fi

echo "run regenie step 1 ......"
mkdir -p ${output_dir}/${COHORT}/${sex}
# Phenotypes adujst with standard covariates and 10 PCs
echo $trait_name
echo $sex

regenie \
 --step 1 \
 --bed /maps/projects/magic-AUDIT/data/regenie/input/genotype/${COHORT}/genotype \
 --phenoFile /projects/magic-AUDIT/data/regenie/input/phenotype/${COHORT}/${COHORT}_MAGIC_${sex}_phenotype.txt \
 --phenoCol ${phenotype} \
 $covariate_command\
 --remove /projects/glostrup-AUDIT/data/supportFiles/Failed_IBD_inter99_health_06_08_10_Danfund \
 --bsize 1000 \
 --qt \
 --write-samples \
 --lowmem \
 --lowmem-prefix tmp_reg_${COHORT}_${sex}_${trait_name} \
 --extract /projects/magic-AUDIT/data/regenie/input/genotype/${COHORT}/snps_pass_hwe_maf.snplist \
 --out ${output_dir}/${COHORT}/${sex}/${sex}_${trait_name} \
 --threads 13

regenie \
   --step 2 \
   --bed /projects/magic-AUDIT/data/${COHORT}/genotype_imputed/michigan.hrc.posids.rsids.plink19/${cohort}.dose \
    $covariate_command\
   --phenoFile /projects/magic-AUDIT/data/regenie/input/phenotype/${COHORT}/${COHORT}_MAGIC_${sex}_phenotype.txt \
   --phenoCol ${phenotype} \
   --remove /projects/glostrup-AUDIT/data/supportFiles/Failed_IBD_inter99_health_06_08_10_Danfund \
   --bsize 400 \
   --qt \
   --write-samples \
   --pred ${output_dir}/${COHORT}/${sex}/${sex}_${trait_name}_pred.list \
   --out ${output_dir}/${COHORT}/${sex}/${sex}_${trait_name} \
   --threads 13

done
