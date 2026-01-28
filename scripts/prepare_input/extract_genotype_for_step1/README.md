# Generate REGENIE Step 1 genotype input (CBMR in-house cohorts)

These scripts prepare REGENIE Step 1 genotype input from CBMR in-house cohort data (Michigan imputation *.dose.vcf.gz), following the current folder structure under /projects/magic-AUDIT/data/<COHORT>/genotype_imputed/.

Goal: extract genotyped variants from the imputation VCFs (instead of using all imputed variants), then create a single PLINK BED dataset for REGENIE Step 1.

Before you start: check how “genotyped sites are labeled

Different cohorts/VCF use differnt column to label genotype or imputed. 
Always inspect the VCF header first and confirm whether the cohort uses FILTER tags or INFO flags.

'FILTER="PASS;GENOTYPED"' : "Inter99" "Danfund" "Health06" "Health08" "Health10" "Vejle" "steno"
'FILTER="GENOTYPED"' : "glostrup"
'INFO/TYPED=1':  additionpro, lofus

# 1. Extract genotyped variants from chromosome 1–22

script: 1.extract_genotype_chr1.22.sh

# 2.  Extract chr23 (X) and combine chr1–23 into one dataset

script: 2.generate_genotype_input.sh

# 3. (optional) Build a common SNP set across multiple cohorts (T2DGGI)

For multi-cohort analyses (e.g. T2DGGI cohorts), you may want a shared set of variants present in all cohorts.

extract common SNP among T2DGGI cohort("Inter99" "Danfund" "Health06" "Health08" "Health10" "steno")

