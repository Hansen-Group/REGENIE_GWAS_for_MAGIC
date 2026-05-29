# REGENIE_GWAS_for_MAGIC

This repository archives the scripts used to run GWAS for the MAGIC consortium (2024) using REGENIE on CBMR internal HPC (Esrum; maintained by DAP).


Data location (internal):
All source data and full outputs are stored under:
/projects/magic-AUDIT/data/regenie/

This repo contains scripts + lightweight metadata only (no individual-level genotype/phenotype data).


# Overview

The pipeline follows a standard REGENIE workflow:

## 1.	Prepare inputs

	•	Extract / convert genotype inputs for REGENIE Step 1 (typed/genotyped variant set)\
	•	Prepare imputed genotype input for REGENIE Step 2 (BED or BGEN depending on cohort setup)
	•	Prepare phenotype and covariate files (MAGIC requirement: sex-stratified files)
	•	Extract imputation quality metrics (Rsq/INFO) when needed for downstream filtering

## 2.	Run GWAS (REGENIE)

	•	Step 1: ridge regression + LOCO predictors
	•	Step 2: association testing on imputed variants
	•	Run Male and Female separately (MAGIC requirement)

## 3.	Post-processing

	•	Filter summary statistics by MAF and imputation quality (Rsq)
	•	Generate QC plots (Manhattan/QQ)


# Quick start (typical workflow)

## 1) Prepare inputs

### Run scripts under:
	•	scripts/prepare_input/

### Typical steps include:
	•	Extract genotyped/typed variants from Michigan imputation VCFs (chr1–22 and chrX/chr23) for Step 1
	•	Concatenate per-chromosome files and convert to PLINK BED
	•	Extract imputation quality metrics (Rsq/INFO) from Minimac/Michigan *.info.gz
	•	(Optional) generate common SNP set across cohorts for multi-cohort consistency

See scripts/prepare_input/README.md (or step-level READMEs) for details.

## 2) Run REGENIE

### Run scripts under:
	•	scripts/regenie/

### A typical REGENIE run is:
	•	Step 1 using typed BED genotype
	•	Step 2 using imputed genotype (BED or BGEN depending on cohort)
	•	Sex-stratified analyses: Male and Female are run separately

Most runs are implemented as SLURM job arrays, where each array task corresponds to one trait from a trait definition file (*_traits.tsv).

## 3) Post-process and qc plot

### After REGENIE completes, post-processing scripts may:
	•	Filter results by MAF (e.g. MAF >= 0.01)
	•	Filter by imputation quality (e.g. Rsq >= 0.3)
	•	Generate Manhattan and QQ plots using Rscript

# Important notes and exclusions

## Sex-stratified requirement (MAGIC)

### MAGIC requires sex-stratified analysis for this project, so phenotypes/covariates are prepared separately for:
	•	Male
	•	Female

### Ensure:
	•	No missing phenotype values in the phenotype files used for each run
	•	Sample IDs match between genotype and phenotype/covariate files

## Sample exclusions: withdrawals and relatedness

### Some cohorts contain:
	•	Participants who requested withdrawal (must be excluded)
	•	Related individuals (should be excluded for standard GWAS)


# Contact 

This repo is an archive of the 2024 MAGIC GWAS scripts and is primarily intended for internal reproducibility and reference within CBMR Hansen group.

Contact people

Shanhsan He <shanshan.he@sund.ku.dk>
