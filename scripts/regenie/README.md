# Run REGENIE (Step 1 + Step 2)

This pipeline runs REGENIE GWAS in parallel across multiple traits using a SLURM job array.
Each SLURM array task corresponds to one row in the trait definition file (*_traits.tsv).

## Parallel runs and missing phenotypes

The scripts support running multiple traits in parallel (recommended: one trait per SLURM array job).

If you attempt to run multiple phenotypes/GWAS in a single REGENIE run, make sure there are no missing values in the phenotype column(s). Otherwise REGENIE may impute missing phenotypes for some individuals, which can affect results.

### Recommended practice
	•	Run each trait separately (one trait per array task), or
	•	Pre-filter samples to remove individuals with missing phenotypes before running REGENIE.

## Sex-stratified analyses (MAGIC requirement)

Due to MAGIC consortium requirements, analyses are run separately in males and females.
Make sure that sex-stratified phenotype/covariate files include only the corresponding samples and contain no missing phenotype values.



# Inputs

For generating genotype variants, phenotype files, and covariate files, refer to scripts under ../prepare_input.

## 1) Trait definition file (*_traits.tsv)

A tab-separated file where each row defines one analysis.

### Minimum required columns:
	•	trait_name: output name for the analysis
	•	phenotype: phenotype column name in the phenotype file
	•	covariates: covariate column names (or no_covariates)

### Example

trait_name    phenotype_col    covariates
BMI           bmi              age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
TG            triglycerides    age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

See Inter99_traits.tsv for a working example.


## 2) Genotype input

	• Step 1 genotype: genotyped variant BED (smaller set used to fit ridge model and LOCO predictors)
	• Step 2 genotype: imputed genotype dataset (BED/BGEN depending on pipeline choice)


## 3) Phenotype and covariate files

### Sex-stratified phenotype/covariate files are typically used, e.g.:
	•	${COHORT}_MAGIC_Male_phenotype.txt, ${COHORT}_MAGIC_Male_covariates.txt
	•	${COHORT}_MAGIC_Female_phenotype.txt, ${COHORT}_MAGIC_Female_covariates.txt



# Sample exclusions and relatedness filtering (IMPORTANT)

### Several CBMR in-house cohorts include:
	•	participants who have requested withdrawal (must be excluded), and/or
	•	related individuals (should be excluded for standard GWAS analyses)

## 1) Read cohort-specific README files

Cohort folders often contain README documentation (provided by DAP) describing known issues and mandatory exclusions.

### Example (Inter99):
	•	Some participants requested withdrawal through Glostrup Hospital and must be excluded.
	•	These individuals are listed in withdrawals.ids under the supportFiles directory.

### Additional cohort-specific documentation may exist, e.g.:
	•	/datasets/inter99-AUDIT/Readme_before_analysis_IMPORTANT.txt

Always review these documents before analysis.


## 2) Remove related individuals 

If your analysis involves multiple cohorts from the larger Glostrup cohort collection, you may need to remove related individuals.

### Relatedness exclusion lists are available under:
	•	/projects/glostrup-AUDIT/data/supportFiles/

### Example:
	•	Failed_IBD_inter99_health_06 (or similar cohort-specific IBD exclusion files)

### Recommended exclusion practice

When running REGENIE, exclusions are typically applied via:
	•	--remove <file_with_FID_IID>

# Notes for multi-cohort workflows

	•	Some cohorts share related individuals across datasets; ensure you use the correct relatedness exclusion list when combining cohorts.
	•	Always confirm chromosome naming conventions (chr23 vs chrX) and sex handling if X chromosome variants are included.
	•	For sex-stratified runs, ensure phenotype/covariate files contain only samples for that sex and have no missing phenotype values.
