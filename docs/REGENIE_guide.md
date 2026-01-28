# REGENIE Usage Guide

This guide provides detailed information on using REGENIE for genome-wide association studies.

## Table of Contents
1. [Introduction to REGENIE](#introduction-to-regenie)
2. [Two-Step Procedure](#two-step-procedure)
3. [Input File Formats](#input-file-formats)
4. [Step-by-Step Tutorial](#step-by-step-tutorial)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

## Introduction to REGENIE

REGENIE is a C++ program for whole genome regression modelling of large genome-wide association studies developed by Mbatchou et al. (2021). It implements a fast and efficient method for association testing that:

- Scales to millions of samples and variants
- Controls for population structure and relatedness
- Handles both binary and quantitative traits
- Supports Firth regression for rare variants
- Uses a two-step procedure for computational efficiency

### Key Features

- **Speed**: Analyzes millions of variants in hours, not days
- **Memory efficiency**: Uses block processing and low-memory modes
- **Statistical power**: LOCO (Leave-One-Chromosome-Out) predictions
- **Flexibility**: Supports multiple phenotypes and covariates
- **Robustness**: Handles case-control imbalance with Firth regression

## Two-Step Procedure

REGENIE uses a two-step approach:

### Step 1: Whole Genome Regression Model
- Uses a subset of high-quality, directly genotyped variants
- Fits a whole genome regression model
- Computes LOCO predictions for each chromosome
- Creates predictive models to control for polygenic effects

**Purpose**: Build a genetic predictor to control for population structure and polygenic effects

### Step 2: Association Testing
- Uses predictions from Step 1 as offset
- Tests all variants (including imputed variants)
- Outputs summary statistics for each variant
- Can apply Firth regression for binary traits

**Purpose**: Test association of each variant with the phenotype

## Input File Formats

### 1. Genotype Files

REGENIE accepts genotype data in several formats:

#### PLINK Binary Format (.bed/.bim/.fam)
```bash
--bed <prefix>
```

Most common format for Step 1 and Step 2.

#### BGEN Format
```bash
--bgen <file>
--sample <file>  # Sample file
```

Commonly used for imputed data in Step 2.

#### PGEN Format (PLINK 2.0)
```bash
--pgen <prefix>
```

Modern PLINK format, more efficient for large datasets.

### 2. Phenotype File

Tab or space-delimited file with header:

```
FID IID phenotype1 phenotype2 ...
FAM001 IND001 95.5 0
FAM002 IND002 102.3 1
```

**Requirements**:
- First two columns: FID (Family ID) and IID (Individual ID)
- Missing values: NA or -9
- Binary traits: 0 (control), 1 (case)
- Quantitative traits: measured values

### 3. Covariate File

Tab or space-delimited file with header:

```
FID IID age sex PC1 PC2 PC3
FAM001 IND001 45 1 0.023 -0.015 0.008
FAM002 IND002 52 0 0.019 -0.012 0.006
```

**Common covariates**:
- Age
- Sex (0/1)
- Principal components (PC1-PC10)
- Assessment center (for multi-center studies)
- Genotyping batch

## Step-by-Step Tutorial

### Prerequisites

1. **Quality-controlled genotype data**
   - Remove low-quality variants
   - Remove samples with high missingness
   - Filter by MAF and HWE

2. **Prepared phenotype and covariate files**
   - Match individual IDs with genotype data
   - Handle missing values
   - Transform quantitative traits if needed

### Running REGENIE

#### Step 1: Whole Genome Regression

```bash
regenie \
  --step 1 \
  --bed data/genotypes \
  --phenoFile data/phenotypes.txt \
  --covarFile data/covariates.txt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --out results/step1
```

**Key parameters**:
- `--bsize`: Block size (larger = more memory, faster)
- `--lowmem`: Enable low-memory mode for large datasets
- `--threads`: Number of CPU threads

**Output**:
- `results/step1_pred.list`: List of LOCO prediction files
- `results/step1_*.loco`: LOCO predictions for each chromosome

#### Step 2: Association Testing

```bash
regenie \
  --step 2 \
  --bed data/genotypes \
  --phenoFile data/phenotypes.txt \
  --covarFile data/covariates.txt \
  --pred results/step1_pred.list \
  --bsize 400 \
  --out results/step2
```

**Key parameters**:
- `--pred`: Predictions from Step 1
- `--firth`: Enable Firth regression (binary traits)
- `--approx`: Use approximate Firth (faster)
- `--minMAC`: Minimum minor allele count

**Output**:
- `results/step2_<phenotype>.regenie`: Association results

### Output Format

REGENIE Step 2 output columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| GENPOS | Genomic position |
| ID | Variant identifier |
| ALLELE0 | Reference allele |
| ALLELE1 | Alternate allele (tested) |
| A1FREQ | Allele 1 frequency |
| N | Sample size |
| TEST | Test performed |
| BETA | Effect size estimate |
| SE | Standard error |
| CHISQ | Chi-square statistic |
| LOG10P | -log10(p-value) |

## Best Practices

### 1. Variant Selection for Step 1

Use high-quality, directly genotyped variants:
- MAF > 1% (or higher for small samples)
- Genotyping rate > 95%
- HWE p-value > 1e-10
- LD-pruned (r² < 0.9)

```bash
# Example LD pruning with PLINK
plink --bfile genotypes \
      --indep-pairwise 1000 100 0.9 \
      --out pruned_snps
```

### 2. Sample Size Considerations

| Sample Size | Recommendations |
|-------------|-----------------|
| < 50K | Standard settings |
| 50K - 200K | Use `--lowmem`, increase `--bsize` |
| > 200K | Use `--lowmem`, consider splitting by chromosome |

### 3. Binary Trait Analysis

For case-control studies:

```bash
regenie \
  --step 2 \
  --bt \  # Binary trait flag
  --firth \  # Firth regression
  --approx \  # Approximate Firth
  --pThresh 0.01 \  # Use Firth for p < 0.01
  [other options]
```

**When to use Firth**:
- Imbalanced case-control ratio
- Rare variants (MAC < 100)
- Prevents inflation of test statistics

### 4. Quantitative Trait Analysis

For continuous phenotypes:

```bash
regenie \
  --step 2 \
  --qt \  # Quantitative trait flag (default)
  [other options]
```

**Recommendations**:
- Consider inverse normal transformation
- Adjust for covariates affecting trait distribution
- Check for outliers

### 5. Multiple Phenotypes

REGENIE can analyze multiple phenotypes simultaneously:

```bash
# Phenotype file with multiple columns
FID IID glucose insulin hba1c
```

REGENIE will output separate results files for each phenotype.

### 6. Computational Optimization

**Memory optimization**:
```bash
--lowmem --lowmem-prefix /tmp/regenie_tmp
```

**Speed optimization**:
```bash
--threads 16  # Use more CPU cores
--bsize 1000  # Larger blocks (more memory)
```

**For very large datasets**:
- Split Step 2 by chromosome
- Use distributed computing
- Consider GPU acceleration (if available)

## Common Use Cases

### Case 1: Standard Quantitative Trait GWAS

```bash
# Step 1
regenie --step 1 --bed genotypes_qc --phenoFile pheno.txt \
        --covarFile covars.txt --bsize 1000 --out step1_out

# Step 2
regenie --step 2 --bed genotypes_qc --phenoFile pheno.txt \
        --covarFile covars.txt --pred step1_out_pred.list \
        --bsize 400 --out step2_out
```

### Case 2: Binary Trait with Firth Regression

```bash
# Step 2 only (with Firth)
regenie --step 2 --bed genotypes_qc --bt \
        --phenoFile pheno.txt --covarFile covars.txt \
        --pred step1_out_pred.list --firth --approx \
        --out step2_binary
```

### Case 3: Imputed Data Analysis

```bash
# Step 2 with BGEN imputed data
regenie --step 2 --bgen imputed_chr{1:22}.bgen \
        --sample imputed.sample \
        --phenoFile pheno.txt --covarFile covars.txt \
        --pred step1_out_pred.list \
        --minINFO 0.6 --minMAC 5 \
        --out step2_imputed
```

## Performance Tips

1. **Use LD-pruned SNPs for Step 1**: Reduces computation without loss of accuracy
2. **Adjust block size**: Balance between memory and speed
3. **Enable low-memory mode**: Essential for large datasets
4. **Use temporary storage**: `--lowmem-prefix` on fast local disk
5. **Parallel processing**: Use `--threads` to utilize multiple cores
6. **Split by chromosome**: For Step 2 on very large datasets

## Interpreting Results

### Genomic Inflation Factor (λ)

Calculate λ to check for inflation:

```r
# In R
pvalues <- results$P
chisq <- qchisq(1 - pvalues, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
```

**Interpretation**:
- λ ≈ 1.0: Good control for population structure
- λ < 1.05: Acceptable
- λ > 1.05: May indicate stratification or other issues

### Significance Thresholds

- **Genome-wide significance**: P < 5×10⁻⁸
- **Suggestive**: P < 1×10⁻⁵
- **Bonferroni correction**: P < 0.05 / (number of variants)

## References

1. Mbatchou et al. (2021). "Computationally efficient whole-genome regression for quantitative and binary traits." Nature Genetics.
2. REGENIE documentation: https://rgcgithub.github.io/regenie/
3. REGENIE GitHub: https://github.com/rgcgithub/regenie
