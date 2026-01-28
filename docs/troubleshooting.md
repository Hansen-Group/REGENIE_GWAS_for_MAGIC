# Troubleshooting Guide

This guide covers common issues and their solutions when running REGENIE GWAS analyses.

## Table of Contents
1. [Installation Issues](#installation-issues)
2. [Input File Problems](#input-file-problems)
3. [Memory and Performance Issues](#memory-and-performance-issues)
4. [Analysis Errors](#analysis-errors)
5. [Result Interpretation](#result-interpretation)

---

## Installation Issues

### Issue: REGENIE binary not found

**Error message:**
```
bash: regenie: command not found
```

**Solution:**
1. Verify REGENIE is installed:
   ```bash
   which regenie
   ```

2. If not installed, download and install (check https://github.com/rgcgithub/regenie/releases for latest version):
   ```bash
   # Version 3.2.5 shown as example:
   wget https://github.com/rgcgithub/regenie/releases/download/v3.2.5/regenie_v3.2.5.gz_x86_64_Linux.zip
   unzip regenie_v3.2.5.gz_x86_64_Linux.zip
   chmod +x regenie
   sudo mv regenie /usr/local/bin/
   ```

3. Verify installation:
   ```bash
   regenie --version
   ```

### Issue: Missing dependencies

**Error message:**
```
error while loading shared libraries: libboost_*.so
```

**Solution:**
Install required libraries:
```bash
# Ubuntu/Debian
sudo apt-get install libboost-all-dev

# CentOS/RHEL
sudo yum install boost-devel
```

---

## Input File Problems

### Issue: Phenotype-genotype ID mismatch

**Error message:**
```
ERROR: No individuals left after matching phenotype and genotype files
```

**Solution:**
1. Check ID format matches between files:
   ```bash
   # Check genotype IDs
   head -5 genotypes.fam
   
   # Check phenotype IDs
   head -5 phenotypes.txt
   ```

2. Ensure FID and IID match exactly (including spaces/tabs)

3. Common fixes:
   ```bash
   # If phenotype file has only IID, duplicate it as FID
   awk '{print $1,$1,$2,$3}' pheno_1col.txt > pheno_2col.txt
   ```

### Issue: Missing phenotype values

**Error message:**
```
WARNING: Phenotype has missing values for X individuals
```

**Solution:**
1. Missing values should be coded as `NA` or `-9`
2. Check for incorrect coding:
   ```bash
   grep -E "nan|NaN|null|NULL" phenotypes.txt
   ```

3. Replace incorrect missing codes:
   ```bash
   sed 's/NaN/NA/g' phenotypes.txt > phenotypes_fixed.txt
   ```

### Issue: Incorrect file format

**Error message:**
```
ERROR: Could not read phenotype file
```

**Solution:**
1. Ensure tab or space-delimited format
2. Check for special characters:
   ```bash
   file phenotypes.txt
   dos2unix phenotypes.txt  # If file is Windows-formatted
   ```

3. Verify header is present:
   ```bash
   head -1 phenotypes.txt
   ```

### Issue: Binary file version mismatch

**Error message:**
```
ERROR: Invalid magic number in .bed file
```

**Solution:**
Regenerate PLINK files with correct version:
```bash
plink --bfile old_format --make-bed --out new_format
```

---

## Memory and Performance Issues

### Issue: Out of memory error

**Error message:**
```
std::bad_alloc
terminate called after throwing an instance of 'std::bad_alloc'
```

**Solution:**
1. Enable low-memory mode:
   ```bash
   regenie --step 1 --lowmem --lowmem-prefix /tmp/regenie_tmp [other options]
   ```

2. Reduce block size:
   ```bash
   --bsize 500  # Instead of default 1000
   ```

3. Free up memory:
   ```bash
   # Check available memory
   free -h
   
   # Clear cache (if needed)
   sync && sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
   ```

### Issue: Process killed by system

**Error message:**
```
Killed
```

**Solution:**
This usually indicates out-of-memory (OOM) killer:

1. Monitor memory usage:
   ```bash
   # In another terminal
   watch -n 1 free -h
   ```

2. Split analysis by chromosome:
   ```bash
   # For Step 2, run one chromosome at a time
   for chr in {1..22}; do
       regenie --step 2 --chr ${chr} [other options]
   done
   ```

3. Use swap space:
   ```bash
   # Check swap
   swapon --show
   
   # Add swap if needed (as root)
   sudo fallocate -l 16G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

### Issue: Very slow Step 1

**Problem:** Step 1 taking days to complete

**Solution:**
1. Reduce number of variants for Step 1:
   ```bash
   # More aggressive LD pruning
   plink --bfile genotypes --indep-pairwise 1000 100 0.95 --out pruned
   # Should have ~100K-500K SNPs
   ```

2. Increase block size (if memory allows):
   ```bash
   --bsize 2000
   ```

3. Use more threads:
   ```bash
   --threads 16
   ```

---

## Analysis Errors

### Issue: No variants pass filters

**Error message:**
```
ERROR: No variants passed the MAC/INFO filters
```

**Solution:**
1. Check filter thresholds:
   ```bash
   # Relax filters if appropriate
   regenie --minMAC 5 --minINFO 0.3 [other options]
   ```

2. Verify variant quality in input data:
   ```bash
   plink --bfile genotypes --freq --out allele_freq
   head allele_freq.frq
   ```

### Issue: Singular matrix error

**Error message:**
```
ERROR: Singular matrix in ridge regression
```

**Solution:**
1. Check for duplicate variants:
   ```bash
   plink --bfile genotypes --list-duplicate-vars
   ```

2. Remove duplicates:
   ```bash
   plink --bfile genotypes --exclude duplicates.txt --make-bed --out cleaned
   ```

3. Check for perfectly correlated covariates:
   ```bash
   # In R, check correlation matrix
   covars <- read.table("covariates.txt", header=TRUE)
   cor(covars[,3:ncol(covars)])
   ```

### Issue: Convergence failure

**Error message:**
```
WARNING: Failed to converge for phenotype X
```

**Solution:**
1. For binary traits, check case-control balance:
   ```bash
   # Count cases and controls
   awk '{print $3}' phenotypes.txt | sort | uniq -c
   ```

2. If imbalanced (<1:10 ratio), ensure Firth is enabled:
   ```bash
   regenie --step 2 --firth --approx [other options]
   ```

3. Check for outliers in quantitative traits:
   ```r
   # In R
   pheno <- read.table("phenotypes.txt", header=TRUE)
   summary(pheno)
   boxplot(pheno$trait)
   ```

### Issue: Prediction file not found

**Error message:**
```
ERROR: Cannot find file in pred list
```

**Solution:**
1. Check Step 1 completed successfully:
   ```bash
   ls -lh step1_output_*.loco
   cat step1_output_pred.list
   ```

2. Verify paths in prediction list file:
   ```bash
   cat step1_output_pred.list
   # Should contain full or relative paths to .loco files
   ```

3. If files moved, update paths:
   ```bash
   # Update paths in pred list
   sed -i 's|old_path|new_path|g' step1_output_pred.list
   ```

---

## Result Interpretation

### Issue: Genomic inflation (λ > 1.1)

**Problem:** High genomic inflation factor indicates systematic bias

**Diagnosis:**
```r
# Calculate lambda in R
pvals <- results$P
chisq <- qchisq(1-pvals, 1)
lambda <- median(chisq, na.rm=TRUE) / 0.456
print(lambda)
```

**Solution:**
1. Check for population stratification:
   - Add more principal components as covariates
   - Verify ancestry of samples

2. Check for cryptic relatedness:
   ```bash
   # Calculate kinship matrix with KING
   king -b genotypes.bed --kinship
   ```

3. Check for genotyping artifacts:
   - Verify genotype calling quality
   - Check for batch effects

### Issue: No significant associations

**Problem:** No variants reach genome-wide significance

**Possible causes:**
1. **Insufficient power:**
   - Sample size too small
   - Effect sizes too small to detect

2. **Wrong phenotype:**
   - Check phenotype coding
   - Verify trait definition

3. **Poor quality control:**
   - Review QC thresholds
   - Check imputation quality

**Solutions:**
1. Examine suggestive hits (P < 1e-5)
2. Check Q-Q plot for evidence of signal
3. Consider gene-based or pathway analyses
4. Validate in external datasets

### Issue: Too many significant associations

**Problem:** Hundreds or thousands of genome-wide significant variants

**Possible causes:**
1. **Population stratification:** Check λ
2. **Sample duplication:** Check for duplicate individuals
3. **Batch effects:** Verify covariates

**Solutions:**
1. Verify QC steps were applied
2. Check Q-Q plot for deviation
3. Review covariate list
4. Consider conditional analysis

### Issue: Inconsistent results between runs

**Problem:** Different results with same input

**Possible causes:**
1. Random seed variation (for some algorithms)
2. Different input file versions
3. Different REGENIE versions

**Solutions:**
1. Set random seed (if applicable):
   ```bash
   --setl0 <seed> --setl1 <seed>
   ```

2. Document software versions:
   ```bash
   regenie --version > analysis_log.txt
   ```

3. Use version control for scripts and parameters

---

## Getting Help

If issues persist:

1. **Check REGENIE documentation:**
   - https://rgcgithub.github.io/regenie/

2. **Search GitHub issues:**
   - https://github.com/rgcgithub/regenie/issues

3. **Post a new issue with:**
   - Complete error message
   - REGENIE version
   - Input file descriptions
   - Command used
   - System information

4. **Contact repository maintainers:**
   - See README for contact information

---

## Quick Reference: Common Fixes

| Problem | Quick Fix |
|---------|-----------|
| Out of memory | Add `--lowmem --lowmem-prefix /tmp/rg` |
| Slow Step 1 | Reduce SNPs: `--extract pruned.snplist` |
| ID mismatch | Check FID/IID format in all files |
| No convergence | Enable Firth: `--firth --approx` |
| High lambda | Add more PCs as covariates |
| Missing predictions | Check `_pred.list` file paths |
| Binary version error | Regenerate with `plink --make-bed` |
