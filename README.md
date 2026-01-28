# REGENIE GWAS for MAGIC Consortium

This repository archives the scripts and workflows used for running Genome-Wide Association Studies (GWAS) for the MAGIC consortium using REGENIE in 2024.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Contributing](#contributing)
- [Citation](#citation)
- [Contact](#contact)

## Overview

REGENIE is a fast and efficient tool for analyzing large-scale genetic data in genome-wide association studies. This repository contains:
- Production-ready scripts for running REGENIE analyses
- Data preparation and quality control pipelines
- Result processing and visualization tools
- Example configurations and documentation

### What is REGENIE?

REGENIE is a C++ program for whole genome regression modelling of large genome-wide association studies. It is designed to handle:
- Binary and quantitative traits
- Multiple testing correction
- Inclusion of covariates
- Analysis of rare variants

### About MAGIC

The Meta-Analyses of Glucose and Insulin-related traits Consortium (MAGIC) investigates the genetic regulation of glucose and insulin-related traits.

## Repository Structure

```
REGENIE_GWAS_for_MAGIC/
├── README.md                     # This file
├── scripts/                      # Main analysis scripts
│   ├── 01_data_preparation.sh   # Data QC and preparation
│   ├── 02_regenie_step1.sh      # REGENIE Step 1: Whole genome regression
│   ├── 03_regenie_step2.sh      # REGENIE Step 2: Association testing
│   └── 04_process_results.sh    # Results processing and filtering
├── config/                       # Configuration files
│   ├── regenie_config.yaml      # REGENIE parameters
│   └── phenotype_example.txt    # Example phenotype file
├── utils/                        # Utility scripts
│   ├── summarize_results.py     # Summarize GWAS results
│   └── plot_manhattan.R         # Create Manhattan plots
└── docs/                         # Additional documentation
    ├── REGENIE_guide.md         # Detailed REGENIE usage guide
    └── troubleshooting.md       # Common issues and solutions
```

## Requirements

### Software Dependencies
- **REGENIE** (v3.0 or higher): [Installation guide](https://rgcgithub.github.io/regenie/)
- **PLINK** (v1.9 or v2.0): For genetic data processing
- **R** (v4.0+): For visualization and result processing
  - Required R packages: `ggplot2`, `data.table`, `dplyr`
- **Python** (v3.7+): For utility scripts
  - Required packages: `pandas`, `numpy`, `scipy`

### Hardware Requirements
- Minimum 32 GB RAM (64 GB+ recommended for large datasets)
- Multi-core CPU (8+ cores recommended)
- Sufficient storage for genotype data and results

## Installation

1. **Clone this repository:**
   ```bash
   git clone https://github.com/Hansen-Group/REGENIE_GWAS_for_MAGIC.git
   cd REGENIE_GWAS_for_MAGIC
   ```

2. **Install REGENIE:**
   ```bash
   # Download and install REGENIE
   wget https://github.com/rgcgithub/regenie/releases/download/v3.2.5/regenie_v3.2.5.gz_x86_64_Linux.zip
   unzip regenie_v3.2.5.gz_x86_64_Linux.zip
   chmod +x regenie
   # Move to a directory in your PATH
   sudo mv regenie /usr/local/bin/
   ```

3. **Install R packages:**
   ```R
   install.packages(c("ggplot2", "data.table", "dplyr"))
   ```

4. **Install Python packages:**
   ```bash
   pip install pandas numpy scipy
   ```

## Usage

### Quick Start

1. **Prepare your data:**
   - Genotype files in PLINK format (`.bed`, `.bim`, `.fam`)
   - Phenotype file with individual IDs and trait values
   - Covariate file (optional)

2. **Edit configuration:**
   ```bash
   # Update paths and parameters in config/regenie_config.yaml
   nano config/regenie_config.yaml
   ```

3. **Run the pipeline:**
   ```bash
   # Step 1: Data preparation
   bash scripts/01_data_preparation.sh
   
   # Step 2: REGENIE Step 1 (whole genome regression)
   bash scripts/02_regenie_step1.sh
   
   # Step 3: REGENIE Step 2 (association testing)
   bash scripts/03_regenie_step2.sh
   
   # Step 4: Process results
   bash scripts/04_process_results.sh
   ```

## Workflow

### Step 1: Data Preparation and QC
- Quality control filtering
- Missing data handling
- Ancestry verification
- Covariate preparation

### Step 2: REGENIE Step 1
- Whole genome regression model
- Builds predictive model using high-quality variants
- Computes LOCO (Leave-One-Chromosome-Out) predictions

### Step 3: REGENIE Step 2
- Association testing at all variants
- Uses predictions from Step 1 to control for population structure
- Outputs summary statistics for each variant

### Step 4: Results Processing
- Filter results by significance threshold
- Generate Manhattan and Q-Q plots
- Summarize top hits
- Export results in standard format

## Contributing

This is an archived repository documenting the 2024 GWAS workflow. For questions or suggestions, please contact the repository maintainers.

## Citation

If you use REGENIE, please cite:
```
Mbatchou et al. (2021) "Computationally efficient whole-genome regression for quantitative and binary traits." 
Nature Genetics, 53, 1097–1103. https://doi.org/10.1038/s41588-021-00870-7
```

For MAGIC consortium publications, please refer to:
https://magicinvestigators.org/publications/

## Contact

For questions about this repository:
- Repository: [REGENIE_GWAS_for_MAGIC](https://github.com/Hansen-Group/REGENIE_GWAS_for_MAGIC)
- MAGIC Consortium: https://magicinvestigators.org/

---

**Note:** This repository contains archived scripts from 2024. Always ensure you are using the latest version of REGENIE and associated tools for new analyses.
