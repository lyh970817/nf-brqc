# BioResource Genetic QC Pipeline

A Nextflow implementation of the BioResource Genetic QC Pipeline, originally written as a bash script (`fullpipe.sh`).

## Overview

This pipeline performs quality control (QC) on genetic data in PLINK format. It includes the following steps:

1. Initial data check and missingness filtering
2. Hardy-Weinberg equilibrium (HWE) filtering
3. Sex check and heterozygosity analysis
4. Identity-by-descent (IBD) analysis
5. Ancestry analysis using 1000 Genomes reference data

## Requirements

- Nextflow (>=22.10.0)
- Docker or Conda
- PLINK 1.9
- PLINK 2.0
- R (>=4.3.0) with required packages

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/bioresource-qc-pipeline.git
   cd bioresource-qc-pipeline
   ```

2. Install Nextflow:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. Choose your execution method:
   - Docker: `docker build -t bioresource-qc:latest .`
   - Conda: No additional steps required, the pipeline will create the environment automatically

## Usage

### Basic usage:

```bash
nextflow run main.nf --data /path/to/plink_files_prefix --name study_name
```

### Required parameters:

- `--data`: Path to the input PLINK files (without the .bed/.bim/.fam extension)
- `--name`: Study name (e.g., GLADv3, GLADv3_eur)

### Optional parameters:

- `--outdir`: Output directory (default: results)
- `--maf`: Minor allele frequency threshold (default: 0.01)
- `--SNP_CR`: SNP call rate threshold (default: 95)
- `--Sample_CR`: Sample call rate threshold (default: 95)
- `--geno`: Missing variant threshold (default: 0.05)
- `--mind`: Missing individual threshold (default: 0.05)
- `--hwe`: Hardy-Weinberg equilibrium threshold (default: 0.0000000001)
- `--ibd`: Pairwise identical-by-descent threshold (default: 0.1875)
- `--ind_ibd`: Individual IBD outlier, standard deviations (default: 3)

### Execution profiles:

- `conda`: Uses Conda for dependency management
  ```bash
  nextflow run main.nf -profile conda --data /path/to/plink_files_prefix --name study_name
  ```

- `docker`: Uses Docker for dependency management
  ```bash
  nextflow run main.nf -profile docker --data /path/to/plink_files_prefix --name study_name
  ```

- `singularity`: Uses Singularity for dependency management
  ```bash
  nextflow run main.nf -profile singularity --data /path/to/plink_files_prefix --name study_name
  ```

## Output

The pipeline creates a directory structure with the following components:

```
results/
└── study_name_timestamp/
    ├── QC_pipeline_review_documents_study_name/
    │   ├── Het_check_hist_study_name.pdf
    │   ├── HetFhat2_SexF_study_name.pdf
    │   ├── IBD_hist_study_name.pdf
    │   ├── Iterative_missingness_table_study_name.txt
    │   ├── callrate_pihat_over_9_study_name.txt
    │   ├── pihat_btw_4_6_study_name.txt
    │   ├── pihat_over_9_study_name.txt
    │   ├── sexmismatch_study_name.txt
    │   └── ... (other QC plots and tables)
    ├── study_name_finalreport.log
    ├── study_name_ancestry_identifier.log
    ├── study_name.PCs_plot_super_pop.png
    └── ... (other output files)
```

## Reference Data

The pipeline requires 1000 Genomes reference data for ancestry analysis. These files should be placed in the `ref/1kg/` directory:

- `1KG_Phase3.chr*.bed/bim/fam`: 1000 Genomes PLINK files per chromosome
- `ref_super_pop_keep_list.txt`: Reference population scale file
- `ref_pop_dat_reduced.txt`: Population data file

## Key Changes from Original Pipeline

This Nextflow implementation includes several important improvements over the original bash script:

1. **Separated PLINK Operations**: All `system()` calls to PLINK/PLINK2 in R scripts (especially in `ancestry_identifier.r`) have been moved to dedicated Nextflow processes, following best practices for workflow management.

2. **Modular Architecture**: Each QC step is now implemented as a separate, reusable Nextflow module with proper input/output channel handling.

3. **Improved Resource Management**: Automatic CPU and memory allocation based on process requirements, with configurable resource limits.

4. **Enhanced Error Handling**: Built-in retry mechanisms and better error reporting.

5. **Containerization Support**: Full support for Docker, Singularity, and Conda environments for reproducible execution.

6. **Execution Tracking**: Comprehensive execution reports, timelines, and resource usage tracking.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- Original bash script developed by the BioResource team
- Nextflow implementation by [Your Name]
