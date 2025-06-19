#!/bin/bash

# Test script for the Nextflow BioResource QC Pipeline
# This script performs basic validation of the pipeline

echo "=========================================="
echo "BioResource QC Pipeline - Test Script"
echo "=========================================="

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed or not in PATH"
    echo "Please install Nextflow: curl -s https://get.nextflow.io | bash"
    exit 1
fi

echo "✓ Nextflow found: $(nextflow -version | head -n1)"

# Check if required files exist
echo ""
echo "Checking pipeline files..."

required_files=(
    "main.nf"
    "nextflow.config"
    "modules/local/initial_qc.nf"
    "modules/local/ld_sex_het.nf"
    "modules/local/ibd_analysis.nf"
    "modules/local/ancestry_plink_ops.nf"
    "modules/local/ancestry_analysis.nf"
    "conda/environment.yml"
)

for file in "${required_files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "✓ $file"
    else
        echo "✗ $file (missing)"
        exit 1
    fi
done

# Check supporting scripts
echo ""
echo "Checking supporting scripts..."

supporting_scripts=(
    "bin/Iterative_Missingness.sh"
    "bin/highLDregions4bim_b38.awk"
    "bin/ancestry_identifier_r_only.r"
)

for script in "${supporting_scripts[@]}"; do
    if [[ -f "$script" ]]; then
        echo "✓ $script"
    else
        echo "✗ $script (missing)"
        exit 1
    fi
done

# Validate Nextflow syntax
echo ""
echo "Validating Nextflow syntax..."
if nextflow run main.nf --help > /dev/null 2>&1; then
    echo "✓ Nextflow syntax validation passed"
else
    echo "✗ Nextflow syntax validation failed"
    echo "Running syntax check..."
    nextflow run main.nf --help
    exit 1
fi

# Check if conda environment can be created
echo ""
echo "Checking conda environment..."
if command -v conda &> /dev/null; then
    echo "✓ Conda found: $(conda --version)"
    
    # Test environment creation (dry run)
    if conda env create -f conda/environment.yml --dry-run > /dev/null 2>&1; then
        echo "✓ Conda environment specification is valid"
    else
        echo "⚠ Conda environment specification may have issues"
    fi
else
    echo "⚠ Conda not found - skipping environment validation"
fi

echo ""
echo "=========================================="
echo "Pipeline validation completed successfully!"
echo "=========================================="
echo ""
echo "To run the pipeline with test data:"
echo "nextflow run main.nf \\"
echo "  --data /path/to/your/plink/files \\"
echo "  --name test_study \\"
echo "  --outdir results \\"
echo "  -profile conda"
echo ""
echo "For more information, see README.md"
