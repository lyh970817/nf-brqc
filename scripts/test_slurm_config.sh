#!/bin/bash

# Test script for Slurm configuration validation
# This script checks if the Slurm configuration is properly set up

set -euo pipefail

echo "=========================================="
echo "Slurm Configuration Test"
echo "=========================================="

# Check if Slurm is available
if ! command -v sinfo &> /dev/null; then
    echo "❌ ERROR: Slurm is not available on this system"
    echo "   Please ensure you're running this on a Slurm-managed cluster"
    exit 1
fi

echo "✅ Slurm is available"

# Check if cpu partition is available
echo ""
echo "Checking for 'cpu' partition:"
if sinfo -h -p cpu &> /dev/null; then
    echo "✅ 'cpu' partition is available"
else
    echo "⚠️  WARNING: 'cpu' partition not found"
    echo "   You may need to update the partition name in the configuration"
fi

# Check available partitions
echo ""
echo "Available Slurm partitions:"
sinfo -h -o "%P" | sort -u | while read partition; do
    echo "  - $partition"
done

# Check if Nextflow is available
echo ""
if ! command -v nextflow &> /dev/null; then
    echo "❌ ERROR: Nextflow is not available"
    echo "   Please install Nextflow: curl -s https://get.nextflow.io | bash"
    exit 1
fi

echo "✅ Nextflow is available: $(nextflow -version | head -n1)"

# Check required files
echo ""
echo "Checking pipeline files:"
required_files=(
    "main.nf"
    "nextflow.config"
    "conf/slurm.config"
    "scripts/submit_slurm.sh"
    "docs/SLURM_USAGE.md"
)

for file in "${required_files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "✅ $file"
    else
        echo "❌ $file (missing)"
        exit 1
    fi
done

# Test Nextflow configuration parsing
echo ""
echo "Testing Nextflow configuration parsing..."
if nextflow config -profile slurm > /dev/null 2>&1; then
    echo "✅ Slurm profile configuration is valid"
else
    echo "❌ ERROR: Slurm profile configuration has errors"
    echo "   Run 'nextflow config -profile slurm' to see details"
    exit 1
fi

# Check if conda/singularity profiles work
echo ""
echo "Testing additional profiles:"

if nextflow config -profile slurm_conda > /dev/null 2>&1; then
    echo "✅ slurm_conda profile is valid"
else
    echo "❌ slurm_conda profile has configuration errors"
fi

if nextflow config -profile slurm_singularity > /dev/null 2>&1; then
    echo "✅ slurm_singularity profile is valid"
else
    echo "❌ slurm_singularity profile has configuration errors"
fi

# Test job submission (dry run)
echo ""
echo "Testing dry run job submission..."
if command -v sbatch &> /dev/null; then
    echo "✅ sbatch command is available"
    
    # Create a simple test job
    cat > test_job.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=test-slurm-config
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:01:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

echo "Test job completed successfully"
EOF
    
    # Test job submission (but don't actually submit)
    if sbatch --test-only test_job.sh > /dev/null 2>&1; then
        echo "✅ Job submission test passed"
    else
        echo "⚠️  WARNING: Job submission test failed"
        echo "   This might be due to partition restrictions or cluster policies"
    fi
    
    rm -f test_job.sh
else
    echo "❌ sbatch command is not available"
fi

echo ""
echo "=========================================="
echo "Configuration Test Summary"
echo "=========================================="
echo "✅ Slurm support has been successfully added to your pipeline!"
echo ""
echo "Next steps:"
echo "1. Customize conf/slurm.config for your cluster if needed"
echo "2. Edit scripts/submit_slurm.sh with your cluster settings"
echo "3. Read docs/SLURM_USAGE.md for detailed usage instructions"
echo ""
echo "Test your setup with:"
echo "  nextflow run main.nf -profile slurm --help"
echo "=========================================="
