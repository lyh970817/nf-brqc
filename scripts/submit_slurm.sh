#!/bin/bash
#SBATCH --job-name=genetic-qc-pipeline
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# Genetic QC Pipeline Slurm Submission Script
# 
# This script submits the Nextflow genetic QC pipeline to a Slurm cluster.
# Customize the SBATCH parameters above according to your cluster requirements.
#
# Usage:
#   sbatch scripts/submit_slurm.sh /path/to/data study_name
#
# Arguments:
#   $1: Path to input PLINK files (without extension)
#   $2: Study name

set -euo pipefail

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: sbatch $0 <data_path> <study_name> [additional_nextflow_args]"
    echo "Example: sbatch $0 /data/mydata GLADv3 --maf 0.05"
    exit 1
fi

DATA_PATH="$1"
STUDY_NAME="$2"
shift 2
ADDITIONAL_ARGS="$@"

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules (customize for your cluster)
# module load nextflow/23.04.0
# module load singularity/3.8.0
# module load conda/4.12.0

# No Slurm account needed for this cluster

# Nextflow configuration
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-"$HOME/.singularity"}

# Pipeline directory
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "=========================================="
echo "Genetic QC Pipeline - Slurm Submission"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"
echo "Pipeline Directory: $PIPELINE_DIR"
echo "Data Path: $DATA_PATH"
echo "Study Name: $STUDY_NAME"
echo "Additional Args: $ADDITIONAL_ARGS"
echo "=========================================="

# Change to pipeline directory
cd "$PIPELINE_DIR"

# Run the pipeline with Slurm profile
echo "Starting Nextflow pipeline..."

nextflow run main.nf \
    -profile slurm \
    -c conf/slurm.config \
    --data "$DATA_PATH" \
    --name "$STUDY_NAME" \
    --outdir "results" \
    -with-report "results/pipeline_info/execution_report_${STUDY_NAME}.html" \
    -with-timeline "results/pipeline_info/execution_timeline_${STUDY_NAME}.html" \
    -with-trace "results/pipeline_info/execution_trace_${STUDY_NAME}.txt" \
    -with-dag "results/pipeline_info/pipeline_dag_${STUDY_NAME}.svg" \
    $ADDITIONAL_ARGS

echo "Pipeline completed with exit code: $?"
echo "Results available in: results/$STUDY_NAME"
echo "Pipeline reports available in: results/pipeline_info/"
