# Slurm Usage Guide for Genetic QC Pipeline

This guide explains how to run the Genetic QC Pipeline on Slurm-managed HPC clusters.

## Overview

The pipeline provides several Slurm-compatible profiles:
- `slurm`: Basic Slurm execution
- `slurm_conda`: Slurm with Conda environment management
- `slurm_singularity`: Slurm with Singularity containers

## Quick Start

### 1. Set Environment Variables (Optional)

```bash
# Optional: Set Singularity cache directory
export NXF_SINGULARITY_CACHEDIR="$HOME/.singularity"
```

### 2. Basic Slurm Execution

```bash
nextflow run main.nf \
    -profile slurm \
    --data /path/to/your/plink_files \
    --name your_study_name
```

### 3. Using Custom Slurm Configuration

```bash
nextflow run main.nf \
    -profile slurm \
    -c conf/slurm.config \
    --data /path/to/your/plink_files \
    --name your_study_name
```

### 4. Submit as Slurm Job

```bash
# Edit scripts/submit_slurm.sh to match your cluster settings
sbatch scripts/submit_slurm.sh /path/to/your/plink_files your_study_name
```

## Configuration Options

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `NXF_SINGULARITY_CACHEDIR` | Singularity cache directory | `$HOME/.singularity` |
| `NXF_OPTS` | Nextflow JVM options | `-Xms1g -Xmx4g` |

### Slurm Profiles

#### `slurm`
Basic Slurm execution with default resource allocations:
- Uses system-installed software
- Suitable for clusters with pre-installed bioinformatics tools
- Uses 'cpu' partition by default

#### `slurm_conda`
Slurm execution with Conda environment:
- Automatically creates and manages Conda environment
- Best for reproducible environments
- Requires Conda/Mamba to be available on compute nodes

#### `slurm_singularity`
Slurm execution with Singularity containers:
- Uses containerized software stack
- Most portable and reproducible option
- Requires Singularity to be available on compute nodes

## Resource Requirements

### Default Resource Allocations

| Process Type | CPUs | Memory | Time | Partition |
|--------------|------|--------|------|-----------|
| Low | 2 | 6 GB | 4h | cpu |
| Medium | 6 | 36 GB | 8h | cpu |
| High | 12 | 72 GB | 16h | cpu |

### Process-Specific Resources

| Process | CPUs | Memory | Time | Notes |
|---------|------|--------|------|-------|
| MAF_FILTER | 4 | 8 GB | 2h | Fast filtering |
| ITERATIVE_MISSINGNESS | 6 | 16 GB | 6h | CPU intensive |
| IBD_CALCULATION | 8 | 32 GB | 12h | Memory intensive |
| ANCESTRY_PCA_WEIGHTS | 16 | 64 GB | 8h | High compute |

## Customization

### 1. Cluster-Specific Settings

Edit `conf/slurm.config` to match your cluster:

```groovy
process {
    // Change default partition
    queue = 'your_partition_name'

    // Add cluster-specific options
    clusterOptions = '--qos=your_qos --constraint=your_constraint'
}
```

### 2. Resource Limits

Adjust maximum resources in `conf/slurm.config`:

```groovy
params {
    max_memory_slurm = '512.GB'  // Your cluster's memory limit
    max_cpus_slurm = 64          // Your cluster's CPU limit
    max_time_slurm = '72.h'      // Your cluster's time limit
}
```

### 3. Partition Configuration

Configure partitions for different job types:

```groovy
params {
    partition_cpu = 'cpu'            // CPU jobs (default)
    // Add other partitions if available on your cluster
}
```

## Troubleshooting

### Common Issues

1. **Out of memory errors**
   - Increase memory allocation in `conf/slurm.config`
   - Check if your cluster has high-memory nodes available

2. **Time limit exceeded**
   - Increase time allocation for specific processes
   - Consider splitting large datasets

3. **Module loading issues**
   - Edit `scripts/submit_slurm.sh` to load required modules
   - Use container profiles (`slurm_singularity`) for better portability

### Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View job details
scontrol show job <job_id>

# Check resource usage
sacct -j <job_id> --format=JobID,JobName,MaxRSS,Elapsed,State

# Monitor Nextflow progress
tail -f .nextflow.log
```

### Performance Optimization

1. **Adjust queue size** in `conf/slurm.config`:
   ```groovy
   executor {
       queueSize = 100  // Increase for better throughput
   }
   ```

2. **Use local scratch space** for temporary files:
   ```groovy
   process {
       scratch = '/tmp'  // Use local SSD if available
   }
   ```

3. **Enable work directory cleanup**:
   ```bash
   nextflow run main.nf -profile slurm --cleanup
   ```

## Example Workflows

### Standard QC Analysis
```bash
nextflow run main.nf \
    -profile slurm_conda \
    --data /data/cohort/genotypes \
    --name COHORT_QC \
    --maf 0.01 \
    --geno 0.05 \
    --mind 0.05 \
    --hwe 0.0000000001
```

### High-throughput Analysis
```bash
nextflow run main.nf \
    -profile slurm \
    -c conf/slurm.config \
    --data /data/large_cohort/genotypes \
    --name LARGE_COHORT_QC \
    --max_memory 512.GB \
    --max_cpus 32
```

### Resume Failed Run
```bash
nextflow run main.nf \
    -profile slurm \
    -resume \
    --data /data/cohort/genotypes \
    --name COHORT_QC
```

## Support

For cluster-specific configuration help:
1. Check your cluster's documentation
2. Contact your system administrator
3. Refer to [Nextflow Slurm documentation](https://www.nextflow.io/docs/latest/executor.html#slurm)
