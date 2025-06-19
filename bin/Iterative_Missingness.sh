#!/bin/bash

# Iterative Missingness script
# This script performs iterative missingness filtering on PLINK data
# Usage: Iterative_Missingness.sh <start> <end> <step> <prefix>

start=$1
end=$2
step=$3
prefix=$4

# Initialize variables
current_snp=$start
current_sample=$start

# Perform iterative missingness filtering
while [ $current_snp -le $end ]; do
    # Convert percentage to decimal for PLINK
    geno_threshold=$(echo "scale=2; 1-($current_snp/100)" | bc)
    
    # First iteration uses the original file
    if [ $current_snp -eq $start ]; then
        plink --bfile $prefix --allow-extra-chr --geno $geno_threshold --make-bed --out ${prefix}.common_SNP${current_snp}
    else
        # Use the previous iteration's output
        plink --bfile ${prefix}.common_sample${prev_sample}.SNP${prev_snp} --allow-extra-chr --geno $geno_threshold --make-bed --out ${prefix}.common_sample${prev_sample}.SNP${current_snp}
    fi

    # Set sample threshold
    mind_threshold=$(echo "scale=2; 1-($current_sample/100)" | bc)

    # Apply sample filtering
    if [ $current_snp -eq $start ]; then
        plink --bfile ${prefix}.common_SNP${current_snp} --allow-extra-chr --mind $mind_threshold --make-bed --out ${prefix}.common_sample${current_sample}.SNP${current_snp}
    else
        plink --bfile ${prefix}.common_sample${prev_sample}.SNP${current_snp} --allow-extra-chr --mind $mind_threshold --make-bed --out ${prefix}.common_sample${current_sample}.SNP${current_snp}
    fi
    
    # Store current values for next iteration
    prev_snp=$current_snp
    prev_sample=$current_sample
    
    # Increment thresholds
    current_snp=$((current_snp + step))
    current_sample=$((current_sample + step))
done
