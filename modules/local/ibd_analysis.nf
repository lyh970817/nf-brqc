/*
 * IBD (Identity by Descent) analysis processes
 */

process IBD_CALCULATION {
    tag "${params.data}"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("*.IBD.bed"), path("*.IBD.bim"), path("*.IBD.fam"), emit: plink_files
    path("*.IBD.genome"), emit: genome
    path("*.IBD.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --bfile ${prefix} \\
        --allow-extra-chr \\
        --genome \\
        --make-bed \\
        --out ${prefix}.IBD \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process IBD_OUTLIER_DETECTION {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(genome_file)
    path(imiss_file)
    val(study_name)
    val(ibd_threshold)

    output:
    path("pihat_over_9_${study_name}.txt"), emit: twin_pairs
    path("pihat_btw_4_6_${study_name}.txt"), emit: sibling_pairs
    path("callrate_pihat_over_9_${study_name}.txt"), emit: twin_pairs_callrate
    path("*.IBD_outliers.txt"), emit: ibd_outliers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.data.split('/').last()}"
    
    """
    # Check >0.9 Pi_Hat individuals: duplicates or twins
    awk '\$10 > 0.9 {print \$0}' ${genome_file} > pihat_over_9_${study_name}.txt
    
    # Check 0.4 < \$10 < 0.6 Pi_Hat individuals: family
    awk '\$10 > 0.4 && \$10 < 0.6 {print \$0}' ${genome_file} > pihat_btw_4_6_${study_name}.txt
    
    # Add call rate info to pihat outputs
    # Get header
    awk 'NR<2 {print \$0 }' pihat_over_9_${study_name}.txt > pihat_over_9_${study_name}_header.txt
    
    # Add Call Rate (F_MISS) to header
    echo "     F_MISS_1     F_MISS_2" > extra_header.txt
    paste pihat_over_9_${study_name}_header.txt extra_header.txt > callrate_pihat_over_9_${study_name}_header.txt
    
    # Sort and join pihat and call rate files
    join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.6 <(sort -k1 pihat_over_9_${study_name}.txt) <(sort -k1 ${imiss_file}) > callrate1_pihat_over_9_${study_name}.txt
    join -1 3 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.6 <(sort -k3 callrate1_pihat_over_9_${study_name}.txt) <(sort -k1 ${imiss_file}) > callrate2_pihat_over_9_${study_name}.txt
    sort -k1 callrate2_pihat_over_9_${study_name}.txt > callrate3_pihat_over_9_${study_name}.txt
    
    # Attach new header
    cat callrate_pihat_over_9_${study_name}_header.txt callrate3_pihat_over_9_${study_name}.txt > callrate_pihat_over_9_${study_name}.txt
    
    # Remove intermediate callrate files
    rm ./extra_header.txt
    rm ./callrate1_pihat_over_9_${study_name}.txt
    rm ./callrate2_pihat_over_9_${study_name}.txt
    rm ./callrate3_pihat_over_9_${study_name}.txt
    rm ./pihat_over_9_${study_name}_header.txt
    
    # Remove one sample from each pair with pi-hat (% IBD) above threshold
    awk -v ibd="${ibd_threshold}" '\$10 >= ibd {print \$1, \$2}' ${genome_file} > ${prefix}.IBD_outliers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version 2>&1 | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}

process INDIVIDUAL_IBD_ANALYSIS {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(genome_file)
    val(study_name)
    val(sd_threshold)

    output:
    path("*.IBD_INDIV_outliers.txt"), emit: individual_outliers
    path("*.IBD_INDIV.txt"), emit: individual_data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.data.split('/').last()}"
    
    """
    R --file=${moduleDir}/../../bin/IndividualIBD.r --args ${prefix} ${sd_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process IBD_HISTOGRAMS {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(genome_file)
    val(study_name)

    output:
    path("*_${study_name}.pdf"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.data.split('/').last()}"
    
    """
    # Plot all IBD Histograms
    Rscript ${moduleDir}/../../bin/IBD_Hist.r ${prefix} ${study_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process INDIVIDUAL_IBD_HISTOGRAMS {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(individual_data_file)
    val(sd_threshold)

    output:
    path("IndvIBD_hist_*.pdf"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.data.split('/').last()}"
    
    """
    # Plot Individual IBD outlier Histograms
    Rscript ${moduleDir}/../../bin/IndvIBD_Hist.r ${prefix} ${sd_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}
