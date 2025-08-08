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
    path("*.IBD.kin0"), emit: kin0
    path("*.IBD.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    # Run KING kinship estimation via PLINK
    plink \\
        --bfile ${prefix} \\
        --allow-extra-chr \\
        --make-king-table \\
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
    path(kin0_file)       // KING kinship table (.kin0)
    path(imiss_file)      // PLINK .imiss file for call rates
    val(study_name)
    val(kinship_threshold)

    output:
    path("kinship_over_354_${study_name}.txt"), emit: twin_pairs
    path("kinship_btw_177_354_${study_name}.txt"), emit: sibling_pairs
    path("callrate_kinship_over_354_${study_name}.txt"), emit: twin_pairs_callrate
    path("*.IBD_outliers.txt"), emit: ibd_outliers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.data.split('/').last()}"
    
    """
    # KING kinship coefficient thresholds:
    #   >0.354 : duplicate/MZ twin
    # 0.177–0.354 : full siblings (1st degree)
    
    # Twins/duplicates
    awk '\$9 > 0.354 {print \$0}' ${kin0_file} > kinship_over_354_${study_name}.txt
    
    # Siblings
    awk '\$9 > 0.177 && \$9 <= 0.354 {print \$0}' ${kin0_file} > kinship_btw_177_354_${study_name}.txt
    
    # Prepare header for call rate table
    awk 'NR<2 {print \$0 }' kinship_over_354_${study_name}.txt > kinship_over_354_${study_name}_header.txt
    echo "     F_MISS_1     F_MISS_2" > extra_header.txt
    paste kinship_over_354_${study_name}_header.txt extra_header.txt > callrate_kinship_over_354_${study_name}_header.txt
    
    # Sort and join to add call rates
    sort -k1,1 kinship_over_354_${study_name}.txt > kinship_over_354_${study_name}_sorted.txt
    sort -k1,1 ${imiss_file} > ${imiss_file}_sorted.txt
    join -1 1 -2 1 -o auto kinship_over_354_${study_name}_sorted.txt ${imiss_file}_sorted.txt > callrate1_kinship_over_354_${study_name}.txt
    
    # Join again on second individual ID (IID2 in KING output is column 3)
    sort -k3,3 callrate1_kinship_over_354_${study_name}.txt > callrate1_kinship_over_354_${study_name}_sorted.txt
    join -1 3 -2 1 -o auto callrate1_kinship_over_354_${study_name}_sorted.txt ${imiss_file}_sorted.txt > callrate2_kinship_over_354_${study_name}.txt
    sort -k1 callrate2_kinship_over_354_${study_name}.txt > callrate3_kinship_over_354_${study_name}.txt
    
    # Final output with header
    cat callrate_kinship_over_354_${study_name}_header.txt callrate3_kinship_over_354_${study_name}.txt > callrate_kinship_over_354_${study_name}.txt
    
    # Clean up intermediates
    rm ./extra_header.txt
    rm ./callrate1_kinship_over_354_${study_name}.txt
    rm ./callrate2_kinship_over_354_${study_name}.txt
    rm ./callrate3_kinship_over_354_${study_name}.txt
    rm ./kinship_over_354_${study_name}_header.txt
    rm ./kinship_over_354_${study_name}_sorted.txt
    rm ./callrate1_kinship_over_354_${study_name}_sorted.txt
    rm ./${imiss_file}_sorted.txt
    
    # Mark outliers above threshold
    awk -v kin="${kinship_threshold}" '\$9 >= kin {print \$1, \$2}' ${kin0_file} > ${prefix}.IBD_outliers.txt

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
