/*
 * LD pruning, sex checks, and heterozygosity analysis
 */

process LD_PRUNING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val(window_size)
    val(step_size)
    val(r2_threshold)

    output:
    tuple val(meta), path("*.LD_Pre.bed"), path("*.LD_Pre.bim"), path("*.LD_Pre.fam"), emit: ld_pre_files
    path("*.LD_Pre.prune.in"), emit: prune_in
    path("*.LD_Pre.prune.out"), emit: prune_out
    path("*.LD_Pre.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --bfile ${prefix} \\
        --allow-extra-chr \\
        --indep-pairwise ${window_size} ${step_size} ${r2_threshold} \\
        --make-bed \\
        --out ${prefix}.LD_Pre \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process EXTRACT_PRUNED_SNPS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(prune_in)

    output:
    tuple val(meta), path("*.LD_Pruned.bed"), path("*.LD_Pruned.bim"), path("*.LD_Pruned.fam"), emit: plink_files
    path("*.LD_Pruned.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --bfile ${prefix} \\
        --allow-extra-chr \\
        --extract ${prune_in} \\
        --make-bed \\
        --out ${prefix}.LD_Pruned \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process HIGH_LD_REGIONS_EXCLUDE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bim)
    val(build)

    output:
    tuple val(meta), path("highLD_and_autosomal_excludes"), emit: exclude_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def awk_script = build == "b38" ? "${moduleDir}/../../bin/highLDregions4bim_b38.awk" : "${moduleDir}/../../bin/highLDregions4bim_b37.awk"
    
    """
    # Make highLDregion and Autosomalexcludes SNP file
    awk -f ${awk_script} ${bim} > highLDexcludes
    awk '(\$1 < 1) || (\$1 > 22) {print \$2}' ${bim} > autosomeexcludes
    cat highLDexcludes autosomeexcludes > highLD_and_autosomal_excludes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version 2>&1 | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}

process EXCLUDE_HIGH_LD_AUTOSOMAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(exclude_list)

    output:
    tuple val(meta), path("*_autosomalchr.bed"), path("*_autosomalchr.bim"), path("*_autosomalchr.fam"), emit: plink_files
    path("*_autosomalchr.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --bfile ${prefix} \\
        --allow-extra-chr \\
        --exclude ${exclude_list} \\
        --make-bed \\
        --out ${prefix}_autosomalchr \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process SEX_CHECK_SPLIT_X {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val(build)

    output:
    tuple val(meta), path("*_split_x.bed"), path("*_split_x.bim"), path("*_split_x.fam"), emit: plink_files
    path("*_split_x.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --allow-extra-chr \\
        --bfile ${prefix} \\
        --split-x 'no-fail' ${build} \\
        --make-bed \\
        --out ${prefix}_split_x \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process SEX_CHECK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*_sexcheck.sexcheck"), emit: sexcheck
    path("*_sexcheck.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --allow-extra-chr \\
        --bfile ${prefix} \\
        --check-sex ycount \\
        --set-hh-missing \\
        --out ${prefix}_sexcheck \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process SEX_CHECK_PLOTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.0' :
        'quay.io/biocontainers/r-base:4.3.0' }"

    input:
    tuple val(meta), path(sexcheck_file)
    val(study_name)

    output:
    tuple val(meta), path("sexcheck_hist_${study_name}*"), emit: plots
    tuple val(meta), path("sexmismatch_${study_name}.txt"), emit: mismatch_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    Rscript ${moduleDir}/../../bin/sex_check_hist.r ${sexcheck_file} ${study_name} \$PWD

    # Extract rows with sex discrepancies and export to a table
    grep PROBLEM ${sexcheck_file} > sexmismatch_${study_name}.txt || touch sexmismatch_${study_name}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process HETEROZYGOSITY_CHECK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*_Hetcheck.bed"), path("*_Hetcheck.bim"), path("*_Hetcheck.fam"), emit: plink_files
    path("*_Hetcheck.ibc"), emit: ibc
    path("*_Hetcheck.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --allow-extra-chr \\
        --bfile ${prefix} \\
        --ibc \\
        --make-bed \\
        --out ${prefix}_Hetcheck \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process HETEROZYGOSITY_PLOTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.0' :
        'quay.io/biocontainers/r-base:4.3.0' }"

    input:
    tuple val(meta), path(ibc_file)
    path(sexcheck_file)
    val(study_name)

    output:
    tuple val(meta), path("Het_check_hist_${study_name}.pdf"), emit: het_hist
    tuple val(meta), path("HetFhat2_SexF_${study_name}.pdf"), emit: het_sex_plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    Rscript ${moduleDir}/../../bin/het_graphs.r ${ibc_file} ${sexcheck_file} ${study_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}
