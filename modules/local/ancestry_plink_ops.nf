/*
 * Ancestry analysis - PLINK operations
 * Separated from R analysis to follow Nextflow best practices
 */

process ANCESTRY_TARGET_SUBSET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(target_keep)
    path(target_fam)

    output:
    tuple val(meta), path("target_subset.bed"), path("target_subset.bim"), path("target_subset.fam"), emit: plink_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keep_arg = target_keep ? "--keep ${target_keep}" : ""
    def fam_arg = target_fam ? "--fam ${target_fam}" : ""
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    plink \\
        --bfile ${prefix} \\
        ${keep_arg} \\
        ${fam_arg} \\
        --allow-extra-chr \\
        --make-bed \\
        --out target_subset \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_TARGET_QC {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple path(bed), path(bim), path(fam)
    val(geno_threshold)
    val(maf_threshold)
    val(hwe_threshold)
    val(sample_size)

    output:
    path("target.QC.snplist"), emit: snplist
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    // Apply different QC based on sample size
    def qc_args = sample_size > 100 ? 
        "--geno ${geno_threshold} --maf ${maf_threshold} --hwe ${hwe_threshold}" :
        "--geno ${geno_threshold}"
    
    """
    plink \\
        --bfile ${prefix} \\
        ${qc_args} \\
        --write-snplist \\
        --allow-extra-chr \\
        --out target.QC \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_REF_INTERSECT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(ref_bed), path(ref_bim), path(ref_fam)
    path(target_snplist)
    val(geno_threshold)
    val(maf_threshold)
    val(hwe_threshold)
    val(chr)

    output:
    tuple val(meta), path("ref_intersect_chr${chr}.bed"), path("ref_intersect_chr${chr}.bim"), path("ref_intersect_chr${chr}.fam"), emit: plink_files
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
        --geno ${geno_threshold} \\
        --maf ${maf_threshold} \\
        --hwe ${hwe_threshold} \\
        --make-bed \\
        --extract ${target_snplist} \\
        --allow-extra-chr \\
        --out ref_intersect_chr${chr} \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_REF_MERGE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(ref_files)
    path(allele_match_snplist)
    path(flip_snplist)

    output:
    tuple val(meta), path("ref_merge.bed"), path("ref_merge.bim"), path("ref_merge.fam"), emit: plink_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    def flip_arg = flip_snplist.name != 'NO_FILE' ? "--flip ${flip_snplist}" : ""
    
    """
    # Create merge list
    ls ref_intersect_chr*.bed | sed 's/.bed//' > ref_mergelist.txt

    plink \\
        --merge-list ref_mergelist.txt \\
        --extract ${allele_match_snplist} \\
        ${flip_arg} \\
        --make-bed \\
        --allow-extra-chr \\
        --out ref_merge \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_LD_PRUNE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1' :
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(long_ld_exclude)

    output:
    tuple val(meta), path("ref_merge.prune.in"), path("ref_merge.prune.out"), emit: prune_files
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
        --exclude ${long_ld_exclude} \\
        --indep-pairwise 1000 5 0.2 \\
        --allow-extra-chr \\
        --out ref_merge \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_EXTRACT_PRUNED {
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
    tuple val(meta), path("ref_merge_pruned.bed"), path("ref_merge_pruned.bim"), path("ref_merge_pruned.fam"), emit: plink_files
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
        --extract ${prune_in} \\
        --make-bed \\
        --allow-extra-chr \\
        --out ref_merge_pruned \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_PCA_WEIGHTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_3' :
        'quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_3' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val(n_pcs)

    output:
    tuple val(meta), path("ref_merge_pruned.eigenvec.allele"), emit: allele_weights
    path("ref_merge_pruned.eigenval"), emit: eigenval
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""

    """
    plink2 \\
        --bfile ${prefix} \\
        --pca ${n_pcs} allele-wts \\
        --allow-extra-chr \\
        --out ref_merge_pruned \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK //; s/ 64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_REF_SCORE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_3' :
        'quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_3' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(allele_weights)
    val(n_pcs)

    output:
    tuple val(meta), path("ref_merge_pruned_score.sscore"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    def score_cols = "6-${n_pcs + 5}"

    """
    plink2 \\
        --bfile ${prefix} \\
        --score ${allele_weights} header-read 2 5 \\
        --score-col-nums ${score_cols} \\
        --allow-extra-chr \\
        --out ref_merge_pruned_score \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK //; s/ 64-bit.*//')
    END_VERSIONS
    """
}

process ANCESTRY_TARGET_SCORE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_3' :
        'quay.io/biocontainers/plink2:2.00a3.7--h4ac6f70_3' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path(allele_weights)
    path(score_snplist)
    val(n_pcs)
    path(target_fam)

    output:
    tuple val(meta), path("profiles.sscore"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    def score_cols = "6-${n_pcs + 5}"
    def fam_arg = target_fam.name != 'NO_FILE' ? "--fam ${target_fam}" : ""

    """
    plink2 \\
        --bfile ${prefix} \\
        ${fam_arg} \\
        --extract ${score_snplist} \\
        --score ${allele_weights} header-read 2 5 no-mean-imputation \\
        --score-col-nums ${score_cols} \\
        --allow-extra-chr \\
        --out profiles \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK //; s/ 64-bit.*//')
    END_VERSIONS
    """
}
