/*
 * Extract SNP list from PCA weights file
 */

process EXTRACT_SNPS_FROM_WEIGHTS {
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(weights_file)

    output:
    path("snp_list.txt"), emit: snp_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract SNP names from weights file (assuming they're in the first column after header)
    awk 'NR>1 {print \$2}' ${weights_file} > snp_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/^.*awk //; s/,.*//')
    END_VERSIONS
    """
}