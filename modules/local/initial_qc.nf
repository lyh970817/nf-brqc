/*
 * Initial QC processes
 */

process MAF_FILTER {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)
    val(maf_threshold)

    output:
    tuple path("*.bed"), path("*.bim"), path("*.fam"), emit: plink_files
    path("*.log"), emit: log
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
        --maf ${maf_threshold} \\
        --make-bed \\
        --out ${prefix}_maf${maf_threshold} \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process MISSINGNESS_CHECK {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("*_check.bed"), path("*_check.bim"), path("*_check.fam"), emit: plink_files
    path("*_check.lmiss"), emit: lmiss
    path("*_check.imiss"), emit: imiss
    path("*_check.hwe"), emit: hwe
    path("*_check.log"), emit: log
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
        --missing \\
        --hardy \\
        --make-bed \\
        --out ${prefix}_check \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process MISSINGNESS_HISTOGRAMS {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(lmiss)
    val(study_name)

    output:
    path("missingness_hist*${study_name}*"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    Rscript ${moduleDir}/../../bin/missingness_histograms.r ${lmiss} ${study_name} \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process ITERATIVE_MISSINGNESS {
    tag "${params.data}"
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)
    val(start_threshold)
    val(end_threshold)
    val(step_size)

    output:
    tuple path("*.common_sample*.SNP*.bed"), path("*.common_sample*.SNP*.bim"), path("*.common_sample*.SNP*.fam"), emit: plink_files
    path("*.common_*.log"), emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    def memory = task.memory ? "--memory ${task.memory.toMega()}" : ""
    
    """
    # Run iterative missingness script
    bash ${moduleDir}/../../bin/Iterative_Missingness.sh ${start_threshold} ${end_threshold} ${step_size} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process ITERATIVE_MISSINGNESS_TABLE {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(log_files)
    val(study_name)

    output:
    path("Iterative_missingness_table_${study_name}.txt"), emit: table
    path("itmiss_table_${study_name}*"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    # Process log files to create iterative missingness table
    # Get first log value
    grep removed *.common_SNP*.log | tr -d ' variants removed due to missing genotype data (--geno).' > geno90_int.txt
    
    # Get every other geno/mind log value
    sed -n '/removed/p' *.common_sample*.SNP*.log | tr -d ' <a-z>(--).' > otherlogs_int.txt
    
    # Concatenate these text files
    cat geno90_int.txt otherlogs_int.txt > logvalues_int.txt
    
    # Get alternate columns for geno/mind values
    sed 'N;s/\\n/ /' logvalues_int.txt > split_cols_int.txt
    
    # Select relevant column (geno/mind) and add 0s to every other row
    awk 'BEGIN { FS=" " }; {print \$1}' split_cols_int.txt | sed '2,\$s/^/0\\n/' > geno0_int.txt
    echo "0" >> geno0_int.txt
    awk 'BEGIN { FS=" " }; {print \$2}' split_cols_int.txt | sed "s/^/0\\n/" > mind0_int.txt
    
    # Rejoin files
    paste geno0_int.txt mind0_int.txt > genomind_int.txt
    
    # Add sample/snp threshold column
    printf '%s\\n' "90" "90/90" "90/91" "91/91" "91/92" "92/92" "92/93" "93/93" "93/94" "94/94" "94/95" "95/95" "95/96" "96/96" "96/97" "97/97" "97/98" "98/98" "98/99" "99/99" > samplesnp_int.txt
    paste samplesnp_int.txt genomind_int.txt > 3col_genomind_int.txt
    
    # Add column names
    awk 'BEGIN { OFS="\\t"; print "Sample/SNP\\t" "Variants\\t" "Participants"}; { print \$0 }' 3col_genomind_int.txt > Iterative_missingness_table_${study_name}.txt
    
    # Remove intermediate files
    rm *_int.txt
    
    # Make nice table
    Rscript ${moduleDir}/../../bin/itmiss_table.r ${study_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process HWE_CHECK {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    tuple path("*_hardycheck.bed"), path("*_hardycheck.bim"), path("*_hardycheck.fam"), emit: plink_files
    path("*_hardycheck.hwe"), emit: hwe
    path("*_hardycheck.log"), emit: log
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
        --hardy \\
        --make-bed \\
        --out ${prefix}_hardycheck \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}

process HARDY_PLOTS {
    tag "${params.data}"
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(hwe_file)
    val(study_name)

    output:
    path("hardy_*${study_name}*"), emit: plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    Rscript ${moduleDir}/../../bin/hardy_plots.r ${hwe_file} ${study_name} \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //; s/ .*//')
    END_VERSIONS
    """
}

process HWE_FILTER {
    tag "${params.data}"
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    tuple path(bed), path(bim), path(fam)
    val(hwe_threshold)

    output:
    tuple path("*.hwe*.bed"), path("*.hwe*.bim"), path("*.hwe*.fam"), emit: plink_files
    path("*.hwe*.log"), emit: log
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
        --hwe ${hwe_threshold} \\
        --make-bed \\
        --out ${prefix}.hwe${hwe_threshold} \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK //; s/64-bit.*//')
    END_VERSIONS
    """
}
