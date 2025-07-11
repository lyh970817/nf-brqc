#!/usr/bin/env nextflow

/*
========================================================================================
    BioResource Genetic QC Pipeline
========================================================================================
    Nextflow implementation of the BioResource Genetic QC Pipeline
    Originally written as a bash script (fullpipe.sh)

    Authors:
    - Original bash script: BioResource team
    - Nextflow implementation: [Your Name]

    Usage:
    nextflow run main.nf --data /path/to/plink/files --name study_name [options]
*/

nextflow.enable.dsl = 2

// Import modules - Initial QC
include { MAF_FILTER } from './modules/local/initial_qc'
include { MISSINGNESS_CHECK } from './modules/local/initial_qc'
include { MISSINGNESS_HISTOGRAMS } from './modules/local/initial_qc'
include { ITERATIVE_MISSINGNESS } from './modules/local/initial_qc'
include { ITERATIVE_MISSINGNESS_TABLE } from './modules/local/initial_qc'
include { HWE_CHECK } from './modules/local/initial_qc'
include { HARDY_PLOTS } from './modules/local/initial_qc'
include { HWE_FILTER } from './modules/local/initial_qc'

// Import modules - LD pruning, sex checks, heterozygosity
include { LD_PRUNING } from './modules/local/ld_sex_het'
include { EXTRACT_PRUNED_SNPS } from './modules/local/ld_sex_het'
include { HIGH_LD_REGIONS_EXCLUDE } from './modules/local/ld_sex_het'
include { EXCLUDE_HIGH_LD_AUTOSOMAL } from './modules/local/ld_sex_het'
include { SEX_CHECK_SPLIT_X } from './modules/local/ld_sex_het'
include { SEX_CHECK } from './modules/local/ld_sex_het'
include { SEX_CHECK_PLOTS } from './modules/local/ld_sex_het'
include { HETEROZYGOSITY_CHECK } from './modules/local/ld_sex_het'
include { HETEROZYGOSITY_PLOTS } from './modules/local/ld_sex_het'

// Import modules - IBD analysis
include { IBD_CALCULATION } from './modules/local/ibd_analysis'
include { IBD_OUTLIER_DETECTION } from './modules/local/ibd_analysis'
include { INDIVIDUAL_IBD_ANALYSIS } from './modules/local/ibd_analysis'
include { IBD_HISTOGRAMS } from './modules/local/ibd_analysis'
include { INDIVIDUAL_IBD_HISTOGRAMS } from './modules/local/ibd_analysis'

// Import modules - Ancestry analysis
include { ANCESTRY_TARGET_SUBSET } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_TARGET_QC } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_REF_INTERSECT } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_REF_MERGE } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_LD_PRUNE } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_EXTRACT_PRUNED } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_PCA_WEIGHTS } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_REF_SCORE } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_TARGET_SCORE } from './modules/local/ancestry_plink_ops'
include { ANCESTRY_ALLELE_MATCHING } from './modules/local/ancestry_analysis'
include { ANCESTRY_LONG_LD_REGIONS } from './modules/local/ancestry_analysis'
include { ANCESTRY_PC_ANALYSIS } from './modules/local/ancestry_analysis'

// Function to help with error handling
def errorMessage() {
    log.info """
    =======================================================
    BioResource Genetic QC Pipeline: ERROR
    =======================================================
    The execution of the pipeline has failed.
    Please check the error message and fix the issue.

    Common issues:
    - Input data not found or incorrect format
    - Missing reference files
    - Software not installed or not in PATH

    For help, please contact the pipeline maintainer.
    =======================================================
    """
}

// Define the workflow
workflow {
    // Print header
    log.info """
    =======================================================
    BioResource Genetic QC Pipeline v1.0
    =======================================================

    Parameters:
    --data             : ${params.data}
    --name             : ${params.name}
    --outdir           : ${params.outdir}

    QC Thresholds:
    --maf              : ${params.maf}
    --SNP_CR           : ${params.SNP_CR}
    --Sample_CR        : ${params.Sample_CR}
    --geno             : ${params.geno}
    --mind             : ${params.mind}
    --hwe              : ${params.hwe}
    --ibd              : ${params.ibd}
    --ind_ibd          : ${params.ind_ibd}

    Software:
    --plink            : ${params.plink}
    --plink2           : ${params.plink2}

    =======================================================
    """

    // Check mandatory parameters
    if (params.data == null) {
        error "Input data not specified. Please provide --data parameter."
    }

    // Validate input files exist and create input channel
    def bed_file = Channel.fromPath("${params.data}.bed", checkIfExists: true)
    def bim_file = Channel.fromPath("${params.data}.bim", checkIfExists: true)
    def fam_file = Channel.fromPath("${params.data}.fam", checkIfExists: true)

    // Create input channel directly with params.data
    def input_ch = bed_file.combine(bim_file).combine(fam_file)

    // Run the pipeline
    // 1. Initial QC Steps
    MAF_FILTER(input_ch, params.maf)
    MISSINGNESS_CHECK(MAF_FILTER.out.plink_files)
    MISSINGNESS_HISTOGRAMS(MISSINGNESS_CHECK.out.lmiss, params.name)

    // 2. Iterative missingness
    ITERATIVE_MISSINGNESS(MISSINGNESS_CHECK.out.plink_files, 90, 99, 1)
    ITERATIVE_MISSINGNESS_TABLE(ITERATIVE_MISSINGNESS.out.logs, params.name)

    // 3. Select the 95% filtered files from iterative missingness
    def selected_files = ITERATIVE_MISSINGNESS.out.plink_files
        // break the (beds,bims,fams) lists into individual triples
        .flatMap { plink_files ->
            def beds = plink_files[0]
            def bims = plink_files[1]
            def fams = plink_files[2]
            beds.indices.collect { i -> [ beds[i], bims[i], fams[i] ] }
        }
        // keep only the triple whose BED file matches the requested sample/SNP
        .filter { bed, _bim, _fam ->
            bed.name.contains("sample${params.Sample_CR}.SNP${params.SNP_CR}")
        }

    // 4. Hardy-Weinberg analysis
    HWE_CHECK(selected_files)
    HARDY_PLOTS(HWE_CHECK.out.hwe, params.name)
    HWE_FILTER(selected_files, params.hwe)

    // 5. LD pruning
    LD_PRUNING(HWE_FILTER.out.plink_files, 1500, 150, 0.2)
    EXTRACT_PRUNED_SNPS(HWE_FILTER.out.plink_files, LD_PRUNING.out.prune_in)

    // 6. High LD regions and autosomal filtering
    HIGH_LD_REGIONS_EXCLUDE(EXTRACT_PRUNED_SNPS.out.plink_files.map { _bed, bim, _fam -> bim }, params.build)
    EXCLUDE_HIGH_LD_AUTOSOMAL(EXTRACT_PRUNED_SNPS.out.plink_files, HIGH_LD_REGIONS_EXCLUDE.out.exclude_list)

    // 7. Sex checks
    SEX_CHECK_SPLIT_X(EXTRACT_PRUNED_SNPS.out.plink_files, params.build)
    SEX_CHECK(SEX_CHECK_SPLIT_X.out.plink_files)
    SEX_CHECK_PLOTS(SEX_CHECK.out.sexcheck, params.name)

    // 8. Heterozygosity analysis
    HETEROZYGOSITY_CHECK(EXCLUDE_HIGH_LD_AUTOSOMAL.out.plink_files)
    HETEROZYGOSITY_PLOTS(HETEROZYGOSITY_CHECK.out.ibc, SEX_CHECK.out.sexcheck, params.name)

    // 9. IBD analysis
    IBD_CALCULATION(EXCLUDE_HIGH_LD_AUTOSOMAL.out.plink_files)
    IBD_OUTLIER_DETECTION(IBD_CALCULATION.out.genome, MISSINGNESS_CHECK.out.imiss, params.name, params.ibd)
    INDIVIDUAL_IBD_ANALYSIS(IBD_CALCULATION.out.genome, params.name, params.ind_ibd)
    IBD_HISTOGRAMS(IBD_CALCULATION.out.genome, params.name)
    INDIVIDUAL_IBD_HISTOGRAMS(INDIVIDUAL_IBD_ANALYSIS.out.individual_data, params.ind_ibd)

    // 10. Ancestry analysis (if reference data is available)
    if (params.ref_1kg_dir && file(params.ref_1kg_dir).exists()) {
        // Create reference channel for chromosomes 1-22
        def ref_ch = Channel.from(1..22).map { chr ->
            def ref_meta = [id: "ref.chr${chr}"]
            def ref_pgen = file("${params.ref_1kg_dir}/ref.chr${chr}.pgen")
            def ref_pvar = file("${params.ref_1kg_dir}/ref.chr${chr}.pvar")
            def ref_psam = file("${params.ref_1kg_dir}/ref.chr${chr}.psam")
            [ref_meta, ref_pgen, ref_pvar, ref_psam, chr]
        }

        // Run ancestry analysis workflow
        runAncestryAnalysis(HWE_FILTER.out.plink_files, ref_ch)
    }
}

// Ancestry analysis workflow
workflow runAncestryAnalysis {
    take:
    target_files
    ref_files_ch

    main:
    // Target QC for ancestry analysis
    target_files.view()
    def sample_size = 1000  // This should be calculated from the actual data
    ANCESTRY_TARGET_QC(target_files, params.geno, params.maf, params.hwe, sample_size)

    // Process reference files for each chromosome
    ANCESTRY_REF_INTERSECT(ref_files_ch, ANCESTRY_TARGET_QC.out.snplist.collect(), params.geno, params.maf, params.hwe)

    // Allele matching analysis
    ANCESTRY_ALLELE_MATCHING(
        ANCESTRY_REF_INTERSECT.out.plink_files.map { meta, bed, bim, fam -> [meta, bim] }.collect() 
            | map { files -> [files[0][0], files.collect { it[1] }] },
        target_files.map { _bed, bim, _fam -> bim }
    )

    // Merge reference files
    ANCESTRY_REF_MERGE(
        ANCESTRY_REF_INTERSECT.out.plink_files.collect() 
            | map { files -> [files[0][0], files.collect { it[1..3] }.flatten()] },
        ANCESTRY_ALLELE_MATCHING.out.match_snplist.map { _meta, snplist -> snplist },
        ANCESTRY_ALLELE_MATCHING.out.flip_snplist.map { _meta, snplist -> snplist }.ifEmpty(file('NO_FILE'))
    )

    // Long LD regions exclusion
    ANCESTRY_LONG_LD_REGIONS(ANCESTRY_REF_MERGE.out.plink_files.map { _meta, _pgen, pvar, _psam -> [_meta, pvar] })

    // LD pruning
    ANCESTRY_LD_PRUNE(ANCESTRY_REF_MERGE.out.plink_files, ANCESTRY_LONG_LD_REGIONS.out.exclude_list.map { _meta, exclude -> exclude })

    // Extract pruned SNPs
    ANCESTRY_EXTRACT_PRUNED(ANCESTRY_REF_MERGE.out.plink_files, ANCESTRY_LD_PRUNE.out.prune_files.map { _meta, prune_in, _prune_out -> prune_in })

    // PCA analysis
    ANCESTRY_PCA_WEIGHTS(ANCESTRY_EXTRACT_PRUNED.out.plink_files, params.n_pcs)
    ANCESTRY_REF_SCORE(ANCESTRY_EXTRACT_PRUNED.out.plink_files, ANCESTRY_PCA_WEIGHTS.out.allele_weights.map { _meta, weights -> weights }, params.n_pcs)

    // Score target samples
    def score_snplist = ANCESTRY_PCA_WEIGHTS.out.allele_weights.map { _meta, weights ->
        // Extract SNP list from weights file
        file(weights.toString().replace('.eigenvec.allele', '_snplist.txt'))
    }
    ANCESTRY_TARGET_SCORE(target_files, ANCESTRY_PCA_WEIGHTS.out.allele_weights.map { _meta, weights -> weights },
                         score_snplist, params.n_pcs, file('NO_FILE'))

    // PC analysis and population assignment
    if (params.pop_data && file(params.pop_data).exists()) {
        ANCESTRY_PC_ANALYSIS(ANCESTRY_REF_SCORE.out.scores, ANCESTRY_TARGET_SCORE.out.scores,
                            file(params.pop_data), file(params.ref_pop_scale),
                            params.n_pcs, params.prob_thresh, params.name)
    }
}


