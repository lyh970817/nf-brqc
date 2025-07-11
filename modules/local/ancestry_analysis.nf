/*
 * Ancestry analysis - R-based analysis
 * Pure R analysis separated from PLINK operations
 */

process ANCESTRY_ALLELE_MATCHING {
    label 'process_medium'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(ref_pvar_files)
    path(target_bim_file)

    output:
    path("ref_allele_match.snplist"), emit: match_snplist
    path("ref_flip.snplist"), emit: flip_snplist, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    
    # Read reference pvar files
    ref_bim <- NULL
    ref_files <- list.files(pattern = ".*\\\\.pvar")
    for(file in ref_files) {
        ref_data <- fread(file, skip = 1)  # Skip header line with ##
        ref_bim <- rbind(ref_bim, ref_data)
    }
    
    # PVAR format: #CHROM POS ID REF ALT
    ref_bim <- ref_bim[,c('V1','V3','V2','V4','V5')]
    names(ref_bim) <- c('CHR','SNP','BP','A1','A2')
    
    # Read target bim files
    targ_bim <- NULL
    targ_files <- list.files(pattern = ".*\\\\.bim")
    targ_files <- targ_files[!grepl("ref_intersect", targ_files)]
    for(file in targ_files) {
        targ_bim <- rbind(targ_bim, fread(file))
    }
    
    targ_bim <- targ_bim[,c('V1','V2','V4','V5','V6')]
    names(targ_bim) <- c('CHR','SNP','BP','A1','A2')
    
    # Create IUPAC codes in target data
    targ_bim\$IUPAC[targ_bim\$A1 == 'A' & targ_bim\$A2 =='T' | targ_bim\$A1 == 'T' & targ_bim\$A2 =='A'] <- 'W'
    targ_bim\$IUPAC[targ_bim\$A1 == 'C' & targ_bim\$A2 =='G' | targ_bim\$A1 == 'G' & targ_bim\$A2 =='C'] <- 'S'
    targ_bim\$IUPAC[targ_bim\$A1 == 'A' & targ_bim\$A2 =='G' | targ_bim\$A1 == 'G' & targ_bim\$A2 =='A'] <- 'R'
    targ_bim\$IUPAC[targ_bim\$A1 == 'C' & targ_bim\$A2 =='T' | targ_bim\$A1 == 'T' & targ_bim\$A2 =='C'] <- 'Y'
    targ_bim\$IUPAC[targ_bim\$A1 == 'G' & targ_bim\$A2 =='T' | targ_bim\$A1 == 'T' & targ_bim\$A2 =='G'] <- 'K'
    targ_bim\$IUPAC[targ_bim\$A1 == 'A' & targ_bim\$A2 =='C' | targ_bim\$A1 == 'C' & targ_bim\$A2 =='A'] <- 'M'
    targ_bim\$SNP_IUPAC <- paste0(targ_bim\$SNP,':',targ_bim\$IUPAC)
    
    # Create IUPAC codes in ref data
    ref_bim\$IUPAC[ref_bim\$A1 == 'A' & ref_bim\$A2 =='T' | ref_bim\$A1 == 'T' & ref_bim\$A2 =='A'] <- 'W'
    ref_bim\$IUPAC[ref_bim\$A1 == 'C' & ref_bim\$A2 =='G' | ref_bim\$A1 == 'G' & ref_bim\$A2 =='C'] <- 'S'
    ref_bim\$IUPAC[ref_bim\$A1 == 'A' & ref_bim\$A2 =='G' | ref_bim\$A1 == 'G' & ref_bim\$A2 =='A'] <- 'R'
    ref_bim\$IUPAC[ref_bim\$A1 == 'C' & ref_bim\$A2 =='T' | ref_bim\$A1 == 'T' & ref_bim\$A2 =='C'] <- 'Y'
    ref_bim\$IUPAC[ref_bim\$A1 == 'G' & ref_bim\$A2 =='T' | ref_bim\$A1 == 'T' & ref_bim\$A2 =='G'] <- 'K'
    ref_bim\$IUPAC[ref_bim\$A1 == 'A' & ref_bim\$A2 =='C' | ref_bim\$A1 == 'C' & ref_bim\$A2 =='A'] <- 'M'
    ref_bim\$SNP_IUPAC <- paste0(ref_bim\$SNP,':',ref_bim\$IUPAC)
    
    # Merge target and reference based on SNP id
    ref_target <- merge(ref_bim, targ_bim, by='SNP')
    
    # Identify SNPs for which alleles need to be flipped
    flip_tmp <- ref_target[(ref_target\$IUPAC.x == 'R' & ref_target\$IUPAC.y == 'Y' |
                            ref_target\$IUPAC.x == 'Y' & ref_target\$IUPAC.y == 'R' |
                            ref_target\$IUPAC.x == 'K' & ref_target\$IUPAC.y == 'M' |
                            ref_target\$IUPAC.x == 'M' & ref_target\$IUPAC.y == 'K'),]
    
    # Identify SNPs which match the reference alleles
    incl <- ref_target[ ref_target\$IUPAC.x == 'R' & ref_target\$IUPAC.y == 'R' |
                        ref_target\$IUPAC.x == 'Y' & ref_target\$IUPAC.y == 'Y' |
                        ref_target\$IUPAC.x == 'K' & ref_target\$IUPAC.y == 'K' |
                        ref_target\$IUPAC.x == 'M' & ref_target\$IUPAC.y == 'M' ]
    
    # If a SNP that needs to be flipped has a duplicate that is on the correct strand, remove it.
    flip <- flip_tmp[!(flip_tmp\$SNP %in% incl\$SNP)]
    
    # Combine SNPs that match and those that need to be flipped.
    incl <- rbind(incl,flip)
    
    if(dim(flip)[1] > 0){
        write.table(flip\$SNP, 'ref_flip.snplist', col.names=F, row.names=F, quote=F)
    }
    
    write.table(incl\$SNP, 'ref_allele_match.snplist', col.names=F, row.names=F, quote=F)
    
    # Write versions
    writeLines(c('"${task.process}":',
                 paste0('    r-base: "', R.version.string, '"')), 
               "versions.yml")
    """
}

process ANCESTRY_LONG_LD_REGIONS {
    label 'process_low'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(bim)

    output:
    path("long_ld.exclude"), emit: exclude_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    
    # Read in the bim file
    ref_bim <- data.frame(fread("${bim}"))
    
    # Create file removing long LD regions
    long_ld_exclude <- ref_bim\$V2[ (ref_bim\$V1 == 1 & ref_bim\$V4 >= 48e6 & ref_bim\$V4 <= 52e6) |
                                    (ref_bim\$V1 == 2 & ref_bim\$V4 >= 86e6 & ref_bim\$V4 <= 100.5e6) |
                                    (ref_bim\$V1 == 2 & ref_bim\$V4 >= 134.5e6 & ref_bim\$V4 <= 138e6) |
                                    (ref_bim\$V1 == 2 & ref_bim\$V4 >= 183e6 & ref_bim\$V4 <= 190e6) |
                                    (ref_bim\$V1 == 3 & ref_bim\$V4 >= 47.5e6 & ref_bim\$V4 <= 50e6) |
                                    (ref_bim\$V1 == 3 & ref_bim\$V4 >= 83.5e6 & ref_bim\$V4 <= 87e6) |
                                    (ref_bim\$V1 == 3 & ref_bim\$V4 >= 89e6 & ref_bim\$V4 <= 97.5e6) |
                                    (ref_bim\$V1 == 5 & ref_bim\$V4 >= 44.5e6 & ref_bim\$V4 <= 50.5e6) |
                                    (ref_bim\$V1 == 5 & ref_bim\$V4 >= 98e6 & ref_bim\$V4 <= 100.5e6) |
                                    (ref_bim\$V1 == 5 & ref_bim\$V4 >= 129e6 & ref_bim\$V4 <= 132e6) |
                                    (ref_bim\$V1 == 5 & ref_bim\$V4 >= 135.5e6 & ref_bim\$V4 <= 138.5e6) |
                                    (ref_bim\$V1 == 6 & ref_bim\$V4 >= 25.5e6 & ref_bim\$V4 <= 33.5e6) |
                                    (ref_bim\$V1 == 6 & ref_bim\$V4 >= 57e6 & ref_bim\$V4 <= 64e6) |
                                    (ref_bim\$V1 == 6 & ref_bim\$V4 >= 140e6 & ref_bim\$V4 <= 142.5e6) |
                                    (ref_bim\$V1 == 7 & ref_bim\$V4 >= 55e6 & ref_bim\$V4 <= 66e6) |
                                    (ref_bim\$V1 == 8 & ref_bim\$V4 >= 8e6 & ref_bim\$V4 <= 12e6) |
                                    (ref_bim\$V1 == 8 & ref_bim\$V4 >= 43e6 & ref_bim\$V4 <= 50e6) |
                                    (ref_bim\$V1 == 8 & ref_bim\$V4 >= 112e6 & ref_bim\$V4 <= 115e6) |
                                    (ref_bim\$V1 == 10 & ref_bim\$V4 >= 37e6 & ref_bim\$V4 <= 43e6) |
                                    (ref_bim\$V1 == 11 & ref_bim\$V4 >= 46e6 & ref_bim\$V4 <= 57e6) |
                                    (ref_bim\$V1 == 11 & ref_bim\$V4 >= 87.5e6 & ref_bim\$V4 <= 90.5e6) |
                                    (ref_bim\$V1 == 12 & ref_bim\$V4 >= 33e6 & ref_bim\$V4 <= 40e6) |
                                    (ref_bim\$V1 == 12 & ref_bim\$V4 >= 109.5e6 & ref_bim\$V4 <= 112e6) |
                                    (ref_bim\$V1 == 20 & ref_bim\$V4 >= 32e6 & ref_bim\$V4 <= 34.5e6)]
    
    write.table(long_ld_exclude, 'long_ld.exclude', col.names=F, row.names=F, quote=F)
    
    # Write versions
    writeLines(c('"${task.process}":',
                 paste0('    r-base: "', R.version.string, '"')), 
               "versions.yml")
    """
}

process ANCESTRY_PC_ANALYSIS {
    label 'process_high'

    conda "${moduleDir}/../../conda/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'bioresource-qc.sif' :
        'bioresource-qc:latest' }"

    input:
    path(ref_scores)
    path(target_scores)
    path(pop_data)
    path(ref_pop_scale_dir)
    val(n_pcs)
    val(prob_thresh)
    val(output_name)

    output:
    path("${output_name}.eigenvec"), emit: ref_eigenvec
    path("${output_name}.scale"), emit: scale_file
    path("${output_name}.*.scale"), emit: pop_scale_files, optional: true
    path("${output_name}.pop_model.rds"), emit: model, optional: true
    path("${output_name}.model_pred"), emit: model_pred, optional: true
    path("${output_name}.*.keep"), emit: keep_files
    path("${output_name}.*.eigenvec"), emit: pop_eigenvec
    path("${output_name}.PCs_plot_*.png"), emit: plots
    path("${output_name}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(caret)
    library(pROC)
    library(verification)
    library(ggplot2)
    library(cowplot)
    
    setDTthreads(${task.cpus})
    
    # Start logging
    sink(file = "${output_name}.log", append = F)
    cat('PC Analysis started at', as.character(Sys.time()), '\\n')
    sink()
    
    # Read in reference PC scores
    PCs_ref <- data.frame(fread("${ref_scores}"))
    # Extract IID (column 1) and PC columns (starting from column 7)
    pc_start_col <- 7
    pc_end_col <- pc_start_col + ${n_pcs} - 1
    PCs_ref <- PCs_ref[,c(1, pc_start_col:pc_end_col)]
    names(PCs_ref) <- c('IID', paste0('PC',1:${n_pcs}))
    # Add FID column (duplicate of IID for compatibility)
    PCs_ref <- data.frame(FID=PCs_ref[,1], PCs_ref)
    
    fwrite(PCs_ref, "${output_name}.eigenvec", sep=' ')
    
    # Scale across all individuals
    PCs_ref_centre_scale <- data.frame(PC=names(PCs_ref[-1:-2]),
                                      Mean=sapply(PCs_ref[,-1:-2], function(x) mean(x)),
                                      SD=sapply(PCs_ref[,-1:-2], function(x) sd(x)),
                                      row.names=seq(1:${n_pcs}))
    
    fwrite(PCs_ref_centre_scale, "${output_name}.scale", sep=' ')
    
    # Process population-specific scaling if ref_pop_scale_dir is provided
    if(dir.exists("${ref_pop_scale_dir}")) {
        # Get all .keep files from the directory
        keep_files <- list.files("${ref_pop_scale_dir}", pattern = "\\\\.keep\$", full.names = TRUE)
        
        for(k in 1:length(keep_files)){
            # Extract population name from filename (remove .keep extension)
            pop <- gsub("\\\\.keep\$", "", basename(keep_files[k]))
            keep <- fread(keep_files[k], header=F)
            PCs_ref_keep <- PCs_ref[(PCs_ref\$FID %in% keep\$V1),]
            
            PCs_ref_centre_scale_pop <- data.frame(PC=names(PCs_ref_keep[-1:-2]),
                                                  Mean=sapply(PCs_ref_keep[,-1:-2], function(x) mean(x)),
                                                  SD=sapply(PCs_ref_keep[,-1:-2], function(x) sd(x)),
                                                  row.names=seq(1:${n_pcs}))
            
            fwrite(PCs_ref_centre_scale_pop, paste0("${output_name}.", pop, ".scale"), sep=' ')
        }
        
        # Build prediction model
        PCs_ref_centre_scale <- fread("${output_name}.scale")
        
        # Scale the reference PCs
        PCs_ref_scaled <- PCs_ref
        for(i in 1:dim(PCs_ref_centre_scale)[1]){
            PCs_ref_scaled[[paste0('PC',i)]] <- PCs_ref[[paste0('PC',i)]]-PCs_ref_centre_scale\$Mean[PCs_ref_centre_scale\$PC == paste0('PC',i)]
            PCs_ref_scaled[[paste0('PC',i)]] <- PCs_ref_scaled[[paste0('PC',i)]]/PCs_ref_centre_scale\$SD[PCs_ref_centre_scale\$PC == paste0('PC',i)]
        }
        
        # Label individuals with ref_pop groups
        pop <- NULL
        for(i in 1:length(keep_files)){
            keep <- fread(keep_files[i], header=F)
            # Extract population name from filename (remove .keep extension)
            pop_name <- gsub("\\\\.keep\$", "", basename(keep_files[i]))
            keep\$pop <- pop_name
            pop <- rbind(pop,keep)
        }
        names(pop) <- c('IID','pop')
        pop\$FID <- pop\$IID

        PCs_ref_scaled_pop <- merge(PCs_ref_scaled,pop, by=c('FID','IID'))
        
        # Build model
        model <- train(y=as.factor(PCs_ref_scaled_pop\$pop), 
                      x=PCs_ref_scaled_pop[grepl('PC',names(PCs_ref_scaled_pop))], 
                      method='glmnet', 
                      metric='logLoss', 
                      trControl=trainControl(method="cv", number=5, classProbs= TRUE, 
                                           savePredictions = 'final', 
                                           summaryFunction = multiClassSummary))
        
        saveRDS(model\$finalModel, "${output_name}.pop_model.rds")
        
        # Generate predictions for target samples
        sink(file = "${output_name}.log", append = T)
        cat('Generating population predictions...\\n')
        sink()
        
        # Read in target PC scores
        target_PCs <- data.frame(fread("${target_scores}"))
        # Target scores from plink2 --score have format: FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM PC1 PC2 ...
        target_PCs <- target_PCs[,c(1:2,5:dim(target_PCs)[2])]
        names(target_PCs) <- c('FID','IID',paste0('PC',1:${n_pcs}))
        
        # Scale target PCs using reference scaling
        target_PCs_scaled <- target_PCs
        for(i in 1:dim(PCs_ref_centre_scale)[1]){
            target_PCs_scaled[[paste0('PC',i)]] <- target_PCs[[paste0('PC',i)]]-PCs_ref_centre_scale\$Mean[PCs_ref_centre_scale\$PC == paste0('PC',i)]
            target_PCs_scaled[[paste0('PC',i)]] <- target_PCs_scaled[[paste0('PC',i)]]/PCs_ref_centre_scale\$SD[PCs_ref_centre_scale\$PC == paste0('PC',i)]
        }
        
        # Generate predictions
        pop_model_pred <- predict(object = model\$finalModel, newx = data.matrix(target_PCs_scaled[grepl('PC',names(target_PCs_scaled))]), type = "response", s=model\$finalModel\$lambdaOpt)
        pop_model_pred <- as.data.frame.table(pop_model_pred)
        pop_model_pred <- data.table(FID=target_PCs_scaled\$FID,
                                    IID=target_PCs_scaled\$IID,
                                    pop=as.character(pop_model_pred\$Var2),
                                    prob=round(pop_model_pred\$Freq,3))
        
        pop_model_pred <- dcast.data.table(pop_model_pred, formula=FID + IID~pop, value.var = "prob")
        fwrite(pop_model_pred, "${output_name}.model_pred", sep='\\t')
        
        # Create keep files based on predictions
        if(${prob_thresh} > 0){
            pop_model_pred\$max_prob <- apply(pop_model_pred[,-1:-2], 1, max)
            pop_model_pred <- pop_model_pred[pop_model_pred\$max_prob > ${prob_thresh},]
            pop_model_pred\$max_prob <- NULL
        }
        
        # Generate population-specific eigenvec files and keep files
        for(i in names(pop_model_pred[,-1:-2])){
            tmp_keep <- pop_model_pred[apply(pop_model_pred[,-1:-2], 1, function(x) x[i] == max(x)),1:2]
            fwrite(tmp_keep, paste0("${output_name}.", i, ".keep"), sep=' ', col.names=F)
            
            # Generate population-specific scaled eigenvec files
            if(length(keep_files) > 0){
                # Find matching keep file for this population
                pop_keep_idx <- which(gsub("\\\\.keep\$", "", basename(keep_files)) == i)
                if(length(pop_keep_idx) > 0){
                    pop_scale_file <- paste0("${output_name}.", i, ".scale")
                    if(file.exists(pop_scale_file)){
                        pop_scale <- fread(pop_scale_file)
                        target_PCs_pop_scaled <- target_PCs
                        for(j in 1:dim(pop_scale)[1]){
                            target_PCs_pop_scaled[[paste0('PC',j)]] <- target_PCs[[paste0('PC',j)]]-pop_scale\$Mean[pop_scale\$PC == paste0('PC',j)]
                            target_PCs_pop_scaled[[paste0('PC',j)]] <- target_PCs_pop_scaled[[paste0('PC',j)]]/pop_scale\$SD[pop_scale\$PC == paste0('PC',j)]
                            target_PCs_pop_scaled[[paste0('PC',j)]] <- round(target_PCs_pop_scaled[[paste0('PC',j)]],3)
                        }
                        fwrite(target_PCs_pop_scaled, paste0("${output_name}.", i, ".eigenvec"), sep='\\t')
                    }
                }
            }
        }
        
        # Create plots
        sink(file = "${output_name}.log", append = T)
        cat('Creating PCA plots...\\n')
        sink()
        
        # Read in population data
        pop_data <- data.frame(fread("${pop_data}"))
        names(pop_data)[1] <- 'IID'
        pop_data\$FID <- pop_data\$IID
        
        # Merge reference PCs with population data
        ref_PCs_plot <- merge(PCs_ref, pop_data, by=c('FID','IID'))
        
        # Add target data
        new_cols <- names(ref_PCs_plot[!grepl('PC|ID', names(ref_PCs_plot))])
        new_cols_df <- data.frame(matrix(rep('Target',length(new_cols)),ncol=length(new_cols)))
        names(new_cols_df) <- new_cols
        target_PCs_plot <- cbind(target_PCs, new_cols_df)
        
        # Combine reference and target
        ref_target_PCs <- rbind(ref_PCs_plot, target_PCs_plot)
        
        # Create plots for each population column
        Label_groups <- names(ref_target_PCs[!grepl('PC|IID|FID',names(ref_target_PCs))])
        
        for(i in Label_groups){
            PC_1_2 <- ggplot(ref_target_PCs[ref_target_PCs[[i]] != 'Target',], aes(x=PC1,y=PC2, colour=get(i))) +
                geom_point() +
                geom_point(data=ref_target_PCs[ref_target_PCs[[i]] == 'Target',], aes(x=PC1,y=PC2), colour='black', shape=21) +
                ggtitle("PCs 1 and 2") +
                labs(colour="")
            
            PC_3_4 <- ggplot(ref_target_PCs[ref_target_PCs[[i]] != 'Target',], aes(x=PC3,y=PC4, colour=get(i))) +
                geom_point() +
                geom_point(data=ref_target_PCs[ref_target_PCs[[i]] == 'Target',], aes(x=PC3,y=PC4), colour='black', shape=21) +
                ggtitle("PCs 3 and 4") +
                labs(colour="")
            
            png(paste0("${output_name}.PCs_plot_", i, ".png"), units='px', res=300, width=4000, height=2500)
            print(plot_grid(PC_1_2,PC_3_4))
            dev.off()
        }
        
        sink(file = "${output_name}.log", append = T)
        cat('Plots created successfully\\n')
        sink()
    }
    
    # Write versions
    writeLines(c('"${task.process}":',
                 paste0('    r-base: "', R.version.string, '"'),
                 paste0('    r-data.table: "', packageVersion("data.table"), '"'),
                 paste0('    r-caret: "', packageVersion("caret"), '"'),
                 paste0('    r-ggplot2: "', packageVersion("ggplot2"), '"')), 
               "versions.yml")
    
    sink(file = "${output_name}.log", append = T)
    cat('PC Analysis completed at', as.character(Sys.time()), '\\n')
    sink()
    """
}
