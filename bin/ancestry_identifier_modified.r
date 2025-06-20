#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
# Modified to separate PLINK/PLINK2 calls from R functionality

start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--target_plink", action="store", default=NA, type='character',
              help="Path to target PLINK files [required]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--n_pcs", action="store", default=10, type='numeric',
              help="Number of PCs (min=4) [optional]"),
  make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
              help="Path for output files [required]"),
  make_option("--ref_pop_scale", action="store", default=NA, type='character',
              help="List of keep files ancestry specific scaling [optional]"),
  make_option("--pop_data", action="store", default=NA, type='character',
              help="Population data for the reference samples [required]"),
  make_option("--model_method", action="store", default='glmnet', type='character',
              help="Method used for generate prediction model [optional]"),
  make_option("--SD_rule", action="store", default=F, type='logical',
              help="Logical indicating whether the 3SD rule should be used to define ancestry, or the model-based approach [optional]"),
  make_option("--prob_thresh", action="store", default='NA', type='numeric',
              help="Indicates whether probability threshold should be used when defining ancestry [optional]"),
  make_option("--profiles_sscore", action="store", default=NA, type='character',
              help="Path to the profiles.sscore file generated by PLINK2 [required]"),
  make_option("--eigenvec_allele", action="store", default=NA, type='character',
              help="Path to the eigenvec.allele file generated by PLINK2 [required]"),
  make_option("--ref_merge_bim", action="store", default=NA, type='character',
              help="Path to the ref_merge.bim file [required]"),
  make_option("--ref_merge_pruned_score_sscore", action="store", default=NA, type='character',
              help="Path to the ref_merge_pruned_score.sscore file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(caret)
library(pROC)
library(verification)
library(ggplot2)
library(cowplot)

setDTthreads(8)

opt$output_dir <- paste0(dirname(opt$output),'/')
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Ancestry_identifier.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

if(is.na(opt$ref_pop_scale)){
  sink(file = paste(opt$output,'.log',sep=''), append = F)
  cat('ref_pop must be specified\n')
  sink()
  q()
}

# Process reference BIM file for allele matching and high LD regions
ref_bim <- data.frame(fread(opt$ref_merge_bim))

# Create file removing high LD regions
long_ld_exclude <- ref_bim$V2[ (ref_bim$V1 == 1 & ref_bim$V4 >= 48e6 & ref_bim$V4 <= 52e6) |
                                (ref_bim$V1 == 2 & ref_bim$V4 >= 86e6 & ref_bim$V4 <= 100.5e6) |
                                (ref_bim$V1 == 2 & ref_bim$V4 >= 134.5e6 & ref_bim$V4 <= 138e6) |
                                (ref_bim$V1 == 2 & ref_bim$V4 >= 183e6 & ref_bim$V4 <= 190e6) |
                                (ref_bim$V1 == 3 & ref_bim$V4 >= 47.5e6 & ref_bim$V4 <= 50e6) |
                                (ref_bim$V1 == 3 & ref_bim$V4 >= 83.5e6 & ref_bim$V4 <= 87e6) |
                                (ref_bim$V1 == 3 & ref_bim$V4 >= 89e6 & ref_bim$V4 <= 97.5e6) |
                                (ref_bim$V1 == 5 & ref_bim$V4 >= 44.5e6 & ref_bim$V4 <= 50.5e6) |
                                (ref_bim$V1 == 5 & ref_bim$V4 >= 98e6 & ref_bim$V4 <= 100.5e6) |
                                (ref_bim$V1 == 5 & ref_bim$V4 >= 129e6 & ref_bim$V4 <= 132e6) |
                                (ref_bim$V1 == 5 & ref_bim$V4 >= 135.5e6 & ref_bim$V4 <= 138.5e6) |
                                (ref_bim$V1 == 6 & ref_bim$V4 >= 25.5e6 & ref_bim$V4 <= 33.5e6) |
                                (ref_bim$V1 == 6 & ref_bim$V4 >= 57e6 & ref_bim$V4 <= 64e6) |
                                (ref_bim$V1 == 6 & ref_bim$V4 >= 140e6 & ref_bim$V4 <= 142.5e6) |
                                (ref_bim$V1 == 7 & ref_bim$V4 >= 55e6 & ref_bim$V4 <= 66e6) |
                                (ref_bim$V1 == 8 & ref_bim$V4 >= 8e6 & ref_bim$V4 <= 12e6) |
                                (ref_bim$V1 == 8 & ref_bim$V4 >= 43e6 & ref_bim$V4 <= 50e6) |
                                (ref_bim$V1 == 8 & ref_bim$V4 >= 112e6 & ref_bim$V4 <= 115e6) |
                                (ref_bim$V1 == 10 & ref_bim$V4 >= 37e6 & ref_bim$V4 <= 43e6) |
                                (ref_bim$V1 == 11 & ref_bim$V4 >= 46e6 & ref_bim$V4 <= 57e6) |
                                (ref_bim$V1 == 11 & ref_bim$V4 >= 87.5e6 & ref_bim$V4 <= 90.5e6) |
                                (ref_bim$V1 == 12 & ref_bim$V4 >= 33e6 & ref_bim$V4 <= 40e6) |
                                (ref_bim$V1 == 12 & ref_bim$V4 >= 109.5e6 & ref_bim$V4 <= 112e6) |
                                (ref_bim$V1 == 20 & ref_bim$V4 >= 32e6 & ref_bim$V4 <= 34.5e6)]

write.table(long_ld_exclude, paste0(opt$output_dir,'long_ld.exclude'), col.names=F, row.names=F, quote=F)

# Read in reference PC scores
PCs_ref <- data.frame(fread(opt$ref_merge_pruned_score_sscore))
PCs_ref <- PCs_ref[,c(1:2,5:dim(PCs_ref)[2])]
names(PCs_ref) <- c('FID','IID',paste0('PC',1:as.numeric(opt$n_pcs)))

fwrite(PCs_ref, paste0(opt$output,'.eigenvec'), sep=' ')

# Scale across all individuals
PCs_ref_centre_scale <- data.frame(PC=names(PCs_ref[-1:-2]),
                                  Mean=sapply(PCs_ref[,-1:-2], function(x) mean(x)),
                                  SD=sapply(PCs_ref[,-1:-2], function(x) sd(x)),
                                  row.names=seq(1:as.numeric(opt$n_pcs)))

fwrite(PCs_ref_centre_scale, paste0(opt$output,'.scale'), sep=' ')

gc()

if(!is.na(opt$ref_pop_scale)){
  # Calculate the mean and sd of scores for each population specified in pop_scale
  pop_keep_files <- read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

  for(k in 1:dim(pop_keep_files)[1]){
    pop <- pop_keep_files$V1[k]
    keep <- fread(pop_keep_files$V2[k], header=F)
    PCs_ref_keep <- PCs_ref[(PCs_ref$FID %in% keep$V1),]

    PCs_ref_centre_scale <- data.frame(PC=names(PCs_ref_keep[-1:-2]),
                                      Mean=sapply(PCs_ref_keep[,-1:-2], function(x) mean(x)),
                                      SD=sapply(PCs_ref_keep[,-1:-2], function(x) sd(x)),
                                      row.names=seq(1:opt$n_pcs))

    fwrite(PCs_ref_centre_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')

    rm(PCs_ref_centre_scale)
    gc()
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Create model predicting ref_pop groups
###

if(!is.na(opt$ref_pop_scale)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Deriving model predicting ref_pop groups...')
  sink()

  # Read in whole sample scale file
  PCs_ref_centre_scale <- fread(paste0(opt$output,'.scale'))

  # Scale the reference PCs
  PCs_ref_scaled <- PCs_ref
  for(i in 1:dim(PCs_ref_centre_scale)[1]){
    PCs_ref_scaled[[paste0('PC',i)]] <- PCs_ref[[paste0('PC',i)]]-PCs_ref_centre_scale$Mean[PCs_ref_centre_scale$PC == paste0('PC',i)]
    PCs_ref_scaled[[paste0('PC',i)]] <- PCs_ref_scaled[[paste0('PC',i)]]/PCs_ref_centre_scale$SD[PCs_ref_centre_scale$PC == paste0('PC',i)]
  }

  # Label individuals with ref_pop groups
  pop <- NULL
  for(i in 1:dim(pop_keep_files)[1]){
    keep <- fread(pop_keep_files$V2[i], header=F)
    keep$pop <- pop_keep_files$V1[i]
    pop <- rbind(pop,keep)
  }
  names(pop) <- c('FID','IID','pop')
  PCs_ref_scaled_pop <- merge(PCs_ref_scaled,pop, by=c('FID','IID'))
  rm(PCs_ref_scaled)
  gc()

  # Build model
  model <- train(y=as.factor(PCs_ref_scaled_pop$pop), x=PCs_ref_scaled_pop[grepl('PC',names(PCs_ref_scaled_pop))], method=opt$model_method, metric='logLoss', trControl=trainControl(method="cv", number=5, classProbs= TRUE, savePredictions = 'final', summaryFunction = multiClassSummary))

  # Save performance information
  sink(file = paste(opt$output,'.pop_model_prediction_details.txt',sep=''), append = F)
  print(model)
  cat('\n')
  obs_pre_tab <- table(model$pred$obs, model$pred$pred)
  dimnames(obs_pre_tab) <- list(paste('obs',dimnames(obs_pre_tab)[[1]]),paste('pred',dimnames(obs_pre_tab)[[2]]))

  # Show confusion matrix before and after applying probability threshold
  cat('Confusion matrix without threshold:\n')
  print(obs_pre_tab)

  if(!is.na(opt$prob_thresh)){
    model$pred$max_prob <- apply(model$pred[,unique(PCs_ref_scaled_pop$pop)], 1, max)
    model$pred <- model$pred[model$pred$max_prob > opt$prob_thresh,]

    obs_pre_tab_thresh <- table(model$pred$obs, model$pred$pred)
    dimnames(obs_pre_tab_thresh) <- list(paste('obs',dimnames(obs_pre_tab_thresh)[[1]]),paste('pred',dimnames(obs_pre_tab_thresh)[[2]]))

    cat('\n')
    cat(paste0('Confusion matrix with ',opt$prob_thresh,' threshold:\n'))
    print(obs_pre_tab_thresh)
  }

  sink()

  saveRDS(model$finalModel, paste0(opt$output,'.pop_model.rds'))

  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Done!\n')
  sink()
}

#####
# Process target sample PCs
#####

# Read target profiles
scores <- fread(cmd=paste0('cut -f 1-2 ',opt$profiles_sscore))
names(scores) <- c('FID','IID')

profile_all <- data.frame(fread(opt$profiles_sscore))
profile_all <- as.matrix(profile_all[,grepl('PC',names(profile_all))])

profile_all <- data.table(profile_all)
names(profile_all) <- paste0('PC',1:as.numeric(opt$n_pcs))
scores <- cbind(scores, profile_all)

targ_PCs <- data.frame(scores)
rm(scores, profile_all)
gc()

###
# Create plot PC scores of target sample compared to the reference
###
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Plotting target sample PCs on reference...')
sink()

# Read in population data
pop_data <- data.frame(fread(opt$pop_data))
names(pop_data)[1] <- 'IID'
pop_data$FID <- pop_data$IID

# Read in reference sample PCs
ref_PCs <- data.frame(fread(paste0(opt$output,'.eigenvec')))
ref_PCs <- merge(ref_PCs, pop_data, by=c('FID','IID'))

# Insert pop_data columns into target PCs
new_cols <- names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
new_cols_2 <- data.frame(matrix(rep('Target',length(new_cols)),ncol=length(new_cols)))
names(new_cols_2) <- names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
targ_PCs <- cbind(targ_PCs,new_cols_2)

# Combine the two sets
ref_PCs_targ_PCs <- rbind(ref_PCs,targ_PCs)

rm(ref_PCs)
gc()

Label_groups <- names(ref_PCs_targ_PCs[!grepl('PC|IID|FID',names(ref_PCs_targ_PCs))])

for(i in Label_groups){
  PC_1_2 <- ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC1,y=PC2, colour=get(i))) +
    geom_point() +
    geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC1,y=PC2), colour='black', shape=21) +
    ggtitle("PCs 1 and 2") +
    labs(colour="")
  PC_3_4 <- ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC3,y=PC4, colour=get(i))) +
    geom_point() +
    geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC3,y=PC4), colour='black', shape=21) +
    ggtitle("PCs 3 and 4") +
    labs(colour="")

  png(paste0(opt$output,'.PCs_plot_',i,'.png'), units='px', res=300, width=4000, height=2500)
  print(plot_grid(PC_1_2,PC_3_4))
  dev.off()

  rm(PC_1_2,PC_3_4)
  gc()

  print(i)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Estimate probability of outcomes in model
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Estimating probability of each population...')
sink()

# Read in the reference scale file
pop_model_scale <- fread(paste0(opt$output,'.scale'))

# Scale the target PCs
targ_PCs_scaled <- targ_PCs
for(i in 1:dim(pop_model_scale)[1]){
  targ_PCs_scaled[[paste0('PC',i)]] <- targ_PCs[[paste0('PC',i)]]-pop_model_scale$Mean[pop_model_scale$PC == paste0('PC',i)]
  targ_PCs_scaled[[paste0('PC',i)]] <- targ_PCs_scaled[[paste0('PC',i)]]/pop_model_scale$SD[pop_model_scale$PC == paste0('PC',i)]
}

# Read in model
pop_model <- readRDS(paste0(opt$output,'.pop_model.rds'))
pop_model_pred <- predict(object = pop_model, newx = data.matrix(targ_PCs_scaled[grepl('PC',names(targ_PCs_scaled))]), type = "response", s=pop_model$lambdaOpt)
pop_model_pred <- as.data.frame.table(pop_model_pred)
pop_model_pred <- data.table(FID=targ_PCs_scaled$FID,
                             IID=targ_PCs_scaled$IID,
                             pop=as.character(pop_model_pred$Var2),
                             prob=round(pop_model_pred$Freq,3))

pop_model_pred <- dcast.data.table(pop_model_pred, formula=FID + IID~pop, value.var = "prob")

fwrite(pop_model_pred, paste0(opt$output,'.model_pred'), sep='\t')

# Create keep files based on the results
if(!is.na(opt$prob_thresh)){
  pop_model_pred$max_prob <- apply(pop_model_pred[,-1:-2], 1, max)
  pop_model_pred <- pop_model_pred[pop_model_pred$max_prob > opt$prob_thresh,]
  pop_model_pred$max_prob <- NULL
}

N_group <- NULL
for(i in names(pop_model_pred[,-1:-2])){
  tmp_keep <- pop_model_pred[apply(pop_model_pred[,-1:-2], 1, function(x) x[i] == max(x)),1:2]
  N_group <- rbind(N_group, data.frame(Group=i, N=dim(tmp_keep)[1]))
  fwrite(tmp_keep, paste0(opt$output,'.model_pred.',i,'.keep'), sep=' ', col.names=F)
}

rm(targ_PCs_scaled,pop_model_pred)
gc()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('----------\n')
cat('N per group based on model:\n')
print(N_group)
cat('----------\n')
sink()

###
# Identify individuals that are within 3SD of population specific mean for all PCs and write out scaled PCs.
###

targ_PCs <- targ_PCs[,grepl('FID|IID|PC',names(targ_PCs))]
# Read in pop_scale_for_keep
pop_scale_for_keep <- paste0(opt$output,'.',pop_keep_files$V1,'.scale')

N_group <- NULL
for(i in 1:length(pop_scale_for_keep)){
  # Idenitfy name of population based on scale file
  pop_name <- gsub('.scale','',substr(pop_scale_for_keep[i], nchar(pop_scale_for_keep[i])-9+1, nchar(pop_scale_for_keep[i])))
  pop_scale_for_keep_i <- fread(pop_scale_for_keep[i])

  # Scale the target based on population scale file
  targ_PCs_scaled_i <- targ_PCs[,grepl('FID|IID|PC', names(targ_PCs))]
  for(j in 1:dim(pop_scale_for_keep_i)[1]){
    targ_PCs_scaled_i[[paste0('PC',j)]] <- targ_PCs[[paste0('PC',j)]]-pop_scale_for_keep_i$Mean[pop_scale_for_keep_i$PC == paste0('PC',j)]
    targ_PCs_scaled_i[[paste0('PC',j)]] <- targ_PCs_scaled_i[[paste0('PC',j)]]/pop_scale_for_keep_i$SD[pop_scale_for_keep_i$PC == paste0('PC',j)]
    targ_PCs_scaled_i[[paste0('PC',j)]] <- round(targ_PCs_scaled_i[[paste0('PC',j)]],3)
  }

  # Remove anyone with a PC value >3 or -3 (i.e. 3SD from the population mean
  # targ_PCs_scaled_i <- targ_PCs_scaled_i[!apply(targ_PCs_scaled_i[,-1:-2], 1, function(x) any(x > 3 | x < -3)),]

  N_group <- rbind(N_group, data.frame(Group=pop_name, N=dim(targ_PCs_scaled_i)[1]))

  # Save keep file of individuals that fit the population
  fwrite(targ_PCs_scaled_i[,1:2], paste0(opt$output,'.',pop_name,'.keep'), col.names=F, sep='\t')

  # Write the scaled PCs
  fwrite(targ_PCs_scaled_i, paste0(opt$output,'.',pop_name,'.eigenvec'), sep='\t')

  rm(pop_name,pop_scale_for_keep_i,targ_PCs_scaled_i)
  gc()
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('----------\n')
cat('N per group based on 3SD rule:\n')
print(N_group)
cat('----------\n')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
