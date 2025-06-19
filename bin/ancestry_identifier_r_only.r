#!/usr/bin/Rscript
# This script contains only the R analysis parts of ancestry_identifier.r
# PLINK operations have been moved to separate Nextflow processes

start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_scores", action="store", default=NA, type='character',
    help="Path to reference PC scores file [required]"),
make_option("--target_scores", action="store", default=NA, type='character',
    help="Path to target PC scores file [required]"),
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
    help="Indicates whether probability threshold should be used when defining ancestry [optional]")
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
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Ancestry_identifier_r_only.R
# R-only analysis for ancestry identification
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

# Read in reference PC scores
PCs_ref <- data.frame(fread(opt$ref_scores))
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
# Process target sample scores
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Processing target sample PCs...')
sink()

# Read target scores
if(file.exists(opt$target_scores)) {
  targ_PCs <- data.frame(fread(opt$target_scores))
  # Adjust column selection based on actual file structure
  if(ncol(targ_PCs) >= (2 + opt$n_pcs)) {
    targ_PCs <- targ_PCs[,c(1:2, (ncol(targ_PCs)-opt$n_pcs+1):ncol(targ_PCs))]
    names(targ_PCs) <- c('FID','IID',paste0('PC',1:as.numeric(opt$n_pcs)))
  } else {
    stop("Target scores file does not have expected number of columns")
  }
} else {
  stop("Target scores file not found")
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

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

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
