# Create Individual IBD outlier Histogram pdf

args <- commandArgs(TRUE)

library(data.table)

#Read data
setDTthreads(8)

# Find the .IBD_INDIV.txt file in the current directory
indiv_files <- list.files(pattern = "\\.IBD_INDIV\\.txt$")
if (length(indiv_files) == 0) {
    stop("No .IBD_INDIV.txt files found in current directory")
}

# Use the first .IBD_INDIV.txt file found
indiv_file <- indiv_files[1]
IndvIBD<-read.table(indiv_file, head=TRUE)

pdf(file=paste0("IndvIBD_hist_",args[2],".pdf"))
hist(IndvIBD$MEAN_PI_HAT, breaks=100)
dev.off()
