#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)

#Read data
setDTthreads(8)

# Find the .IBD.genome file in the current directory
genome_files <- list.files(pattern = "\\.IBD\\.genome$")
if (length(genome_files) == 0) {
    stop("No .IBD.genome files found in current directory")
}

# Use the first .IBD.genome file found
genome_file <- genome_files[1]
IBD<-fread(genome_file, data.table=F, select=c("IID1", "IID2", "PI_HAT"))

pdf(file=paste0("IBD_hist_",args[2],".pdf"))
hist(IBD$PI_HAT, 100, ylim=c(0,100))
