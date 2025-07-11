#Create Sex Check Histogram pdfs

args <- commandArgs(TRUE)

setwd(args[3])

Sexcheck <-read.table(args[1], head=T)
Sexcheck_2 <-Sexcheck[Sexcheck$PEDSEX == 2, ]
Sexcheck_1 <-Sexcheck[Sexcheck$PEDSEX == 1, ]
Sexcheck_0 <-Sexcheck[Sexcheck$PEDSEX == 0, ]

# Check if there's any data to plot
if (nrow(Sexcheck_2) > 0 || nrow(Sexcheck_1) > 0 || nrow(Sexcheck_0) > 0) {
  pdf(file=paste0("sexcheck_hist_",args[2],".pdf"))
  
  # Plot first non-empty group to establish the plot
  first_plot <- TRUE
  if (nrow(Sexcheck_2) > 0) {
    hist(Sexcheck_2$F, 100, col=rgb(1,0,0,0.5))
    first_plot <- FALSE
  }
  if (nrow(Sexcheck_1) > 0) {
    hist(Sexcheck_1$F, 100, col=rgb(0,1,0,0.5), add=!first_plot)
    first_plot <- FALSE
  }
  if (nrow(Sexcheck_0) > 0) {
    hist(Sexcheck_0$F, 100, col=rgb(0,0,1,0.5), add=!first_plot)
  }
  
  while (!is.null(dev.list()))  dev.off()

  pdf(file=paste0("sexcheck_hist_y20_",args[2],".pdf"))
  
  # Plot first non-empty group to establish the plot with ylim
  first_plot <- TRUE
  if (nrow(Sexcheck_2) > 0) {
    hist(Sexcheck_2$F, 100, col=rgb(1,0,0,0.5), ylim=c(0,20))
    first_plot <- FALSE
  }
  if (nrow(Sexcheck_1) > 0) {
    hist(Sexcheck_1$F, 100, col=rgb(0,1,0,0.5), add=!first_plot, ylim=if(first_plot) c(0,20) else NULL)
    first_plot <- FALSE
  }
  if (nrow(Sexcheck_0) > 0) {
    hist(Sexcheck_0$F, 100, col=rgb(0,0,1,0.5), add=!first_plot, ylim=if(first_plot) c(0,20) else NULL)
  }
  
  while (!is.null(dev.list()))  dev.off()
}
