
arg <- commandArgs(TRUE)

#Add al r scripts here - make one big pdf of all outputs, instead of separate ones?


pdf(paste(arg[2],'.pdf',sep=""),onefile=T,width = 8.3, height = 11.7)
layout(matrix(c(1,1,1,2,2,2,3,4,4), 3, 3, byrow = TRUE),widths=c(2,1,1), heights=c(1,1,1))
