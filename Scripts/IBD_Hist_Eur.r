library(data.table)

#Read data
setDTthreads(8)
IBD<-fread("COPINGR1R2R3_v1v2_b38_EUR_maf1_sample95.SNP95.hwe10.LD_Pruned_autosomalchr.IBD.genome", data.table=F, select=c("IID1", "IID2", "PI_HAT"))
