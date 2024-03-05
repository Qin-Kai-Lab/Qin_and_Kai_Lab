
library(AUCell)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


# get count matrix
source('../../programs/renameClusters.R')

curDate<-Sys.Date()

levels(RNA.combined.norm)

DefaultAssay(RNA.combined.norm)


countsMat<-as.matrix(GetAssayData(RNA.combined.norm, slot= 'counts'))

rm(RNA.combined.norm)

gc()

# get named list of genes

kegg<-read.csv('../data/keggMouse_clustprof_db.csv')

pathList<-unique(kegg$Kegg_path)

keggList<-list()

for (i in pathList){
  kPath<-kegg[(kegg$Kegg_path==i),]
  kgenes<-list(unique(kPath$Gene))
  names(kgenes)<-i
  keggList<-c(keggList, kgenes)
}

cells_AUC <- AUCell_run(countsMat, keggList)


save(cells_AUC, file="cellsAUC_keggClustProf.RData")