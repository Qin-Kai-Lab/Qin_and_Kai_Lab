library(ggplot2)
library(stringr)
library(AUCell)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

curDate<-Sys.Date()

levels(RNA.combined.norm)

DefaultAssay(RNA.combined.norm)

# list of kegg genes 
kegg<-read.delim('../data/kegg_mouse_2019.tsv', T, sep='\t')

# get gene matrix from  seurat

dimNames<-RNA.combined.norm@assays$RNA@data@Dimnames
all_genes<-dimNames[[1]]

#counts<-GetAssayData(RNA.combined.norm, slot= 'counts')

countsMat<-as.matrix(GetAssayData(RNA.combined.norm, slot= 'counts'))

# workflow

kegg_path<-kegg

kegg_list<-list()

for ( i in 1:nrow(kegg_path)){
  genes<-strsplit(kegg_path$Genes[i], "\t")
  
  genesFromat<-list(str_to_title(genes[[1]]))
  
  kegg_list<-c(kegg_list, genesFromat)
}

names(kegg_list)<-kegg_path$Pathway

#cell_ranking<- AUCell_buildRankings(countsMat)

cells_AUC <- AUCell_run(countsMat, kegg_list)


save(cells_AUC, file="cellsAUC_kegg2019Mouse.RData")


# threshholds

#set.seed(13)
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE, nCores =4) 

#save(cells_assignment, file="cellAssignmet_kegg2019Mouse.RData")


#selectedThresholds <- getThresholdSelected(cells_assignment)

AUCmat <- AUCell::getAUC(cells_AUC)
