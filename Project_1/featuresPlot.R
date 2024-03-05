library(Seurat)
library(ggplot2)
library(stringr)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

RNA.combined.norm<-readRDS('integRNADoublFilt')
DefaultAssay(RNA.combined.norm) <- "RNA"
DefaultAssay(RNA.combined.norm)
Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

genes<-as.character(t(read.table('../data/features_10-20-22.txt', F, sep=',')))

genes<-gsub(' ', '', genes)

dimNames<-RNA.combined.norm@assays$RNA@data@Dimnames
all_genes<-dimNames[[1]]

markersPresent<-genes[(genes%in%all_genes)]

missGenes<-genes[!(genes%in%all_genes)]

"Mtmr3"%in%all_genes

testMarkers<-c(markersPresent, "Aldh1l1", "Olig2", "Grm3", "Mtmr3", "Arhgap45")

markersList<-split(testMarkers, ceiling(seq_along(testMarkers)/6))

curDate<-Sys.Date()

# plot function
exprPlot<-function(x) {
  for (i in 1:length(x)) {
    genes<-x[[i]]
    fPlot<-FeaturePlot(RNA.combined.norm, features = genes, min.cutoff = "q9")
    ggsave(paste0('featurePlot_list_', i,"_", curDate, ".jpeg"), plot=fPlot, height = 10, width = 12, units = 'in', dpi = 300)
  }
}

exprPlot(markersList)

# additional genes

extraGenes<- c('Fn1', 'Slc17a6', 'Meis2', 'Rxfp1', 'Nts')

extraPresent<-extraGenes[(extraGenes%in%all_genes)]

fPlot<-FeaturePlot(RNA.combined.norm, features = extraGenes, min.cutoff = "q9")
ggsave(paste0('featurePlot_list_', "12" ,"_", curDate, ".jpeg"), plot=fPlot, height = 10, width = 12, units = 'in', dpi = 300)

# one gene
fPlot<-FeaturePlot(RNA.combined.norm, features = 'Sema6d', min.cutoff = "q9")
ggsave(paste0('featurePlot_Sema6d_', curDate, ".jpeg"), plot=fPlot, height = 10, width = 12, units = 'in', dpi = 300)