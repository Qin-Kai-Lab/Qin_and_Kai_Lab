library(AUCell)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

# upload data and set variables
source('../../programs/renameClusters.R')

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

# load enrichment results

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AUCmat)
  DefaultAssay(RNA.combined.norm) <- 'AUC'
  RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))
  return(RNA.combined.norm)
}

# add AUcell data
RNA.combined.norm<-addEnrich(x='./cellsAUC_keggClustProf.RData')

DefaultAssay(RNA.combined.norm)

# add differential expression data

deDat<-read.csv('./Enrichment/AUcell/Kegg_ClustProfiler/pos_and_neg/All_ContrVsStress_2022-11-28.csv')

deDat$Log_abs<-abs(deDat$avg_log2FC)

getTopMatch<-function(x){
  combTop<-data.frame(matrix(nrow=0, ncol = 0))
  for ( i in clusters){
    dfSel<-x[(x$Cell_Type==i),]
    dfOrder<-dfSel[order(dfSel$Log_abs, decreasing = T),]
    dfTop<-head(dfOrder, 10)
    combTop<-rbind(combTop, dfTop)
  }
  return(combTop)
}

topMatch<-getTopMatch(x=deDat)

testMarkers<-unique(topMatch$Genes)

markersList<-split(testMarkers, ceiling(seq_along(testMarkers)/4))

exprPlot<-function(x) {
  for (i in 1:length(x)) {
    genes<-x[[i]]
    fPlot<-FeaturePlot(RNA.combined.norm, features = genes, min.cutoff = "q10", split.by = 'group', ncol = 4)
    ggsave(paste0('FP_AuClustKeg_list_', i,"_", curDate, ".jpeg"), plot=fPlot, height = 20, width = 14, units = 'in', dpi = 300)
  }
}

exprPlot(markersList)


violPlot<-function(x) {
  for (i in 1:length(x)) {
    genes<-x[[i]]
    fPlot<-VlnPlot(RNA.combined.norm, features = genes, pt.size = 0, split.by = 'group')
    ggsave(paste0('VP_byGr_AuClustKeg_list_', i,"_", curDate, ".jpeg"), plot=fPlot, height = 10, width = 22, units = 'in', dpi = 300)
  }
}

violPlot(markersList)