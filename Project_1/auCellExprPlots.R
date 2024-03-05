library(AUCell)
library(limma)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

curDate<-Sys.Date()

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

# load enrichment results

load('./cellsAUC_keggClustProf.RData')

AUCmat <- AUCell::getAUC(cells_AUC)

RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AUCmat)
DefaultAssay(RNA.combined.norm) <- 'AUC'

RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))

kegg<-rownames(AUCmat)

#allHeat<-DoHeatmap(RNA.combined.norm, features = kegg)  + NoLegend()

#ggsave(paste0('./heatmap_aucKeggClustProf_', curDate, '.jpeg'), plot = allHeat, height = 35, width = 35, units = 'in', dpi = 300)



# select top pathways for each cluster 

allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = T, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox")

#write.csv(allMarkers, paste0('auCellPosMarkersClustProf_',curDate, '.csv'), row.names = F)

allMarkersSign<-allMarkers[(allMarkers$p_val_adj<0.05),]

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')


topMarkers<-character()

extractTopMarkers<-function(x){
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$AbsLog2FC, decreasing = T),]
    top10<-head(dfSort, 10)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}
top10MarkersDf<-extractTopMarkers(allMarkersSign)

topMarkers<-unique(top10MarkersDf$gene)

allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers) + NoLegend()

ggsave(paste0('./auCellTop10PosMarkersClustProf_', curDate, '.jpeg'), plot = allHeat, height = 16, width = 18, units = 'in', dpi = 300)