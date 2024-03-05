library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)
library(limma)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig3/'

dir.create(targDir, recursive = T, showWarnings = F)

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DimPlot(RNA.combined.norm, group.by = "MonocClust")

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"

Idents(RNA.combined.norm) = RNA.combined.norm$newMonocClust

levels(x =RNA.combined.norm) <- c('OPC', 'Intermideate', 'ODC')

allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

#write.csv(allMarkers, paste0(targDir, "All_RNA_Markers_", curDate, ".csv"), row.names = F)


clusters = c('OPC', 'Intermideate', 'ODC')
# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$cluster)
  markers = character()
  dfPct = df[(df$pct.1 >0.25) | (df$pct.2 >0.25), ]
  for (cluster in clusters) {
    dfSub =  dfPct[(dfPct[['cluster']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$avg_log2FC, decreasing = T),]
    topMark = head(dfOrd$gene, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}

topMark = getTopMarkers(df = allMarkers, topNumb = 10)

makePlot<-function(x, genesSet, proj) {
  dotPlot<-DotPlot(object = proj, features = x, scale.max = 100, dot.scale = 16)+
    scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(text = element_text(size = 24))+ # all text elements size
    theme(axis.text = element_text(size = 24)) 
  #scale_y_discrete(limits=rev)
  
  ggsave(paste0(targDir,'DotPlot', genesSet, '2023-12-14.jpeg'), plot =  dotPlot, height = 6, width = 20, units = 'in', dpi = 300)
}

makePlot(x = topMark, genesSet = "_Top10PosMarkers_", proj = RNA.combined.norm)

# AU cell

# load enrichment results
opcOdcCells = rownames(RNA.combined.norm@meta.data)

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  AucSub = AUCmat[, c(opcOdcCells)]
  RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AucSub)
  DefaultAssay(RNA.combined.norm) <- 'AUC'
  RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))
  return(RNA.combined.norm)
}

# add AUcell data
RNA.combined.norm<-addEnrich(x='./cellsAUC_keggClustProf.RData')
DefaultAssay(RNA.combined.norm)

allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox")

write.csv(allMarkers, paste0(targDir, "All_AUCell_KeggClustProf_Markers_", curDate, ".csv"), row.names = F)

topMark = getTopMarkers(df = allMarkers, topNumb = 8)

makePlot<-function(x, genesSet, proj) {
  dotPlot<-DotPlot(object = proj, features = x , scale.max = 100, dot.scale = 16)+
    scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(text = element_text(size = 16))+ # all text elements size
    theme(axis.text = element_text(size = 16)) + 
    theme(axis.text.x = element_text(size = 14)) +
    scale_x_discrete(limits=rev)
  
  ggsave(paste0(targDir,'DotPlot', genesSet, curDate, '.jpeg'), plot =  dotPlot, height = 12, width = 24, units = 'in', dpi = 300)
}

makePlot(x = topMark, genesSet = "_Top8PosMarkers_AUCell_KeggClustProf_", proj = RNA.combined.norm)
