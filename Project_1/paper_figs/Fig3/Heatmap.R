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
RNA.combined.norm = ScaleData(RNA.combined.norm)

levels(x =RNA.combined.norm) <- c('OPC', 'Intermideate', 'ODC')

allMarkers = read.csv(paste0(targDir, "All_RNA_Markers_2023-05-08.csv" ))


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
  curHeat<-DoHeatmap(proj, features = x)+
    theme(text = element_text(size = 20)) +
    guides(color="none") +
    theme(legend.text = element_text(size = 18)) +
    theme(axis.text.y = element_text(size = 18)) +
    theme(legend.title = element_text(size = 18))
  
  ggsave(paste0(targDir,'Heatmap', genesSet, curDate, '.jpeg'), plot =  curHeat, height = 16, width = 20, units = 'in', dpi = 300)
}

makePlot(x = topMark, genesSet = "_Top10Markers_", proj = RNA.combined.norm)
