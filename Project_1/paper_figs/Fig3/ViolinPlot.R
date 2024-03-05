library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)
library(limma)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig3/ViolinPlots/'

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

genes<-c('Plp1', 'Ptgds', 'Trf', 'Calcrl', 'Cspg4', 'Vcan')

vnPlot<-VlnPlot(
  object = RNA.combined.norm,
  features = genes,
  ncol = 3,
  pt.size = 0.2
)

ggsave(paste0(targDir,'VnPlots_monoc_3_Clust_CustGenes_', curDate, '.jpeg'), plot = vnPlot, height = 14, width = 16, units = 'in', dpi = 300)

# top genes folowwing pseudotime

allMarkers = read.csv('OPC_ODC/Monocle3/86PC/OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv')

# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$cluster)
  markers = character()
  minNumb = ncol(RNA.combined.norm) * 0.1
  dfPct = df[(df$num_cells_expressed > minNumb), ]
  dfSort = dfPct[order(dfPct$q_value, decreasing = F),]
  markers = head(dfSort, topNumb)
  unMark = unique(markers$gene_id)
  return(unMark)
}

topMark = getTopMarkers(df = allMarkers, topNumb = 12)

levels(x =RNA.combined.norm) <- c('OPC', 'Intermideate', 'ODC')

RNA.combined.norm$newMonocClust = factor(RNA.combined.norm$newMonocClust, levels = c('OPC', 'Intermideate', 'ODC'))

Idents(RNA.combined.norm) = RNA.combined.norm$newMonocClust

vnPlot<-VlnPlot(
  object = RNA.combined.norm,
  features = topMark,
  ncol = 3,
  pt.size = 0.1,
  group.by = 'newMonocClust'
)

ggsave(paste0(targDir,'VnPlots_monoc_3_Clust_Top12TimeGenes_', curDate, '.jpeg'), plot = vnPlot, height = 16, width = 18, units = 'in', dpi = 300)