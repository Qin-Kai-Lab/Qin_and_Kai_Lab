library(Seurat)
library(ggplot2)
library(monocle3)

targDir <- 'OPC_ODC/Monocle3/86PC/'

rnaDat = readRDS('integ_OPC_ODC')
rnaDat  = ScaleData(rnaDat)

genes = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv'))

genes$fdr_p = p.adjust(genes$p_value)

genes = genes[order(genes$estimate, decreasing = T),]

genesFilt = unique(genes$gene_id[genes$q_value < 0.05])

renameClusters = function(df, clustCol) {
  newClust = character()
  for ( i in 1:nrow(df)) {
    if (df[i, clustCol] == 1) {
      curClust = "ODC"
    } else if (df[i, clustCol] == 2) {
      curClust = "OPC"
    } else if ((df[i, clustCol] == 3) | (df[i, clustCol] == 4)) {
      curClust = "Intermideate"
    }
    newClust = c(newClust, curClust)
  }
  df$newMonocClust = newClust
  return(df)
}

rnaDat@meta.data = renameClusters(df = rnaDat@meta.data, clustCol = "MonocClust")

Idents(rnaDat) = rnaDat$newMonocClust

levels(x =rnaDat) <- c('OPC', 'Intermideate', 'ODC')

curPlot = DoHeatmap(rnaDat,  features = genesFilt, size = 10)  +
  theme(text = element_text(size = 24)) +
  #guides(color="none") +
  theme(legend.text = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 1)) +
  theme(legend.title = element_text(size = 24))


curPlot

targDir = './Paper_figs/Fig3/'

topN = "AllSign"
typeGenes = "RNA_Linear"
png(filename = paste0(targDir,"TimeExpression_Heatmap_", topN, "_", typeGenes, "_2023-12-19.png"), width = 22, height = 22, units = "in", res = 300)
print(curPlot )
dev.off()


curPlot<-DoHeatmap(rnaDat, features = genesFilt,  label = F)+
  theme(text = element_text(size = 24)) +
  #guides(color="none") +
  theme(legend.text = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 1)) +
  theme(legend.title = element_text(size = 24))

topN = "AllSign"
typeGenes = "RNA_Linear"
png(filename = paste0(targDir,"TimeExpression_Heatmap_Legend_", topN, "_", typeGenes, "_2023-12-19.png"), width = 22, height = 22, units = "in", res = 300)
print(curPlot )
dev.off()