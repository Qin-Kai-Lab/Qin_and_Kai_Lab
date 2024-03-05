library(Seurat)
library(Signac)
library(ggplot2)


atacFilt = readRDS("atacIntegrated_macs2_2_RNA")
DefaultAssay(atacFilt) = "RNA"
corDf = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/PredictGene_RNA_spearman_AllClusters_2023-10-10.csv")

makeHeat = function(atacFilt, cluster,  topN, corDf) {
  corDfSel = corDf[corDf$Cluster == cluster,]
  corDfSel = corDfSel[corDfSel$Estimate > 0 & corDfSel$fdr_p < 0.05,]
  corDfSel = corDfSel[order(corDfSel$Estimate, decreasing = T),]
  
  subset_seurat_object <- subset(atacFilt, Annotations == cluster)
  # subset_seurat_object = NormalizeData(subset_seurat_object)
  # subset_seurat_object = ScaleData(subset_seurat_object)
  
  if (topN == "all") {
    genesList = corDfSel$Gene
    curHeat<-DoHeatmap(subset_seurat_object, features = genesList,  size = 12)+
      guides(color="none") +
      theme(legend.text = element_text(size = 18)) +
      theme(legend.title = element_text(size = 18)) +
      theme(axis.text.y = element_text(size = 3.5))
    
  } else {
    genesList = head(corDfSel$Gene, topN)
    curHeat<-DoHeatmap(subset_seurat_object, features = genesList,  size = 12)+
      theme(text = element_text(size = 18)) +
      guides(color="none") +
      theme(legend.text = element_text(size = 18)) +
      theme(axis.text.y = element_text(size = 18)) +
      theme(legend.title = element_text(size = 18))
  }
  
  #print(curHeat)
  return(curHeat)
  
}

targDir = "atacRna/PredictGeneActCor/Heatmaps/"
dir.create(targDir, recursive = T, showWarnings = F)

clusters = unique(atacFilt$Annotations)

topN = 50
for (cluster in clusters) {
  try({
    curHeatMap = makeHeat(atacFilt=atacFilt, cluster=cluster,  topN=topN, corDf=corDf)
    png(filename = paste0(targDir, cluster, "_", topN, "_PredGene_RNA_spearman_2023-10-10.png"), width = 20, height = 16, units = "in", res = 300)
    print(curHeatMap)
    dev.off()
  })
}


topN = "all"
for (cluster in clusters) {
  try({
    curHeatMap = makeHeat(atacFilt=atacFilt, cluster=cluster,  topN=topN, corDf=corDf)
    png(filename = paste0(targDir, cluster, "_", topN, "_PredGene_RNA_spearman_2023-10-10.png"), width = 20, height = 16, units = "in", res = 300)
    print(curHeatMap)
    dev.off()
  })
}
