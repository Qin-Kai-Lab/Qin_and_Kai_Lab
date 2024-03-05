library(Seurat)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig3/StressVsContr_DimAnalis/'

dir.create(targDir, recursive = T, showWarnings = F)

RNA.combined.norm<-readRDS('integ_OPC_ODC')

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"

clusters = unique(RNA.combined.norm$newMonocClust)

fgetNPCs <- function(obj,MIN_CUM_SD){
  sds <- Stdev(obj,reduction="pca")
  cumsds <- 0
  for(i in 1:length(sds)){
    cumsds <- cumsds + sds[i]
    if(cumsds >= MIN_CUM_SD){
      return(i)
    }
  }
}

dimAnalis = function(x, cluster, dataSet, dimNum) {
  cellType = subset(x = x, subset = newMonocClust == cluster)
  if (dimNum == "Rescale") {
    cellType = NormalizeData(cellType)
    cellType = FindVariableFeatures(cellType)
    cellType = ScaleData(cellType)
    cellType = RunPCA(cellType, verbose = FALSE)
    minPC = fgetNPCs(obj=cellType, MIN_CUM_SD=95)
    pcaPlot = DimPlot(cellType, reduction = 'pca', group.by = 'group') + 
      theme(text = element_text(size = 22)) +  
      ggtitle(cluster) + 
      theme(plot.title = element_text(size = 25))
    
    cellType <- FindNeighbors(cellType, dims = 1:minPC)
    cellType <- RunUMAP(cellType, reduction = "pca", dims = 1:minPC)
    uPlot = DimPlot(cellType, reduction = "umap", group.by = "group", label = F) + 
      theme(text = element_text(size = 22)) + 
      ggtitle(cluster) + 
      theme(plot.title = element_text(size = 25))


  } else {
    pcaPlot = DimPlot(cellType, reduction = 'pca', group.by = 'group') + 
      theme(text = element_text(size = 22)) +  ggtitle(cluster) + 
      theme(plot.title = element_text(size = 25))
    uPlot = DimPlot(cellType, reduction = "umap", group.by = "group", label = F) + 
      theme(text = element_text(size = 22)) + 
      ggtitle(cluster) + 
      theme(plot.title = element_text(size = 25))
  }
  
  pcaPlot
  uPlot
  
  #png(paste0(targDir, 'PCA_', dataSet, "_", cluster, "_", dimNum, "_", curDate, '.png'), width = 16, height = 12, units = 'in', res = 300)
  #print(pcaPlot)
  #dev.off()
  
  ggsave(plot = pcaPlot, file = paste0(targDir, 'PCA_', dataSet, "_", cluster, "_", dimNum, "_", curDate, '.png'), width = 16, height = 12, dpi = 300, units = 'in')
  
  #png(paste0(targDir, 'UMAP_', dataSet, "_", cluster, "_", dimNum, "_", curDate, '.png'), width = 16, height = 12, units = 'in', res = 300)
  #print(uPlot)
  #dev.off()
  
  ggsave(plot = uPlot, file = paste0(targDir, 'UMAP_', dataSet, "_", cluster, "_", dimNum, "_", curDate, '.png'), width = 16, height = 12, dpi = 300, units = 'in')
  
}

for (cluster in clusters) {
  dimAnalis(x =RNA.combined.norm , cluster = cluster , dataSet = "RNA", dimNum = "Rescale")
}

for (cluster in clusters) {
  dimAnalis(x =RNA.combined.norm , cluster = cluster , dataSet = "RNA", dimNum = "Default")
}

