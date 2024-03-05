library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(clustree)



# get seoarate cell names for OPC | ODC, control and stress

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

source('../../programs/renameClusters.R')

metadat<-RNA.combined.norm@meta.data

metaSel<-metadat[(metadat$Annotations == 'OPC') | (metadat$Annotations== 'ODC'),]

controlNames<-metaSel$CellName[(metaSel$group=='Control')]

stressNames<-metaSel$CellName[(metaSel$group=='Stress')]

write.table(stressNames, 'OPC_ODC_stressCells.txt', row.names = F, col.names = F, quote = F)

write.table(controlNames, 'OPC_ODC_contrCells.txt', row.names = F, col.names = F, quote = F)

rm(list=ls())

controlNames<-read.table('OPC_ODC_contrCells.txt', F)

stressNames<-read.table('OPC_ODC_stressCells.txt', F)

#
controlCells<-controlNames$V1

stressCells<-stressNames$V1

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

# import seurat Control data
control = Read10X(data.dir = "../data/control")
RNA_control = CreateSeuratObject(counts = control$`Gene Expression`)
RNA_control$group = "Control"
RNA_control[["CellName"]] <- colnames(RNA_control)
# filter
controlFiltered <- subset(RNA_control, subset = CellName %in% controlCells )
rm(RNA_control, control)

# import and filter stress group
stress = Read10X(data.dir = "../data/stress")
RNA_stress = CreateSeuratObject(counts = stress$`Gene Expression`)
RNA_stress$group = "Stress"
RNA_stress[["CellName"]] <- colnames(RNA_stress)
# filter
stressFiltered <- subset(RNA_stress, subset = CellName %in% stressCells)
rm(RNA_stress, stress)

#you can use RNA_combined.list as "ifnb.list" in their example
RNA_combined.list = list(control = controlFiltered, stress = stressFiltered)

# normalize and identify variable features for each dataset independently
RNA_combined.list <- lapply(X = RNA_combined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = RNA_combined.list)
immune.anchors <- FindIntegrationAnchors(object.list = RNA_combined.list, anchor.features = features)
RNA.combined.norm <- IntegrateData(anchorset = immune.anchors)
rm(list = ls()[!ls() %in% c("RNA.combined.norm")])

DefaultAssay(RNA.combined.norm)<-"integrated"
RNA.combined.norm  <- ScaleData(RNA.combined.norm , verbose = FALSE)
RNA.combined.norm <- RunPCA(RNA.combined.norm, verbose = FALSE)

###
ElbowPlot(RNA.combined.norm, ndims = 50)
#ggsave("ElbowPlot_DupFiltr.jpeg", width = 16, height = 10)

# function borrowed from Sergei Bombin
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

minPCA<-fgetNPCs(obj=RNA.combined.norm, MIN_CUM_SD=95)

RNA.combined.norm <- JackStraw(RNA.combined.norm, num.replicate = 100, dims = 50)
RNA.combined.norm <- ScoreJackStraw(RNA.combined.norm, dims = 1:50)
JackStrawPlot(RNA.combined.norm, dims = 1:50)
ggsave("JackStrawPlot_OPC_ODC.jpeg", width = 16, height = 10)

###
# clustering
RNA.combined.norm <- FindNeighbors(RNA.combined.norm, dims = 1:37)
RNA.combined.norm <- RunUMAP(RNA.combined.norm, reduction = "pca", dims = 1:37)
resols<-c(0,0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
RNA.combined.norm <- FindClusters(RNA.combined.norm, resolution = resols)


setwd('../output')

# terst clusters
clustree <- FindClusters(RNA.combined.norm , resolution = resols)
p <- clustree(clustree)
ggsave("Clustree.jpeg", plot=p, width = 8, height = 10)

curDate<-Sys.Date()

plotUmap<-function(x){
  for (i in x){
    current_resol<-paste0("integrated_snn_res.", i)
    Idents(RNA.combined.norm)<-current_resol
    p1 <- DimPlot(RNA.combined.norm, reduction = "umap", group.by = "group", label = F)+ggtitle(current_resol)
    p2 <- DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
    p3<-grid.arrange(p1,p2, ncol=2, nrow=1)
    ggsave(paste0('OPC_ODC_umapClustering_',current_resol, "_", curDate, ".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)
  }
}

plotUmap(resols)





saveRDS(RNA.combined.norm, file = 'integ_OPC_ODC')