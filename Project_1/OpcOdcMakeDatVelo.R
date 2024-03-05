#library(monocle3)
library(Seurat)
#library(SeuratWrappers)
library(gridExtra)
library(SeuratData)
library(SeuratDisk)


curDate<-Sys.Date()

source('~/Wang/programs/renameClusters.R')

selCells<-subset(x = RNA.combined.norm, idents = c("OPC", "ODC"))

selCells$rna_rename_clust<-selCells@meta.data$Annotations

DefaultAssay(selCells)<-'integrated'

rm(RNA.combined.norm)

RNA.combined.norm<-selCells

RNA.combined.norm <- RunPCA(RNA.combined.norm, verbose = FALSE)

#RNA.combined.norm <- JackStraw(RNA.combined.norm, num.replicate = 100, dims = 30)
#RNA.combined.norm <- ScoreJackStraw(RNA.combined.norm, dims = 1:30)
#JackStrawPlot(RNA.combined.norm, dims = 1:30)

RNA.combined.norm <- FindNeighbors(RNA.combined.norm, dims = 1:15)
RNA.combined.norm <- RunUMAP(RNA.combined.norm, reduction = "pca", dims = 1:15)
resols<-c(0,0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
RNA.combined.norm <- FindClusters(RNA.combined.norm, resolution = resols)

i<-0.2
current_resol<-paste0("integrated_snn_res.", i)
Idents(RNA.combined.norm)<-current_resol
p1 <- DimPlot(RNA.combined.norm, reduction = "umap", group.by = "group", label = F)+ggtitle(current_resol)
p2 <- DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
p3<-grid.arrange(p1,p2, ncol=2, nrow=1)

#ggsave(paste0('OPC_ODC_umapClustering_OrigInt_',current_resol, "_", curDate, ".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)


DefaultAssay(RNA.combined.norm)<-'RNA'

RNA.combined.norm@meta.data$scvelo_name<-NA

for (i in 1:nrow(RNA.combined.norm@meta.data)){
  if (RNA.combined.norm@meta.data$group[i]=='Control'){
    RNA.combined.norm@meta.data$scvelo_name[i]<-paste0('control:', gsub('-1', 'x', RNA.combined.norm@meta.data$CellName[i]))
  } else if (RNA.combined.norm@meta.data$group[i]=='Stress'){
    RNA.combined.norm@meta.data$scvelo_name[i]<-paste0('stress:', gsub('-1', 'x', RNA.combined.norm@meta.data$CellName[i]))
  }
}

DefaultAssay(RNA.combined.norm)

RNA.combined.norm<- RenameCells(RNA.combined.norm, new.names = RNA.combined.norm@meta.data$scvelo_name)

DimPlot(RNA.combined.norm, reduction = 'umap')


SaveH5Seurat(RNA.combined.norm, filename = "opc_odc.h5Seurat")
Convert("opc_odc.h5Seurat", dest = "h5ad")