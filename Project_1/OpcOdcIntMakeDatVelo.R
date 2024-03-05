#library(monocle3)
library(Seurat)
#library(SeuratWrappers)
library(gridExtra)
library(SeuratData)
library(SeuratDisk)


curDate<-Sys.Date()


source('~/Wang/programs/renameClusters.R')

selCells<-subset(x = RNA.combined.norm, idents = c("OPC", "ODC"))

opcOdc<-readRDS('integ_OPC_ODC')

identical(colnames(selCells), colnames(colnames(opcOdc)))

refRNA<-

selCells$rna_rename_clust<-selCells@meta.data$Annotations