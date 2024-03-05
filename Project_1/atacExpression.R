library(Seurat)
library(Signac)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

atacInt<-readRDS('atacIntegrated')

DefaultAssay(atacInt)

#gene.activity<-GeneActivity(atacInt)


atacCells<-rownames(atacInt@meta.data)

source('../../programs/renameClusters.R')

rnaCells<-rownames(RNA.combined.norm@meta.data)

length(atacCells)
length(rnaCells)

length(atacCells[(atacCells%in%rnaCells)])

# need to add peak_region_fragments column
# differential expression
DefaultAssay(atacInt)
da_peaks<-FindMarkers(
  object=atacInt,
  ident.1 = rownames(atacInt[[]][atacInt$dataset=='Stress',]),
  ident.2 = rownames(atacInt[[]][atacInt$dataset=='Control',]),
  min.pct=0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

# test which barcodes match with second metadata and if barcodes are unique among the groups (no)
q1<-atacInt@meta.data[(atacInt@meta.data$dataset=='Control'),]

cellId<-q1$CellName

length(cellId[cellId%in%per_barcode_metrics$gex_barcode])

allCellId<-atacInt$CellName

length(allCellId[allCellId%in%per_barcode_metrics$gex_barcode])
##
ClosestFeature(atacInt, regions = rownames(atacInt))

DefaultAssay(atacInt)
# check correlationm of components
RunUMAP(atacInt, reduction = "integrated_lsi", dims = 2:30)

DepthCor(atacInt, reduction='integrated_lsi')