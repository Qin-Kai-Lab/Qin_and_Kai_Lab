library(Seurat)
library(Signac)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

stress <- readRDS('atacRnaStress')
control <- readRDS('atacRnaControl')

identical(control, stress)

DefaultAssay(control) <- 'ATAC'

DefaultAssay(stress) <- 'ATAC'


contrGenes<-control@assays$ATAC@data@Dimnames
contrGenes<-contrGenes[[1]]

stressGenes<-stress@assays$ATAC@data@Dimnames
stressGenes<-stressGenes[[1]]

length(contrGenes)
length(stressGenes)

length(stressGenes[(stressGenes%in%contrGenes)])

# stress cell IDs 
stressCellId<-rownames(stress@meta.data)

# try to get overlapping regions
fragpath <- "./stress/atac_fragments.tsv.gz"
fragcounts <- CountFragments(fragments = fragpath)

# checks
atac.cells <- fragcounts[fragcounts$frequency_count > 2000, "CB"]
atac.cells<-subset(fragcounts, subset = CellName %in% stressCellId)

length(stressCellId)
length(stressCellId[stressCellId%in%atac.cells])
length(atac.cells[(atac.cells%in%stressCellId)])
#
atac.frags <- CreateFragmentObject(path = fragpath, cells = stressCellId)
length(atac.frags@cells)

counts <- FeatureMatrix(
  fragments = atac.frags,
  features = granges(control),
  cells = stressCellId
)


atac.assay <- CreateChromatinAssay(
  counts = counts,
  fragments = atac.frags,
  sep=c(":", "-"),
  min.cells=1,
  genome = "mm10"
)

stress[['ATAC']]<-atac.assay
DefaultAssay(stress) <- 'ATAC'

contrGenes<-control@assays$ATAC@data@Dimnames
contrGenes<-contrGenes[[1]]

stressGenes<-stress@assays$ATAC@data@Dimnames
stressGenes<-stressGenes[[1]]

length(contrGenes)
length(stressGenes)

length(stressGenes[(stressGenes%in%contrGenes)])

#stress[['ATAC2']] <- NULL
# assays have the same ranges
DefaultAssay(control)
DefaultAssay(stress)

control$dataset <- "Control"
stress$dataset <- "Stress"


control <- FindTopFeatures(control, min.cutoff = 10)
control <- RunTFIDF(control)
control <- RunSVD(control)

stress <- FindTopFeatures(stress, min.cutoff = 10)
stress <- RunTFIDF(stress)
stress <- RunSVD(stress)

#
pbmc.combined <- merge(control, stress)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")
p1

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(control, stress),
  anchor.features = rownames(control),
  reduction = "rlsi",
  dims = 2:30
)


# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)


integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "dataset")

(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))

p3<-(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))

ggsave('../output/atac_integrated_allDim.jpeg', plot = p3, height = 6, width = 10, units = 'in', dpi = 300)


saveRDS(integrated, file = 'atacIntegrated')