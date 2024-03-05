library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(future)

setwd("/home/flyhunter/Wang/data")

stress <- readRDS('atacStress_2_macs')
control <- readRDS('atacControl_2_macs')

DefaultAssay(stress)<-'macs2'
DefaultAssay(control)<-'macs2'

# find common peaks
stress_peak<-granges(stress)

control_peak<-granges(control)

identical(stress_peak, control_peak)

combined.peaks <- reduce(x = c(control_peak, stress_peak))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

control_frag<-Fragments(control)

stress_frag<-Fragments(stress)

identical(control_frag, stress_frag)

control_count<- FeatureMatrix(fragments = control_frag, features = combined.peaks)
#
stress_count<- FeatureMatrix(fragments = stress_frag, features = combined.peaks)

identical(control_count, stress_count)

###
###
control_chrom<-CreateChromatinAssay(control_count, fragments = control_frag)

stress_chrom<-CreateChromatinAssay(stress_count, fragments = stress_frag)

control[['Combined_peaks']]<-control_chrom
DefaultAssay(control)<-'Combined_peaks'

stress[['Combined_peaks']]<-stress_chrom
DefaultAssay(stress)<-'Combined_peaks'

control$dataset <- "Control"
stress$dataset <- "Stress"

rm(list = ls()[!ls() %in% c("control", "stress")])

gc()

control <- FindTopFeatures(control, min.cutoff = 10)
control <- RunTFIDF(control)
control <- RunSVD(control)

stress <- FindTopFeatures(stress, min.cutoff = 10)
stress <- RunTFIDF(stress)
stress <- RunSVD(stress)

#saveRDS(stress, file = 'atacAllFiltMacs2_stress')

#saveRDS(control, file = 'atacAllFiltMacs2_control')

# start integrations


stress[['macs2']]<-NULL

control[['macs2']]<-NULL


DefaultAssay(stress)
DefaultAssay(control)

pbmc.combined <- merge(control, stress)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")
p1


##

# rownames control and stress are identical

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(control, stress),
  anchor.features = rownames(stress),
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

setwd('../output')


saveRDS(integrated, file = 'atacIntegrated_macs2_2')
saveRDS(pbmc.combined, file = 'atacMerged_macs2_2')