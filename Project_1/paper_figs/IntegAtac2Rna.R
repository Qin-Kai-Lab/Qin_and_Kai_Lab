library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)
library(rtracklayer)


setwd("~/Wang/output")

source('../programs/renameClusters.R')

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv")

curDate<-Sys.Date()

atacInt<-readRDS('atacIntegrated_macs2_2')
DefaultAssay(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacInt) <- annotations

#combine rna and atac
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)
rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )
atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

identical(colnames(rnaFilt), colnames(atacFilt))
identical(rownames(rnaFilt@meta.data), rownames(atacFilt@meta.data))
identical(rnaFilt$Merged_CellName, atacFilt$Merged_CellName)

atacFilt[['RNA']]<-rnaFilt@assays$RNA
DefaultAssay(atacInt)<-'Combined_peaks'
atacFilt$Annotations = rnaFilt$Annotations

rm(RNA.combined.norm, rnaFilt, atacInt)
gc()

Idents(atacFilt) = atacFilt$Annotations

#
cond<-'Control_Stress'
targetDir<-paste0('Atac_Rna_coverPlots/macs2_Integ/', cond, '/')
dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(atacFilt))

atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)

table(atacFilt$Annotations)

gene.activities <- GeneActivity(atacFilt, assay='Combined_peaks',  max.width = NULL)

saveRDS(gene.activities, "PredGeneActivities_atacIntegrated_macs2_2")

atacFilt[['PredictActivity']] <- CreateAssayObject(counts = gene.activities)
atacFilt <- NormalizeData(
  object = atacFilt,
  assay = 'PredictActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacFilt$nCount_RNA)
)

saveRDS(atacFilt, "atacIntegrated_macs2_2_RNA")

