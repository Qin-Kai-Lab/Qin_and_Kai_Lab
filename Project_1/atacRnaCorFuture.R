library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(foreach)
library(doParallel)
library(future)
plan("multicore", workers = 5)

options(future.globals.maxSize = 58 * 1024 ^ 3)

options(future.globals.seed=TRUE)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

source('../programs/renameClusters.R')

identical(RNA.combined.norm@assays$RNA@counts, RNA.combined.norm@assays$RNA@data)

atacInt = readRDS("atacIntegrated_macs2_2")

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacInt) <- annotations

# combine
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)

rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )

atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

atacFilt[['RNA']]<-rnaFilt@assays$RNA

rm(RNA.combined.norm, atacInt)
gc()

atacFilt = RegionStats(object = atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)

identical(rownames(atacFilt@meta.data), rownames(rnaFilt@meta.data))
atacFilt$Annotations = rnaFilt$Annotations


rm(rnaFilt)
gc()

clusters = unique(atacFilt$Annotations)
clusters = c("CA3", "CA2", "GABA", "SUB", "C-R")

getPeakRna = function(dat, test, targetDir) {
  dir.create(targetDir, recursive = T)
  for (cluster in clusters) {
    try({
      atacSub = subset(x = dat, subset = Annotations == cluster)
      peakRna = LinkPeaks(
        object = atacSub,
        peak.assay = 'Combined_peaks',
        expression.assay = "RNA",
        peak.slot = "counts",
        expression.slot = "data",
        method = test,
        distance = 5e+05,
        min.cells = 10,
        pvalue_cutoff = 0.05,
        score_cutoff = 0.05,
        verbose = TRUE)
      p2g = data.frame(Links(peakRna))
      p2g$Cluster = cluster
      write.csv(p2g, file = paste0(targetDir, "peakGeneCor_", cluster, "_", curDate, ".csv"), row.names = F)
      print(paste0(cluster, "_____Done!"))
    })
    
  }
}

print("Starting finding correlations")

getPeakRna(dat = atacFilt, test = "pearson", targetDir = "atacRna/peaksGenesCor/atacIntegMacs2_2/pearson/")