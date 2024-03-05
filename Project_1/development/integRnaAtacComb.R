library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)
source('../programs/renameClusters.R')

curDate<-Sys.Date()


atacInt<-readRDS('atacIntegrated_macs2')
DefaultAssay(atacInt)

DefaultAssay(atacInt)
#Annotation(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(atacInt) <- annotations
# add RNA data to ATAC
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)

rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )

atacInt[['RNA']]<-rnaFilt@assays$RNA

DefaultAssay(atacInt)<-'Combined_peaks'

### per condition

cond<-'Control'

targetDir<-paste0('Atac_Rna_coverPlots/macs2_Integ/', cond, '/')


dir.create(targetDir, recursive = T)

atacSub<-subset(x = atacInt, subset = dataset == cond)

annotGenes<-data.frame(Annotation(atacSub))

newMarkers<-read.csv(paste0('RNA_cluster_markers/allPosMarkers_CA1_', cond, '.csv'))
newOrd<-newMarkers[order(newMarkers$avg_log2FC, decreasing = T),]
newGenes<-head(newOrd$Gene, 12)
presentGenes<-newGenes[(newGenes%in%annotGenes$gene_name)]

atacSub <- RegionStats(atacSub, genome = BSgenome.Mmusculus.UCSC.mm10)

atacSub <- LinkPeaks(
  object = atacSub,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = newGenes
)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

for ( i in presentGenes){
  p1<-CoveragePlot(
    object = atacSub ,
    region = i,
    features = i,
    expression.assay = "RNA",
    idents = clusters,
    extend.upstream = 10000,
    extend.downstream = 10000
  )
  outFile<-paste0(targetDir, 'CovPlot_AtacRna_integr_macs2_ext10k_', i, '_', curDate, '_', cond, '.jpeg')
  ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
}
