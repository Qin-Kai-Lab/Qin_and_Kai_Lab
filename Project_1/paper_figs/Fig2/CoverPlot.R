library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MAST)
library(ggplot2)
library(qlcMatrix)

source('../programs/renameClusters.R')
curDate<-Sys.Date()

targDir = './Paper_figs/Fig2/CoverPlots/'
dir.create(targDir, recursive = T, showWarnings = F)

atacInt<-readRDS('atacMerged_macs2_2')
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
atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

atacFilt[['RNA']]<-rnaFilt@assays$RNA

DefaultAssay(atacFilt)<-'Combined_peaks'

rnaFilt[['Combined_peaks']] <- atacFilt@assays$Combined_peaks

#
#identical(colnames(atacInt), colnames(rnaFilt))
identical(rownames(atacFilt@meta.data), rownames(rnaFilt@meta.data))
identical(atacFilt$Merged_CellName, rnaFilt$Merged_CellName)
identical(colnames(atacFilt), colnames(rnaFilt))
#
atacFilt$Annotations = rnaFilt$Annotations

rm(annotations, atacInt, RNA.combined.norm, rnaFilt)
gc()

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$Cell_Type)
  markers = character()
  dfPct = df[(df$pct.1 >0.25) | (df$pct.2 >0.25), ]
  for (cluster in clusters) {
    dfSub =  dfPct[(dfPct[['Cell_Type']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$AbsLog, decreasing = T),]
    topMark = head(dfOrd$Genes, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}

topMark = getTopMarkers(df = groupMarkers, topNumb = 3)

annotGenes<-data.frame(Annotation(atacFilt))
presentGenes<-topMark[(topMark%in%annotGenes$gene_name)]


atacFilt$Cluster_Group = paste(atacFilt$Annotations, atacFilt$dataset, sep = "_")
Group = unique(atacFilt$dataset)
clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

Cluster_Group <- lapply(clusters, function(x) {
  lapply(Group, function(y) {
    paste(x, y, sep = "_")
  })
})

Cluster_Group = unlist(Cluster_Group)

atacFilt$Cluster_Group <- factor(atacFilt$Cluster_Group , levels = Cluster_Group)

Idents(atacFilt) = atacFilt$Cluster_Group
atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)

runCovPlot = function(obj, genes, datType){
  atacInt <- LinkPeaks(
    object = obj,
    peak.assay = "Combined_peaks",
    expression.assay = "RNA",
    genes.use = genes
  )
  
  clusters<-unique(obj$Cluster_Group)
  Idents(obj)<- obj$Cluster_Group
  
  for ( i in presentGenes){
    try({
      p1<-CoveragePlot(
        object = obj ,
        region = i,
        features = i,
        expression.assay = "RNA",
        idents =  clusters,
        extend.upstream = 100000,
        extend.downstream = 100000
      )
      outFile<-paste0(targDir, 'CovPlot_', i, '_ext100K_', datType, "_",  curDate, '.png')
      ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
    })
  }
  
}

runCovPlot(obj = atacFilt, genes = presentGenes, datType = "Top3Genes")