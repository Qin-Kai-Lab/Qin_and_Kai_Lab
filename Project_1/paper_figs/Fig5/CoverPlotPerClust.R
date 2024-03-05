library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MAST)
library(ggplot2)
library(qlcMatrix)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig5/CoverPlot/PerClust/'
dir.create(targDir, recursive = T, showWarnings = F)

atacInt<-readRDS('../data/OpcOdc_atacMerged_macs2')
DefaultAssay(atacInt)

DefaultAssay(atacInt)
#Annotation(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacInt) <- annotations

# add RNA data to ATAC
RNA.combined.norm<-readRDS('integ_OPC_ODC')
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)


rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )
atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

atacFilt[['RNA']]<-rnaFilt@assays$RNA
atacInt[['RNA']]<-rnaFilt@assays$RNA
DefaultAssay(atacFilt)<-'Combined_peaks'
rnaFilt[['Combined_peaks']] <- atacFilt@assays$Combined_peaks

identical(rownames(atacInt@meta.data), rownames(rnaFilt@meta.data))
identical(atacInt$Merged_CellName, rnaFilt$Merged_CellName)
identical(colnames(atacFilt), colnames(rnaFilt))
atacInt$MonocClust <- rnaFilt$MonocClust
rm(RNA.combined.norm, rnaFilt, annotations, atacFilt)
gc()

renameClusters = function(df, clustCol) {
  newClust = character()
  for ( i in 1:nrow(df)) {
    if (df[i, clustCol] == 1) {
      curClust = "ODC"
    } else if (df[i, clustCol] == 2) {
      curClust = "OPC"
    } else if ((df[i, clustCol] == 3) | (df[i, clustCol] == 4)) {
      curClust = "Intermideate"
    }
    newClust = c(newClust, curClust)
  }
  df$newMonocClust = newClust
  return(df)
}

atacInt@meta.data = renameClusters(df = atacInt@meta.data, clustCol = "MonocClust")
atacInt <- RegionStats(atacInt, genome = BSgenome.Mmusculus.UCSC.mm10)

# select genes
groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
groupMarkers = renameClusters(df = groupMarkers, clustCol = "Cell_Type")

# plots
runCovPlot = function(obj, genes, datType){
  atacInt <- LinkPeaks(
    object = obj,
    peak.assay = "Combined_peaks",
    expression.assay = "RNA",
    genes.use = genes
  )
  
  clusters<-unique(obj$dataset)
  Idents(obj)<- obj$dataset
  
  for ( i in genes){
    try({
      p1<-CoveragePlot(
        object = obj ,
        region = i,
        features = i,
        expression.assay = "RNA",
        idents = clusters,
        extend.upstream = 100000,
        extend.downstream = 100000
      )
      outFile<-paste0(targDir, 'CovPlot_', i, '_ext100K_', datType, "_",  curDate, '.png')
      ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
    })
  }
  
}



covPlotClust = function(obj, markDf, cluster, topNumb) {
  objFilt = subset(obj, subset = newMonocClust == cluster)
  markDfFilt = markDf[(markDf$newMonocClust == cluster),]
  markDfFilt$AbsLog = abs(markDfFilt$avg_log2FC)
  dfOrd = markDfFilt[order(markDfFilt$AbsLog, decreasing = T),]
  topMark = head(dfOrd$Genes, topNumb)
  
  runCovPlot(obj = objFilt, genes = topMark, datType = cluster)
}

allClust = unique(atacInt$newMonocClust)
for (curClust in allClust) {
  covPlotClust(obj = atacInt, markDf = groupMarkers, cluster = curClust, topNumb = 5)
}
