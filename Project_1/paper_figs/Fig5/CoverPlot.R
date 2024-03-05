library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MAST)
library(ggplot2)
library(qlcMatrix)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig5/CoverPlot/'
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

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

annotGenes<-data.frame(Annotation(atacInt))
presentGenes<-topMark[(topMark%in%annotGenes$gene_name)]

atacInt$Group_Cluster = paste(atacInt$dataset, atacInt$newMonocClust , sep = "_")
atacInt$Group_Cluster <- factor(atacInt$Group_Cluster, levels = c("Control_OPC", "Stress_OPC", "Control_Intermideate", "Stress_Intermideate", "Control_ODC", "Stress_ODC"))

# plots
runCovPlot = function(obj, genes, datType){
  atacInt <- LinkPeaks(
    object = obj,
    peak.assay = "Combined_peaks",
    expression.assay = "RNA",
    genes.use = genes
  )
  
  clusters<-unique(obj$Group_Cluster)
  Idents(obj)<- obj$Group_Cluster
  
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

runCovPlot(obj = atacInt, genes = presentGenes, datType = "Top5Genes")