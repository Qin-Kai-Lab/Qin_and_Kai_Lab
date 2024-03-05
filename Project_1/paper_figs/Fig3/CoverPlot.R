library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MAST)
library(ggplot2)
library(qlcMatrix)

curDate<-Sys.Date()

targDir = './Paper_figs/Fig3/CoverPlots/'
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

atacInt[['RNA']]<-rnaFilt@assays$RNA

DefaultAssay(atacInt)<-'Combined_peaks'

rnaFilt[['Combined_peaks']] <- atacInt@assays$Combined_peaks

#
#identical(colnames(atacInt), colnames(rnaFilt))
identical(rownames(atacInt@meta.data), rownames(rnaFilt@meta.data))
identical(atacInt$Merged_CellName, rnaFilt$Merged_CellName)
atacInt$MonocClust <- rnaFilt$MonocClust
rm(RNA.combined.norm, rnaFilt, annotations)
#
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

# custom genes
genes<-c('Plp1', 'Ptgds', 'Trf', 'Calcrl', 'Cspg4', 'Vcan')

annotGenes<-data.frame(Annotation(atacInt))
presentGenes<-genes[(genes%in%annotGenes$gene_name)]
#

atacInt$newMonocClust = factor(atacInt$newMonocClust, levels = c("OPC", "Intermideate", "ODC"))

runCovPlot = function(obj, genes, datType){
  atacInt <- LinkPeaks(
    object = obj,
    peak.assay = "Combined_peaks",
    expression.assay = "RNA",
    genes.use = genes
  )
  
  clusters<-unique(obj$newMonocClust)
  Idents(obj)<- obj$newMonocClust
  
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

runCovPlot(obj = atacInt, genes = presentGenes, datType = "CustomGenes")

# top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$cluster)
  markers = character()
  dfPct = df[(df$pct.1 >0.25) | (df$pct.2 >0.25), ]
  for (cluster in clusters) {
    dfSub =  dfPct[(dfPct[['cluster']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$avg_log2FC, decreasing = T),]
    topMark = head(dfOrd$gene, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}

allMarkers = read.csv("Paper_figs/Fig3/All_RNA_Markers_2023-05-08.csv")
topMark = getTopMarkers(df = allMarkers, topNumb = 5)
presentGenes<-topMark[(topMark%in%annotGenes$gene_name)]
runCovPlot(obj = atacInt, genes = presentGenes, datType = "Top5PosGenes")

# pseudotiming genes
allMarkers = read.csv('OPC_ODC/Monocle3/86PC/OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv')

# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$cluster)
  markers = character()
  minNumb = ncol(atacInt) * 0.1
  dfPct = df[(df$num_cells_expressed > minNumb), ]
  dfSort = dfPct[order(dfPct$q_value, decreasing = F),]
  markers = head(dfSort, topNumb)
  unMark = unique(markers$gene_id)
  return(unMark)
}

topMark = getTopMarkers(df = allMarkers, topNumb = 10)
presentGenes<-topMark[(topMark%in%annotGenes$gene_name)]
runCovPlot(obj = atacInt, genes = presentGenes, datType = "Top10TimeGenes")
