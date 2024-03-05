library(MAST)
library(EnsDb.Mmusculus.v79)
source('../programs/renameClusters.R')

curDate<-Sys.Date()


### find all Markers for CA1 separately for stress and control
cluster<-c('CA1')

DefaultAssay(RNA.combined.norm)

conditions<-unique(RNA.combined.norm$group)


targDir<-'RNA_cluster_markers/'
dir.create(targDir)
  
for ( i in conditions){
  df<-subset(x = RNA.combined.norm, subset = group == i)
  allMarkers<-FindMarkers(object=df, ident.1= cluster, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
  allMarkers$Gene<-rownames(allMarkers)
  outFile<-paste0(targDir, 'allPosMarkers_CA1_', i, '.csv')
  
  write.csv(allMarkers, outFile, row.names = F)
}

contrMark<-read.csv('RNA_cluster_markers/allPosMarkers_CA1_Control.csv')
contrOrd<-contrMark[order(contrMark$avg_log2FC, decreasing = T),]
contrGenes<-head(contrOrd$Gene, 10)

stressMark<-read.csv('RNA_cluster_markers/allPosMarkers_CA1_Stress.csv')
stressOrd<-stressMark[order(stressMark$avg_log2FC, decreasing = T),]
stressGenes<-head(stressOrd$Gene, 10)

identical(contrGenes, stressGenes)

# build atacseq coverage plots 

atacInt<-readRDS('atacMerged_macs2')

DefaultAssay(atacInt)
Annotation(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(atacInt) <- annotations

atacContr<-subset(x = atacInt, subset = dataset == 'Stress')

CoveragePlot(
  object = atacContr,
  region = c("Cntnap2"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)

DimPlot(atacContr, label=T)

CoveragePlot(
  object = atacInt,
  region = c("Cntnap2"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)