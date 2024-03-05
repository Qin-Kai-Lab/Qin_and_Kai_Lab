library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)

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

covPlotPerClust = function(seurObj, curCluster, genesDf, annotGenes, clusters, targetDir, curDate) {
  posGenes = genesDf[genesDf$avg_log2FC > 0,]
  posGenesPct = posGenes[(posGenes$pct.1 >0.25) | (posGenes$pct.2 >0.25), ]
  posGenesCl = posGenesPct[posGenesPct$cluster == curCluster,]
  posGenesClOrd = posGenesCl[order(posGenesCl$avg_log2FC, decreasing = T),]
  top10 = head(posGenesClOrd$gene, 12)
  presentGenes<- top10[(top10%in%annotGenes$gene_name)]
  
  #seurObj = subset(seurObj, subset = Annotations %in% clusters )
  seurObj<- LinkPeaks(
    object = seurObj,
    peak.assay = "Combined_peaks",
    expression.assay = "RNA",
    genes.use = presentGenes
  )
  
  for ( i in presentGenes){
    p1<-CoveragePlot(
      object = seurObj ,
      region = i,
      features = i,
      expression.assay = "RNA",
      idents = clusters,
      extend.upstream = 20000,
      extend.downstream = 20000
    ) & theme(text = element_text(size = 22)) & scale_fill_grey()
    
  
    outFile<-paste0(targetDir,'/CovPlot_AtacRna_',curCluster, '_', i, '_', curDate, '.png')
    ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
  }
  
}


unique(atacFilt$Annotations)
# excluded C-R cluster as it has few cells
clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC')

for (curCluster in clusters) {
  targetDir<-paste0('Atac_Rna_coverPlots/macs2_Integ/', cond, '/Grey_Col/', curCluster, "/")
  dir.create(targDir, recursive = T, showWarnings = F)
  covPlotPerClust(seurObj=atacFilt, curCluster=curCluster, genesDf=allMarkers, annotGenes=annotGenes, clusters=clusters, targetDir=targetDir, curDate=curDate)
}
