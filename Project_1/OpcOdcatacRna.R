library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)

curDate<-Sys.Date()

targDir<-'./OPC_ODC/atacRna/'
dir.create(targDir)

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

atacInt$MonoocClust <- rnaFilt$MonocClust

DimPlot(atacInt, group.by = 'MonoocClust')
# find RNA markers

Idents(rnaFilt)<-'MonocClust'

DefaultAssay(rnaFilt) <- 'RNA'

rnaMarkers<-FindMarkers(object=rnaFilt, ident.1= '1', only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

rnaMarkers$Genes<-row.names(rnaMarkers)
outRna<-paste0(targDir, 'rnaPosMark_1vsAll_monocClust_macs2.csv')
write.csv(rnaMarkers, outRna, row.names = F)

rnaMarkOrd <- rnaMarkers[order(rnaMarkers$avg_log2FC, decreasing = T),]
rnaMarkTop <- head(rnaMarkOrd, 15)
genes<-c('Plp1', 'Ptgds', 'Trf', 'Calcrl', 'Cspg4', 'Vcan')
rnaMarkComb <- unique(c(rnaMarkTop$Genes, genes))
#
DefaultAssay(atacInt)<-'Combined_peaks'
annotGenes<-data.frame(Annotation(atacInt))
presentGenes<-rnaMarkComb[(rnaMarkComb%in%annotGenes$gene_name)]
#
atacInt <- RegionStats(atacInt, genome = BSgenome.Mmusculus.UCSC.mm10)

atacInt <- LinkPeaks(
  object = atacInt,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = presentGenes
)

clusters<-unique(atacInt$MonoocClust)
Idents(atacInt)<- atacInt$MonoocClust

for ( i in presentGenes){
  try({
    p1<-CoveragePlot(
      object = atacInt ,
      region = i,
      features = i,
      expression.assay = "RNA",
      idents = clusters,
      extend.upstream = 100000,
      extend.downstream = 100000
    )
    outFile<-paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_ext100k_', i, '_', curDate, '.jpeg')
    ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
  })
}

