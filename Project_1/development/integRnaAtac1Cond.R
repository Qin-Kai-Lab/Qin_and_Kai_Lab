library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
curDate<-Sys.Date()
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)

setwd("/home/flyhunter/Wang/output")

stressPeak <- readRDS('../data/atacAllFiltMacs2_control_1.87gensize_control')

atacInt<-readRDS('atacMerged_macs2')

atacStress<-subset(x = atacInt, subset = dataset == 'Control')

stressSub<-subset(stressPeak, subset = gex_barcode %in% atacStress$gex_barcode )

identical(colnames(stressSub), colnames(atacStress))

atacStress@meta.data$Cluster<-atacStress@meta.data$Annotations
ind = match(atacStress$gex_barcode, stressSub$gex_barcode)
#stressSub$Cluster<-NA
stressSub@meta.data[ind, "Cluster"] = atacStress@meta.data[, "Cluster"]
stressSub@meta.data[ind, "Annotations"] = atacStress@meta.data[, "Annotations"]

Idents(stressSub)<-stressSub$Annotations

stressSub<- FindTopFeatures(stressSub, min.cutoff = 10)
stressSub<- RunTFIDF(stressSub)
stressSub <- RunSVD(stressSub)
stressSub <- RunUMAP(stressSub, reduction = "lsi", dims = 2:30)

DimPlot(stressSub)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(stressSub) <- annotations


stressMark<-read.csv('RNA_cluster_markers/allPosMarkers_CA1_Control.csv')
stressOrd<-stressMark[order(stressMark$avg_log2FC, decreasing = T),]
stressGenes<-head(stressOrd$Gene, 10)

stressGenes

CoveragePlot(
  object = stressSub,
  region = c("Galntl6"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)

#allMarkers <- FindMarkers(stressSub , only.pos = T, min.pct = 0.05, test.use = 'LR', 
#                          latent.vars = 'atac_peak_region_fragments', 
#                          ident.1 = "CA1")


#allMarkers$query_region<-row.names(allMarkers)

#closGenes<-ClosestFeature(atacInt, regions = rownames(allMarkers))

#allMarkEdit<-plyr::join(allMarkers, closGenes, by='query_region', type='left', match='all')

#write.csv(allMarkEdit, 'atacControl_macs2_markersCA1.csv', row.names = F)

stressSub <- RegionStats(stressSub, genome = BSgenome.Mmusculus.UCSC.mm10)



# get RNA data
# import seurat Control data
control = Read10X(data.dir = "../data/control/outs")
RNA_control = CreateSeuratObject(counts = control$`Gene Expression`)
RNA_control$group = "Control"
RNA_control[["CellName"]] <- colnames(RNA_control)
# filter
controlFiltered <- subset(RNA_control, subset = CellName %in% stressSub$gex_barcode )
controlFiltered<-NormalizeData(controlFiltered)
controlFiltered <- ScaleData(controlFiltered , verbose = FALSE)

stressSub[['RNA']]<-controlFiltered@assays$RNA

rm(atacInt, atacStress, controlFiltered, stressPeak)
gc()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

# trey to find other specific set of Markers
DefaultAssay(stressSub)<-'RNA'

newMarkers<-FindMarkers(object=stressSub, ident.1= 'CA1', ident.2 = NULL, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
newMarkers$Genes<-rownames(newMarkers)
newOrd<-newMarkers[order(newMarkers$avg_log2FC, decreasing = T),]
newGenes<-head(newOrd$Genes, 12)

newGenes
# macs2 1 condition genes
DefaultAssay(stressSub)<-'macs2'

Annotation(stressSub) <- annotations

stressSub <- RegionStats(stressSub, genome = BSgenome.Mmusculus.UCSC.mm10)

stressSub <- LinkPeaks(
  object = stressSub,
  peak.assay = "macs2",
  expression.assay = "RNA",
  genes.use = newGenes
)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

targetDir<-'Atac_Rna_coverPlots/Control_only/macs2/ext50kb/'

dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(stressSub))

newGenes<-c(newGenes, 'Meis2', 'Mpped1', 'Satb2')

presentGenes<-newGenes[(newGenes%in%annotGenes$gene_name)]

for ( i in presentGenes){
  p1<-CoveragePlot(
    object = stressSub ,
    region = i,
    features = i,
    expression.assay = "RNA",
    idents = clusters,
    extend.upstream = 50000,
    extend.downstream = 50000
  )
  outFile<-paste0(targetDir, 'CovPlot_AtacRna_macs2_contrOnly_ext50k_', i, '_', curDate, '.jpeg')
  ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
}


# try with original peaks
DefaultAssay(stressSub)<-'peaks'
stressSub<- FindTopFeatures(stressSub, min.cutoff = 10)
stressSub<- RunTFIDF(stressSub)
stressSub <- RunSVD(stressSub)
stressSub <- RunUMAP(stressSub, reduction = "lsi", dims = 2:30)

Annotation(stressSub) <- annotations

stressSub <- RegionStats(stressSub, genome = BSgenome.Mmusculus.UCSC.mm10)

stressSub <- LinkPeaks(
  object = stressSub,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = newGenes
)


targetDir<-'Atac_Rna_coverPlots/Control_only/original_peaks/ext50kb/'

dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(stressSub))

presentGenes<-newGenes[(newGenes%in%annotGenes$gene_name)]

for ( i in presentGenes){
  p1<-CoveragePlot(
    object = stressSub ,
    region = i,
    features = i,
    expression.assay = "RNA",
    idents = clusters,
    extend.upstream = 50000,
    extend.downstream = 50000
  )
  outFile<-paste0(targetDir, 'CovPlot_AtacRna_origPeak_contrOnly_ext50k_', i, '_', curDate, '.jpeg')
  ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
}