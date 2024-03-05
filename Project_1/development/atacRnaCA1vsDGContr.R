library(Seurat)
library(Signac)
library(MAST)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)


curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir<-'./atacRnaContr/macs2/'

dir.create(targDir, recursive = T)

atacInt<-readRDS('atacAllFiltMacs2_control_1.87')

DefaultAssay(atacInt)<-'macs2'



atacMarkers <- FindMarkers(atacInt , only.pos = T, min.pct = 0.05, test.use = 'LR', 
                          latent.vars = 'atac_peak_region_fragments', logfc.threshold = 0.25,
                          ident.1 = "CA1", ident.2 = "DG")

closGenes<-ClosestFeature(atacInt, regions = rownames(atacMarkers))
allMarkEdit<-cbind(atacMarkers, closGenes)

outAtac<-paste0(targDir, 'atacPosMarkCA1vsDG_macs2_1.87_control.csv')
write.csv(allMarkEdit, outAtac, row.names = F)


# RNA markers

DefaultAssay(atacInt)<-'RNA'

rnaMarkers<-FindMarkers(object=atacInt, ident.1= 'CA1', ident.2 = 'DG', only.pos = T, min.pct = 0.05, logfc.threshold = 0.25, test.use = "MAST")

rnaMarkers$Genes<-row.names(rnaMarkers)
outRna<-paste0(targDir, 'rnaPosMarkCA1vsDG_macs2_1.87_control.csv')
write.csv(rnaMarkers, outRna, row.names = F)

# compare lists of genes 
rnaMarkList<-unique(rnaMarkers$Genes)
atacMarkList<-unique(allMarkEdit$gene_name)
rnaInAtac<-rnaMarkList[(rnaMarkList%in%atacMarkList)]
length(rnaInAtac) # 475 out of 525

# infer gene expression from peaks
DefaultAssay(atacInt)<-'macs2'

gene.activities <- GeneActivity(atacInt)
atacInt[['inferAct']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(atacInt)<-'inferAct'
atacInt <- NormalizeData(object = atacInt, assay = 'inferAct')
atacInt <- ScaleData(object = atacInt, assay = 'inferAct')
# find markers
inferMarkers<-FindMarkers(object=atacInt, ident.1= 'CA1', ident.2 = 'DG', only.pos = T, min.pct = 0.05, logfc.threshold = 0.25, test.use = "MAST")
inferMarkers$Genes<-row.names(inferMarkers)
outRna<-paste0(targDir, 'inferRnaPosMarkCA1vsDG_macs2_1.87_control.csv')
write.csv(rnaMarkers, outRna, row.names = F)
# check infered markers vs RNA and ATAC
inferMarkList<-inferMarkers$Genes
# vs RNA
infInRNA<-inferMarkList[(inferMarkList%in%rnaMarkList)]
length(infInRNA) # 164 out of 369
# vs ATAC
infInAtac<-inferMarkList[(inferMarkList%in%atacMarkList)]
length(infInAtac) # 363 out of 369

DefaultAssay(atacInt)<-'macs2'
saveRDS(atacInt, file = 'atacAllFiltMacs2_control_1.87')

# make cover plots based on all marker types

#atacMarkers
DefaultAssay(atacInt)<-'macs2'

#Annotation(stressSub) <- annotations

atacInt <- RegionStats(atacInt, genome = BSgenome.Mmusculus.UCSC.mm10)

atacMarkFilt<-allMarkEdit[(allMarkEdit$distance < 60000) & (allMarkEdit$p_val_adj < 0.05),]
newOrd<-atacMarkFilt[order(atacMarkFilt$avg_log2FC, decreasing = T),] # change this line for RNA and infer RNA
newGenes<-head(newOrd$gene_name, 12) # change this line for RNA and infer RNA
newGenes<-unique(c(newGenes, 'Meis2', 'Mpped1', 'Satb2'))

atacInt <- LinkPeaks(
  object = atacInt,
  peak.assay = "macs2",
  expression.assay = "RNA",
  genes.use = newGenes
)

clusters<-c('CA1','DG')

targetDir<-'Atac_Rna_coverPlots/Control_only/macs2_1.87/atacMarkers/'

dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(atacInt))

presentGenes<-newGenes[(newGenes%in%annotGenes$gene_name)]

for ( i in presentGenes){
  try({
  p1<-CoveragePlot(
    object = atacInt ,
    region = i,
    features = i,
    expression.assay = "RNA",
    idents = clusters,
    extend.upstream = 60000,
    extend.downstream = 60000
  )
  outFile<-paste0(targetDir, 'CovPlot_AtacRna_CA1_DG_macs2_1.87_contrOnly_ext50k_', i, '_', curDate, '.jpeg')
  ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
  })
}

### marker based on infered genes activity
atacInt<-readRDS('atacAllFiltMacs2_control_1.87')
inferMarkers<-read.csv('./atacRnaContr/macs2/inferRnaPosMarkCA1vsDG_macs2_1.87_control.csv')

DefaultAssay(atacInt)<-'macs2'
newOrd<-inferMarkers[order(inferMarkers$avg_log2FC, decreasing = T),] # change this line for RNA and infer RNA
newGenes<-head(newOrd$Genes, 12) # change this line for RNA and infer RNA
newGenes<-unique(c(newGenes, 'Meis2', 'Mpped1', 'Satb2'))

atacInt <- LinkPeaks(
  object = atacInt,
  peak.assay = "macs2",
  expression.assay = "RNA",
  genes.use = newGenes
)

clusters<-c('CA1','DG')

targetDir<-'Atac_Rna_coverPlots/Control_only/macs2_1.87/inferRnaMarkers_CA1_DG/'

dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(atacInt))

presentGenes<-newGenes[(newGenes%in%annotGenes$gene_name)]

for ( i in presentGenes){
  try({
    p1<-CoveragePlot(
      object = atacInt ,
      region = i,
      features = i,
      expression.assay = "RNA",
      idents = clusters,
      extend.upstream = 60000,
      extend.downstream = 60000
    )
    outFile<-paste0(targetDir, 'CovPlot_AtacRna_CA1_DG_macs2_1.87_contrOnly_ext60k_', i, '_', curDate, '.jpeg')
    ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)
  })
}