library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

set.seed(13)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")


atacInt = readRDS('atacMerged_macs2_2')

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)
rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )

atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

identical(colnames(atacFilt), colnames(rnaFilt))
identical(rownames(atacFilt@meta.data), rownames(rnaFilt@meta.data))

atacFilt$MonocClust = rnaFilt$MonocClust

##
atacFilt

atacFilt <- FindTopFeatures(atacFilt, min.cutoff = 10)
atacFilt <- RunTFIDF(atacFilt)
atacFilt <- RunSVD(atacFilt)
atacFilt <- RunUMAP(atacFilt, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(atacFilt, group.by = "MonocClust")
p1
#
table(atacFilt@meta.data$MonocClust)
#
Idents(atacFilt) <- 'MonocClust'
atacFilt$MonocClust <- factor(atacFilt$MonocClust, levels=c("1", "2", "3", "4"))
DimPlot(atacFilt)
#
atacFilt[['RNA']]<-rnaFilt@assays$RNA
DefaultAssay(atacFilt)<- 'Combined_peaks'
# annotate
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacFilt) <- annotations

#saveRDS(atacFilt, 'OpcOdc_allCelIntMacs2_2')

rm(atacInt, RNA.combined.norm)
gc()

targDir <- './OPC_ODC/Markers_Comp/Signac/'
dir.create(targDir, recursive = T)
#

atacMarkers <- FindAllMarkers(atacFilt , only.pos = T, min.pct = 0.05, test.use = 'LR', 
                           latent.vars = c('TSS.enrichment'), logfc.threshold = 0.25)

closGenes<-ClosestFeature(atacFilt, regions = rownames(atacMarkers))
allMarkEdit<-cbind(atacMarkers, closGenes)

outAtac<-paste0(targDir, 'atacPosMarkAll_macs2_2.csv')
write.csv(allMarkEdit, outAtac, row.names = F)

# split results table by cluster
splitResults <- function(x, var) {
  splitVars <- unique(x[, var])
  for (i in splitVars) {
    selDf <- x[(x[, var] == i),]
    outfile <- paste0(targDir, 'atacPosMarkAll_macs2_2_', i, '_', curDate, '.csv')
    write.csv(selDf, outfile, row.names = F)
  }
}

splitResults(x = allMarkEdit, var = 'cluster')

## check how many of the top RNA 
rnaMark = read.csv('./OPC_ODC/atacRna/rnaPosMark_1vsAll_monocClust_macs2.csv')
rnaOrd <- rnaMark[order(rnaMark$avg_log2FC, decreasing = T),]
rnaMarkAdj <- rnaOrd[(rnaOrd$p_val_adj < 0.05),]
top15 <- head(rnaMarkAdj, 15)
# atacMark
clust1 <- allMarkEdit[(allMarkEdit$cluster == 1),]
clust1Ord <- clust1[order(clust1$avg_log2FC, decreasing = T),]
clustAdj <- clust1Ord[(clust1Ord$p_val_adj < 0.05),]

tabMatchTop <- clustAdj[(clustAdj$gene_name %in% top15$Genes),]
nrow(tabMatchTop)
tabMatch <- clustAdj[(clustAdj$gene_name %in% rnaMarkAdj$Genes),]
nrow(tabMatch)

# for cluster 1 43/338 sign adj p, out of them 1 in top 15 RNA and 17 in significant adjusted RNA (557 sign adj)

# not adjusted p atac

clustAdj <- clust1Ord[(clust1Ord$p_val < 0.05),]
nrow(clustAdj)

tabMatchTop <- clustAdj[(clustAdj$gene_name %in% top15$Genes),]
nrow(tabMatchTop)
tabMatch <- clustAdj[(clustAdj$gene_name %in% rnaMarkAdj$Genes),]
nrow(tabMatch)

# for cluster 1 338/338 sign  p, out of them 6 in top 15 RNA and 91 in significant adjusted RNA (557 sign adj)