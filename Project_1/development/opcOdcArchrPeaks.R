library(EnsDb.Mmusculus.v79)
library(archR)


targDir <- './OPC_ODC/Markers_Comp/archR/Macs2Peaks/'
dir.create(targDir, recursive = T)

curDate<-Sys.Date()


proj = loadArchRProject(path = "./archContrStress")

proj$PeakGr <- 'All'

proj<-addGroupCoverages(proj, groupBy = 'PeakGr')

## not using now
projPeaks <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "PeakGr",
  peakMethod = "Tiles",
  method = "p"
)
##

macs2Path <- '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'

projMacs2 <-addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "PeakGr", 
  pathToMacs2 = macs2Path,
  peakMethod = 'Macs2',
  genomeSize = 1.87e+09
)

projMacs2Peak <- addPeakMatrix(projMacs2)

getAvailableMatrices(projMacs2Peak)


rm(proj, projMacs2)
gc()
###

proj <- projMacs2Peak
getAvailableMatrices(proj)

###

opcStress <- read.table('OPC_ODC_stressCells.txt', T)
opcControl <- read.table('OPC_ODC_contrCells.txt', T)
opcControl$newId <- paste0('control#', opcControl$CellName)
opcStress$newId <- paste0('stress#', opcStress$CellName)

clustInfo <- rbind(opcStress, opcControl)

###

opcOdcPres <- clustInfo[(clustInfo$newId %in%proj$cellNames),]
opcOdcArch = subsetCells(ArchRProj = proj, cellNames = opcOdcPres$newId)
identical(opcOdcPres$newId, opcOdcArch$cellNames)
opcOdcPres$MonocClust <- as.factor(opcOdcPres$MonocClust)
opcOdcArch$MonocClust <- opcOdcPres$MonocClust

opcOdcArch <- addIterativeLSI(ArchRProj = opcOdcArch, useMatrix = "PeakMatrix", name = "IterativeLSI", force = TRUE)
opcOdcArch <- addClusters(input = opcOdcArch, reducedDims = "IterativeLSI")
##
opcOdcArch <- addUMAP(ArchRProj = opcOdcArch, reducedDims = "IterativeLSI")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "MonocClust", embedding = "UMAP")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
###

# save
saveArchRProject(ArchRProj = opcOdcArch, outputDirectory = "opcOdcArchMacs2")

##
# find markers
opcOdcArch$MonocClust <- as.character(opcOdcArch$MonocClust)

#opcOdcArch<-addPeakAnnotations(opcOdcArch)

atacMark <- getMarkerFeatures(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", useMatrix = 'PeakMatrix'
)

markerList <- getMarkers(atacMark, cutOff = "Pval < 0.05 & Log2FC >  0 ", returnGR = T)
clust1Mark <- data.frame(markerList$`1`)
granges<-markerList$`1`

regions<- paste(clust1Mark$seqnames, clust1Mark$start, clust1Mark$end, sep ='-')


# annotate
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

nearest_feature <- distanceToNearest(x = granges, subject = annotations)
feature_hits <- annotations[subjectHits(x = nearest_feature)]
df <- as.data.frame(x = mcols(x = feature_hits))
df$closest_region <- GRangesToString(grange = feature_hits)
df$query_region <- GRangesToString(grange = granges)
df$distance <- mcols(x = nearest_feature)$distance


write.csv(df, paste0(targDir, 'atacPosMarkAll_', '1', '_', curDate, '.csv'), row.names = F)
#

atacGenes<- unique(df$gene_name)

## check how many of the top RNA 
rnaMark = read.csv('./OPC_ODC/atacRna/rnaPosMark_1vsAll_monocClust_macs2.csv')
rnaOrd <- rnaMark[order(rnaMark$avg_log2FC, decreasing = T),]
rnaMarkAdj <- rnaOrd[(rnaOrd$p_val_adj < 0.05),]
top15 <- head(rnaMarkAdj, 15)

tabMatchTop <- df[(df$gene_name %in% top15$Genes),]
length(unique(tabMatchTop$gene_name))
tabMatch <- df[(df$gene_name %in% rnaMarkAdj$Genes),]
length(unique(tabMatch$gene_name))
