library(ArchR)
library(MAST)

proj = loadArchRProject(path = "./archContrStress")

targDir <- './OPC_ODC/Markers_Comp/archR/'
dir.create(targDir, recursive = T)

curDate<-Sys.Date()

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

opcOdcArch <- addIterativeLSI(ArchRProj = opcOdcArch, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
opcOdcArch <- addClusters(input = opcOdcArch, reducedDims = "IterativeLSI")
##
opcOdcArch <- addUMAP(ArchRProj = opcOdcArch, reducedDims = "IterativeLSI", force = TRUE)
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "MonocClust", embedding = "UMAP")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

condGr <- gsub("\\#.*", "", opcOdcArch$cellNames)

opcOdcArch$Group <- condGr
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "Group", embedding = "UMAP")


opcOdcArch <- addHarmony(
  ArchRProj = opcOdcArch,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Group"
)

#opcOdcArch <- addUMAP(ArchRProj = opcOdcArch, reducedDims = "Harmony", force = T)
#plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "MonocClust", embedding = "UMAP")

saveArchRProject(ArchRProj = opcOdcArch, outputDirectory = "opcOdcArch")

# find markers
opcOdcArch$MonocClust <- as.character(opcOdcArch$MonocClust)

atacMark <- getMarkerFeatures(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(atacMark, cutOff = "FDR < 0.05 & Log2FC >  0 ")
clust1Mark <- data.frame(markerList$`1`)

## check how many of the top RNA 
rnaMark = read.csv('./OPC_ODC/atacRna/rnaPosMark_1vsAll_monocClust_macs2.csv')
rnaOrd <- rnaMark[order(rnaMark$avg_log2FC, decreasing = T),]
rnaMarkAdj <- rnaOrd[(rnaOrd$p_val_adj < 0.05),]
top15 <- head(rnaMarkAdj, 15)

tabMatchTop <- clust1Mark[(clust1Mark$name %in% top15$Genes),]
nrow(tabMatchTop)
tabMatch <- clust1Mark[(clust1Mark$name %in% rnaMarkAdj$Genes),]
nrow(tabMatch)
# for cluster 1  sign adj p 1056, out of them 11 in top 15 RNA and 126 in significant adjusted RNA (557 sign adj)

# split results table by cluster
splitResults <- function(x) {
  for (i in 1:4) {
    selDf <- as.data.frame(x[[i]])
    outfile <- paste0(targDir, 'atacPosMarkAll_', i, '_', curDate, '.csv')
    write.csv(selDf, outfile, row.names = F)
  }
}

splitResults(x = markerList)



#markerListP <- getMarkers(atacMark, cutOff = "Pval < 0.05 & Log2FC >  0 ")
#clust1MarkP <- data.frame(markerListP$`1`)