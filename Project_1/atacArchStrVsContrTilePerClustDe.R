library(ArchR)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Signac)

curDate = Sys.Date()

addArchRThreads(threads = 14) 

projFilt = readRDS("atacArchRnaFilt")
addArchRGenome("mm10")

#plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "newClust", embedding = "UMAP")+
#  theme_ArchR(baseSize = 20,legendTextSize = 20)


projFilt<-addGroupCoverages(projFilt, groupBy = 'Annotations')

TilesPath <- '/home/flyhunter/miniconda3/envs/Tiles/bin/Tiles'

projTiles <-addReproduciblePeakSet(
  ArchRProj = projFilt, 
  groupBy = "Annotations",
  peakMethod = 'Tiles',
  cutOff = 0.001,
  genomeSize = 1.87e+09, force = T
)

projFilt <- addPeakMatrix(projTiles, force = T)

getAvailableMatrices(projFilt)

#saveRDS(projFilt, 'atacArchRnaFiltTilesPerClust')

rm(projTiles)

# differential expression
targDir <- './archContrStress/atacRna/StressVsControl/Tiles_0.001_PerCluster/'
dir.create(targDir)

atacMark = function() {
  allClustSum = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(projFilt$Annotations)
  for ( cluster in clusters) {
    cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
    cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
    cellsClust = addIterativeLSI(ArchRProj = cellsClust, useMatrix = "PeakMatrix", name = "IterativeLSITiles", force = TRUE)
    atacMark <- getMarkerFeatures(
      ArchRProj = cellsClust, 
      groupBy = "Group",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      useGroups = "stress",
      bgdGroups = "control",
      maxCells = 20000,
      useMatrix = "PeakMatrix"
    )
    manSum = data.frame(cbind(data.frame(atacMark@elementMetadata), atacMark@assays@data$Pval, atacMark@assays@data$FDR, atacMark@assays@data$Log2FC)) # checked manual sumary is identical with package
    colnames(manSum)[5:7] = c('Pval', 'FDR_def', 'Log2FC')
    manSum$bonferonni_p = p.adjust(p=manSum$Pval, method = 'bonferroni')
    manSum$fdr_p = p.adjust(p=manSum$Pval, method = 'fdr')
    manSum$Cluster = cluster
    outName = paste0(targDir, 'StressVsControl_', cluster, '_', curDate, '_TilesPerCluster.csv')
    write.csv(manSum, outName, row.names = F)
    allClustSum = rbind(allClustSum, manSum)
  }
  return(allClustSum)
}

allClustMark = atacMark()
write.csv(allClustMark, paste0(targDir, 'StressVsControl_TilesPerClust_allClusters', curDate, '.csv'), row.names = F)

#allClustMark = read.csv(paste0(targDir, 'StressVsControl_Tiles_allClusters', curDate, '.csv'))

# get annotations
targDir <- './archContrStress/atacRna/StressVsControl/Tiles_0.001_PerCluster/Annotations/'
dir.create(targDir ,recursive = T)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

atacAnnot = function() {
  allClustSum = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(projFilt$Annotations)
  for ( cluster in clusters) {
    cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
    cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
    cellsClust = addIterativeLSI(ArchRProj = cellsClust, useMatrix = "PeakMatrix", name = "IterativeLSITiles", force = TRUE)
    atacMark <- getMarkerFeatures(
      ArchRProj = cellsClust, 
      groupBy = "Group",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      useGroups = "stress",
      bgdGroups = "control",
      maxCells = 20000,
      useMatrix = "PeakMatrix"
    )
    markerList <- getMarkers(atacMark, cutOff = "Pval <= 1", returnGR = T)
    granges<-markerList$stress
    nearest_feature <- distanceToNearest(x = granges, subject = annotations)
    feature_hits <- annotations[subjectHits(x = nearest_feature)]
    df <- as.data.frame(x = mcols(x = feature_hits))
    df$closest_region <- GRangesToString(grange = feature_hits)
    df$query_region <- GRangesToString(grange = granges)
    df$distance <- mcols(x = nearest_feature)$distance
    df$Log2FC <- granges$Log2FC
    df$Cluster <- cluster
    allClustSum = rbind(allClustSum, df)
    outName = paste0(targDir, 'StressVsControl_', cluster, '_', curDate, '_TilesPerCluster.csv')
    write.csv(df, outName, row.names = F)
  }
  return(allClustSum)
}

allMarkAnnot = atacAnnot()

write.csv(allMarkAnnot, paste0(targDir, 'StressVsControl_Tiles_allClusters', curDate, '.csv'), row.names = F)

#allMarkAnnot = read.csv(paste0(targDir, 'StressVsControl_Tiles_allClusters', curDate, '.csv'))

allMarkAnnot$Reg_Clust = paste(allMarkAnnot$query_region, allMarkAnnot$Cluster, sep = "-")

length(unique(allMarkAnnot$Reg_Clust))

allClustMark$Reg_Clust <- paste(allClustMark$seqnames, allClustMark$start, allClustMark$end, allClustMark$Cluster, sep = "-")

mergeMark = merge(allClustMark, allMarkAnnot, by = "Reg_Clust")

targDir <- './archContrStress/atacRna/StressVsControl/Tiles_0.001_PerCluster/Annotated/'
dir.create(targDir)
write.csv(mergeMark, paste0(targDir, 'StressVsControl_Tiles_allClusters', curDate, '.csv'), row.names = F)

clusters = unique(mergeMark$Cluster.x)

for ( cluster in clusters ) {
  curDf = mergeMark[(mergeMark$Cluster.x == cluster),]
  outName = paste0(targDir, 'StressVsControl_Tiles_', cluster, '_', curDate, '.csv')
  write.csv(curDf, outName, row.names = F)
}