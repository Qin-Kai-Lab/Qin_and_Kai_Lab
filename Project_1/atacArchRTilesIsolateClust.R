library(ArchR)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Signac)

curDate = Sys.Date()

addArchRThreads(threads = 14)
TilesPath <- '/home/flyhunter/miniconda3/envs/Tiles/bin/Tiles'

projFilt = readRDS("atacArchRnaFilt")
clusters = unique(projFilt$Annotations)

addArchRGenome("mm10")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

getTilesPeaks = function(x) {
  x <-addGroupCoverages(x, groupBy = 'Annotations')
  projTiles <-addReproduciblePeakSet(
    ArchRProj = x, 
    groupBy = "Annotations",
    peakMethod = 'Tiles',
    cutOff = 0.001,
    genomeSize = 1.87e+09, force = T
  )
  x = addPeakMatrix(projTiles, force = T)
  return(x)
}

atacMark = function(x, targDir) {
  projFilt = x
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
  }
}

atacAnnot = function(x, targDir) {
  projFilt = x
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
    outName = paste0(targDir, 'StressVsControl_', cluster, '_', curDate, '_TilesPerCluster.csv')
    write.csv(df, outName, row.names = F)
  }
}

runAll = function() {
  for (cluster in clusters) {
    try({
      cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
      cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
      
      cellsClust = getTilesPeaks(x = cellsClust)
      # differential expression
      targDir <- './archContrStress/atacRna/StressVsControl/Tiles_0.001_IsolateClust/'
      dir.create(targDir)
      atacMark(x = cellsClust, targDir = targDir)
      # get annotations
      targDir <- './archContrStress/atacRna/StressVsControl/Tiles_0.001_IsolateClust/Annotations/'
      dir.create(targDir ,recursive = T)
      atacAnnot(x = cellsClust, targDir = targDir)
      
    })
  }
}

# rune pipeline
runAll()