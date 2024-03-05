library(ArchR)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Signac)

curDate = Sys.Date()

addArchRThreads(threads = 14) 

projFilt = readRDS("atacArchRnaFilt")
addArchRGenome("mm10")

plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "newClust", embedding = "UMAP")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)

projFilt$PeakGr <- 'All'
projFilt<-addGroupCoverages(projFilt, groupBy = 'PeakGr', force = T)

macs2Path <- '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'

projMacs2 <-addReproduciblePeakSet(
  ArchRProj = projFilt, 
  groupBy = "PeakGr", 
  pathToMacs2 = macs2Path,
  peakMethod = 'Macs2',
  genomeSize = 1.87e+09
)

projFilt <- addPeakMatrix(projMacs2)

getAvailableMatrices(projFilt)

rm(projMacs2)

projFilt <- addIterativeLSI(ArchRProj = projFilt, useMatrix = "PeakMatrix", name = "IterativeLSIMacs2", force = TRUE)
projFilt <- addUMAP(ArchRProj = projFilt, reducedDims = "IterativeLSIMacs2", name = 'UMAP_macs2', force = T)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "newClust", embedding = "UMAP_macs2")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)

targDir <- './archContrStress/atacRna/UMAPs/'

# save Umaps not finished
fName <- paste0(targDir, 'atacFiltRna_Macs2_Umap_annot', '_', curDate, '.png')
png(file= fName, width=14, height=14, units = 'in', res = 300)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "newClust", embedding = "UMAP")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)
dev.off()

fName <- paste0(targDir, 'atacFiltRna_Macs2_Umap_Group', '_', curDate, '.png')
png(file= fName, width=14, height=14, units = 'in', res = 300)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "Group", embedding = "UMAP")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)
dev.off()

# stress vs control
targDir <- './archContrStress/atacRna/StressVsControl/'

atacMark = function() {
  allClustSum = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(projFilt$Annotations)
  for ( cluster in clusters) {
    cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
    cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
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
    #clust1Mark <- data.frame(markerList$stress)
    #clust1Mark$fdr_p <- atacMark@assays@data$FDR
    manSum = data.frame(cbind(data.frame(atacMark@elementMetadata), atacMark@assays@data$Pval, atacMark@assays@data$FDR, atacMark@assays@data$Log2FC)) # checked manual sumary is identical with package
    colnames(manSum)[5:7] = c('Pval', 'FDR_def', 'Log2FC')
    #sumSel = manSum[(manSum$Log2FC > 0),]
    manSum$bonferonni_p = p.adjust(p=manSum$Pval, method = 'bonferroni')
    manSum$fdr_p = p.adjust(p=manSum$Pval, method = 'fdr')
    manSum$Cluster = cluster
    outName = paste0(targDir, 'StressVsControl_Macs2_', cluster, '_', curDate, '.csv')
    write.csv(manSum, outName, row.names = F)
    allClustSum = rbind(allClustSum, manSum)
  }
  return(allClustSum)
}

# check than manual summary is the same as archR
#q1 <- manSum[order(manSum$x.2),]
#q2 <- clust1Mark[order(clust1Mark$Log2FC),]
#identical(q1$start, q2$start)

allClustMark = atacMark()
write.csv(allClustMark, paste0(targDir, 'StressVsControl_Macs2_allClusters', curDate, '.csv'), row.names = F)

saveRDS(projFilt, 'atacArchRnaFilt')

allClustMark = read.csv(paste0(targDir, 'StressVsControl_Macs2_allClusters', curDate, '.csv'))

# get annotations
targDir <- './archContrStress/atacRna/StressVsControl/Annotations/'
dir.create(targDir ,recursive = T)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'

atacAnnot = function() {
  allClustSum = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(projFilt$Annotations)
  for ( cluster in clusters) {
    cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
    cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
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
    outName = paste0(targDir, 'StressVsControl_Macs2_', cluster, '_', curDate, '.csv')
    write.csv(df, outName, row.names = F)
  }
  return(allClustSum)
}

allMarkAnnot = atacAnnot()

write.csv(allMarkAnnot, paste0(targDir, 'StressVsControl_Macs2_allClusters', curDate, '.csv'), row.names = F)

allMarkAnnot = read.csv(paste0(targDir, 'StressVsControl_Macs2_allClusters', curDate, '.csv'))

allMarkAnnot$Reg_Clust = paste(allMarkAnnot$query_region, allMarkAnnot$Cluster, sep = "-")

length(unique(allMarkAnnot$Reg_Clust))

allClustMark$Reg_Clust <- paste(allClustMark$seqnames, allClustMark$start, allClustMark$end, allClustMark$Cluster, sep = "-")

mergeMark = merge(allClustMark, allMarkAnnot, by = "Reg_Clust")

targDir <- './archContrStress/atacRna/StressVsControl/Annotated/'
dir.create(targDir)
write.csv(mergeMark, paste0(targDir, 'StressVsControl_Macs2_allClusters', curDate, '.csv'), row.names = F)

clusters = unique(mergeMark$Cluster.x)

for ( cluster in clusters ) {
  curDf = mergeMark[(mergeMark$Cluster.x == cluster),]
  outName = paste0(targDir, 'StressVsControl_Macs2_', cluster, '_', curDate, '.csv')
  write.csv(curDf, outName, row.names = F)
}