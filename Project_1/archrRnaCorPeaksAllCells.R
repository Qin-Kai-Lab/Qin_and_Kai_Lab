library(ArchR)

setwd("/home/flyhunter/Wang/output/")

curDate <- Sys.Date()

addArchRGenome("mm10")

addArchRThreads(threads = 10) 

projFilt = readRDS("atacArchRMacs2IntRna")

getAvailableMatrices(projFilt)

archMetadat = data.frame(cbind(projFilt$cellNames, projFilt$Sample, projFilt$Annotations))
colnames(archMetadat) = c("CellName", "Group", "Annotations")

clusters = unique(archMetadat$Annotations)

getPeakGenes = function(projFilt, cluster, archMetadat, targDir, curDate) {
  try({
    subCell = archMetadat$CellName[archMetadat["Annotations"] == cluster]
    projSub = subsetCells(ArchRProj = projFilt, cellNames = subCell)
    # connect peaks and rna
    projSub <- addPeak2GeneLinks(
      ArchRProj = projSub,
      reducedDims = "IterativeLSI")
    
    p2g <- getPeak2GeneLinks(
      ArchRProj = projSub  ,
      corCutOff = 0,
      resolution = 1,
      returnLoops = F,
      FDRCutOff =1
    )
    
    p2geneDF <- metadata(projSub @peakSet)$Peak2GeneLinks
    p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
    p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
    p2GeneDFComplete = data.frame(p2geneDF)
    p2GeneDFComplete$Clusters = cluster
    write.csv(p2GeneDFComplete, paste0(targDir, cluster, "_peaksRnaCor_", curDate, ".csv"), row.names = F)
  })
}

p2GAllClust = function(projFilt, clusters, archMetadat, targDir, curDate) {
  dir.create(targDir, recursive = T)
  for (cluster in clusters) {
    getPeakGenes(projFilt = projFilt, cluster = cluster, archMetadat = archMetadat, targDir = targDir, curDate = curDate)
  }
}

targDir = "atacRna/peaksGenesCor/archRMacs2PeaksAllCells/"
p2GAllClust(projFilt = projFilt, clusters = clusters, archMetadat = archMetadat, targDir = targDir, curDate = curDate)