library(ArchR)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRGenome("mm10")

setwd("/home/flyhunter/Wang/output/")

curDate <- Sys.Date()

addArchRGenome("mm10")

addArchRThreads(threads = 10) 

projFilt = readRDS("atacArchRMacs2IntRna")

getAvailableMatrices(projFilt)

archMetadat = data.frame(cbind(projFilt$cellNames, projFilt$Sample, projFilt$Annotations))
colnames(archMetadat) = c("CellName", "Group", "Annotations")

clusters = unique(archMetadat$Annotations)

saveArchRProject(ArchRProj = projFilt, outputDirectory = "archrRnaPeakPerClust/", load = T)

getPeakGenes = function(projFilt, cluster, archMetadat, targDir, curDate) {
  macs2Path <- '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
  try({
    subCell = archMetadat$CellName[archMetadat["Annotations"] == cluster]
    projSub = subsetCells(ArchRProj = projFilt, cellNames = subCell)
    
    projSub$MacsGroup = "Group1" # to call peaks per group probably need to remove C-R
    projSub<-addGroupCoverages(projSub, groupBy = 'MacsGroup', force =T)
    
    projMacs2 <-addReproduciblePeakSet(
      ArchRProj = projSub, 
      groupBy = "MacsGroup", 
      pathToMacs2 = macs2Path,
      peakMethod = 'Macs2',
      cutOff = 0.05,
      genomeSize = 1.87e+09, force = T
    )
    
    projSub <- addPeakMatrix(projMacs2, force = T)
    
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

targDir = "atacRna/peaksGenesCor/archRMacs2PeaksPerClust/"
p2GAllClust(projFilt = projFilt, clusters = clusters, archMetadat = archMetadat, targDir = targDir, curDate = curDate)