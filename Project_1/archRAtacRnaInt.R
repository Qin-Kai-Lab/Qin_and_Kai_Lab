library(ArchR)
library(MAST)
library(BSgenome.Mmusculus.UCSC.mm10)

curDate <- Sys.Date()

addArchRGenome("mm10")

addArchRThreads(threads = 10) 

inputFiles <- c("/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz", "/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz")

names(inputFiles)<- c('control', 'stress')

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archContrStress",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

saveArchRProject(ArchRProj = proj, outputDirectory = "archContrStressRna")

curDate<-Sys.Date()

addArchNames = function(x) {
  x$archName <- NA
  for (i in 1:nrow(x)) {
    if (x$group[i] == 'Control') {
      x$archName[i] = paste0('control#', x$CellName[i])
    } else if (x$group[i] == 'Stress') {
      x$archName[i] = paste0('stress#', x$CellName[i])
    }
  }
  return(x)
}

rnaMetadat = read.csv('RnaSeqMetadat.csv')
archMeta = data.frame(proj$cellNames)
colnames(archMeta)[1] = 'archName'

rnaMetadat = addArchNames(rnaMetadat)

atacRnaPres = merge(archMeta, rnaMetadat, by = 'archName', sort = F)
nrow(atacRnaPres)
length(archMeta$archName[archMeta$archName %in% rnaMetadat$archName])

projFilt = subsetCells(ArchRProj = proj, cellNames = atacRnaPres$archName)

identical(projFilt$cellNames, atacRnaPres$archName)

projFilt$Annotations = atacRnaPres$Annotations

getAvailableMatrices(ArchRProj = projFilt)

projFilt <- addIterativeLSI(ArchRProj = projFilt, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
projFilt <- addClusters(input = projFilt, reducedDims = "IterativeLSI", force = T)

rm(proj)
gc()
## add RNA
source('../programs/renameClusters.R')
archCellNames = function(x) {
  allNames = character()
  for ( i in 1:nrow(x@meta.data)) {
    if (x@meta.data$group[i] == "Control") {
      newName = paste0('control#', x@meta.data$CellName[i])
    } else if (x@meta.data$group[i] == "Stress") {
      newName = paste0('stress#', x@meta.data$CellName[i])
    }
    allNames = c(allNames, newName)
  }
  return(allNames)
}

archNames = archCellNames(x = RNA.combined.norm)

RNA.combined.norm = RenameCells(RNA.combined.norm, new.names = archNames )
RNA.combined.norm$archName = rownames(RNA.combined.norm@meta.data)

archPresent = archNames[(archNames%in%projFilt$cellNames)]
length(archPresent)
length(projFilt$cellNames)

rnaSub = subset(x =RNA.combined.norm, subset =  archName %in% archPresent )


rna_data <- as.matrix(GetAssayData(rnaSub, assay="RNA"))

rnaSub$Integ_Clust = rnaSub$Annotations

projFilt <-addGeneIntegrationMatrix(projFilt, seRNA = rnaSub, groupRNA = "Integ_Clust", groupATAC = "Annotations", force = TRUE)


#projFilt = readRDS("atacArchRMacs2IntRna")
# add peaks, need to check how to increase maximum number of cells for peak calling
projFilt$MacsGroup = "Group1" # to call peaks per group probably need to remove C-R
projFilt<-addGroupCoverages(projFilt, groupBy = 'MacsGroup', force =T)

macs2Path <- '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'

#rm(archMeta, archCellNames, archPresent, atacRnaPres, doubScores, rnaSub)
#gc()


getArchRThreads()

projMacs2 <-addReproduciblePeakSet(
  ArchRProj = projFilt, 
  groupBy = "MacsGroup", 
  pathToMacs2 = macs2Path,
  peakMethod = 'Macs2',
  genomeSize = 1.87e+09, force = T
)

projFilt <- addPeakMatrix(projMacs2, force = T)

getAvailableMatrices(projFilt)

# connect peaks and rna
projFilt <- addPeak2GeneLinks(
  ArchRProj = projFilt,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projFilt ,
  corCutOff = 0,
  resolution = 1,
  returnLoops = F,
  FDRCutOff =1
)




#saveRDS(projFilt, "atacArchRMacs2IntRna")

#saveArchRProject(ArchRProj = projFilt)


p2geneDF <- metadata(projFilt@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2GeneDFComplete = data.frame(p2geneDF)

targDir = "atacRna/peaksGenesCor/archRMacs2PeaksAllCells/"
write.csv(p2GeneDFComplete, paste0(targDir, "AllCelss_peaksRnaCor_", curDate, ".csv"), row.names = F)