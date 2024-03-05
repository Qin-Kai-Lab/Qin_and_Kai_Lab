library(ArchR)
library(Signac)


source('../programs/renameClusters.R')


atac = readRDS('atacArchRnaFilt')



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

archPresent = archNames[(archNames%in%atac$cellNames)]
length(archPresent)
length(atac$cellNames)

rnaSub = subset(x =RNA.combined.norm, subset =  archName %in% archPresent )

rnaSub@meta.data$Annotations = as.character(rnaSub@meta.data$Annotations)

atacFilt = subsetCells(ArchRProj = atac, cellNames = archPresent)

rm(RNA.combined.norm, atac)
gc()

atacFilt <- addIterativeLSI(ArchRProj = atacFilt, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)

atacFilt <- addGeneIntegrationMatrix(
  ArchRProj = atacFilt,
  seRNA = rnaSub, groupRNA = 'archName', force = T, addToArrow = TRUE)

getAvailableMatrices(ArchRProj = atacFilt)

atacFilt <- addPeak2GeneLinks(
  ArchRProj = atacFilt,
  reducedDims = "IterativeLSI"
)


p2g = getPeak2GeneLinks(
  ArchRProj = atacFilt,
  corCutOff = 0.45,
  FDRCutOff = 1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = T
)
