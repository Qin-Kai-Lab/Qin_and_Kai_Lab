library(ArchR)
library(MAST)

curDate <- Sys.Date()

addArchRGenome("mm10")

addArchRThreads(threads = 14) 

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

saveArchRProject(ArchRProj = proj, outputDirectory = "archContrStress")

proj_filt <- filterDoublets(ArchRProj = proj) # doublets filter removes too many cells, over 200 for OPC

###

opcStress <- read.table('OPC_ODC_stressCells.txt', T)
opcControl <- read.table('OPC_ODC_contrCells.txt', T)
opcControl$newId <- paste0('control#', opcControl$CellName)
opcStress$newId <- paste0('stress#', opcStress$CellName)

clustInfo <- rbind(opcStress, opcControl)

###

opcOdcPres <- clustInfo[(clustInfo$newId %in%proj_filt$cellNames),]
opcOdcArch = subsetCells(ArchRProj = proj_filt, cellNames = opcOdcPres$newId)
identical(opcOdcPres$newId, opcOdcArch$cellNames)
opcOdcPres$MonocClust <- as.factor(opcOdcPres$MonocClust)
opcOdcArch$MonocClust <- opcOdcPres$MonocClust

opcOdcArch <- addIterativeLSI(ArchRProj = opcOdcArch, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
opcOdcArch <- addClusters(input = opcOdcArch, reducedDims = "IterativeLSI")
##
opcOdcArch <- addUMAP(ArchRProj = opcOdcArch, reducedDims = "IterativeLSI")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "MonocClust", embedding = "UMAP")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
###
targDir<-'./OPC_ODC/atacRna/'
rnaMarkers<- read.csv(paste0(targDir, 'rnaPosMark_1vsAll_monocClust_macs2.csv'))
rnaMarkOrd <- rnaMarkers[order(rnaMarkers$avg_log2FC, decreasing = T),]
rnaMarkTop <- head(rnaMarkOrd, 12)
genes<-c('Plp1', 'Ptgds', 'Trf', 'Calcrl', 'Cspg4', 'Vcan')
rnaMarkComb <- unique(c(rnaMarkTop$Genes, genes))
length(rnaMarkComb)
#
archGenes <- getGenes(ArchRProj = opcOdcArch)
allArchGenes <- archGenes$symbol
presMark <- rnaMarkComb[rnaMarkComb%in%allArchGenes]
length(presMark)
##
# manual top markers
presMark
i = "Prr5l"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 180000,
  downstream = 60000
)

targDir<-'./OPC_ODC/atacRna/archR/DoubFilt/'
dir.create(targDir, recursive = T)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()
##
presMark
i = "Rnf220"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 290000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Slc24a2"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 290000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Plp1"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 10000,
  downstream = 25000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Fnbp1"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 40000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "St18"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 550000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Ptgds"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 10000,
  downstream = 10000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Trf" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 60000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

presMark
i = "Pex5l"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 250000,
  downstream = 80000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Mobp"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 30000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Grm3"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 0,
  geneSymbol = i,
  upstream = 300000,
  downstream = 90000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Mag"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 0,
  geneSymbol = i,
  upstream = 40000,
  downstream = 30000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Calcrl" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Cspg4" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 40000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Vcan" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 140000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
## No doublet filtering
opcOdcPres <- clustInfo[(clustInfo$newId %in%proj$cellNames),]
opcOdcArch = subsetCells(ArchRProj = proj, cellNames = opcOdcPres$newId)
identical(opcOdcPres$newId, opcOdcArch$cellNames)
opcOdcPres$MonocClust <- as.factor(opcOdcPres$MonocClust)
opcOdcArch$MonocClust <- opcOdcPres$MonocClust

opcOdcArch <- addIterativeLSI(ArchRProj = opcOdcArch, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
opcOdcArch <- addClusters(input = opcOdcArch, reducedDims = "IterativeLSI")
##
opcOdcArch <- addUMAP(ArchRProj = opcOdcArch, reducedDims = "IterativeLSI")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "MonocClust", embedding = "UMAP")
plotEmbedding(ArchRProj = opcOdcArch, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#

length(presMark)

##
# manual top markers
##
# manual top markers
targDir<-'./OPC_ODC/atacRna/archR/'

presMark
i = "Prr5l"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 180000,
  downstream = 60000
)


fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()
##
presMark
i = "Rnf220"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 290000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Slc24a2"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 290000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Plp1"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 10000,
  downstream = 25000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Fnbp1"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 40000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "St18"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 550000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Ptgds"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 10000,
  downstream = 10000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Trf" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 60000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

presMark
i = "Pex5l"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 250000,
  downstream = 80000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Mobp"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 30000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Grm3"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 0,
  geneSymbol = i,
  upstream = 300000,
  downstream = 90000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()


##
presMark
i = "Mag"
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 0,
  geneSymbol = i,
  upstream = 40000,
  downstream = 30000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Calcrl" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 150000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Cspg4" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 40000,
  downstream = 60000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()

##
presMark
i = "Vcan" 
p3<-plotBrowserTrack(
  ArchRProj = opcOdcArch, 
  groupBy = "MonocClust", 
  minCells = 1,
  geneSymbol = i,
  upstream = 140000,
  downstream = 50000
)

fName <- paste0(targDir, 'OpcOdc_CovPlot_AtacRna_MonocClust_1vsAll_macs2_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()