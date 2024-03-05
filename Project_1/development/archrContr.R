#BiocManager::install("DirichletMultinomial")
#BiocManager::install("motifmatchr")
#BiocManager::install("chromVAR")

#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR)
library(MAST)
ArchR::installExtraPackages()

curDate <- Sys.Date()

addArchRGenome("mm10")

addArchRThreads(threads = 14) 

inputFiles <- c("/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = c('control'),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archControl",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj_filt <- filterDoublets(ArchRProj = proj)

proj <- proj_filt

#atacInt<-readRDS('atacAllFiltMacs2_control_1.87')

#clustInfo <- atacInt@meta.data[, c('gex_barcode', 'Annotations')]
#clustInfo$newId <- paste0('control#', clustInfo$gex_barcode)

source('../programs/renameClusters.R')

rnaContr<-subset(x = RNA.combined.norm, subset = group == 'Control')
clustInfo <- rnaContr@meta.data[, c('CellName', 'Annotations')]
clustInfo$newId <- paste0('control#', clustInfo$CellName)

idxPass <- which(proj$cellNames %in% clustInfo$newId)

cellsPass <- proj$cellNames[idxPass]
proj_rna <- proj[cellsPass, ]

proj_rna <- addIterativeLSI(ArchRProj = proj_rna, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
proj_rna <- addClusters(input = proj_rna, reducedDims = "IterativeLSI")

clustInfFiltr<- clustInfo[(clustInfo$newId %in% proj_rna$cellNames),]

identical(clustInfFiltr$newId, proj_rna$cellNames)

clustInfSort<- clustInfFiltr[match(proj_rna$cellNames, clustInfFiltr$newId),]

identical(clustInfSort$newId, proj_rna$cellNames)

proj_rna$Annotations <- clustInfSort$Annotations

proj_rna <- addUMAP(ArchRProj = proj_rna, reducedDims = "IterativeLSI")

plotEmbedding(ArchRProj = proj_rna, colorBy = "cellColData", name = "Annotations", embedding = "UMAP")

plotEmbedding(ArchRProj = proj_rna, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

# plot genomic regions 

p1<-plotBrowserTrack(
  ArchRProj = proj_rna, 
  groupBy = "Annotations", 
  geneSymbol = c('Meis2', 'Mpped1', 'Satb2'),
  upstream = 300000,
  downstream = 100000
)

grid::grid.draw(p1$Mpped1)

plotPDF(plotList = p1, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj_rna, 
        addDOC = FALSE, width = 18, height = 18)


p2<-plotBrowserTrack(
  ArchRProj = proj_rna, 
  groupBy = "Annotations", 
  geneSymbol = c('Meis2', 'Mpped1', 'Satb2'),
  upstream = 20000,
  downstream = 100000
)

grid::grid.draw(p2$Mpped1)

proj_rna <- saveArchRProject(ArchRProj = proj_rna)

# save ready plots

png(file="./archControl/atacRnaFilt_Mpped1.png", width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p2$Mpped1)
dev.off()

png(file="./archControl/atacRnaFilt_Meis2.png", width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p1$Meis2)
dev.off()

png(file="./archControl/atacRnaFilt_Satb2.png", width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p1$Satb2)
dev.off()

# find top RNA markers for control group
atacInt<-readRDS('atacAllFiltMacs2_control_1.87')
DefaultAssay(atacInt)<-'RNA'
rnaMarkers<-FindMarkers(object=atacInt, ident.1= 'CA1', only.pos = T, min.pct = 0.05, logfc.threshold = 0.25, test.use = "MAST")
rnaMarkers$Genes<-row.names(rnaMarkers)
outRna<-paste0('./atacRnaContr/macs2/rnaPosMarkCA1vsAll_macs2_1.87_control.csv')
write.csv(rnaMarkers, outRna, row.names = F)

rm(atacInt)
gc()

rnaMarkord = rnaMarkers[order(rnaMarkers$avg_log2FC, decreasing = T),]

topMark <- head(rnaMarkord$Genes, 15)

archGenes <- getGenes(ArchRProj = proj_rna)
allArchGenes <- archGenes$symbol
presMark <- topMark[topMark%in%allArchGenes]

# automatic plot making is not good in archR
p3<-plotBrowserTrack(
  ArchRProj = proj_rna, 
  groupBy = "Annotations", 
  geneSymbol = presMark,
  upstream = 900000,
  downstream = 400000
)

identical(names(p3), presMark)

for ( i in presMark) {
  fName <- paste0('./archControl/atacRnaFilt_', i, '_', curDate, '.png')
  png(file= fName, width=18, height=14, units = 'in', res = 300)
  grid::grid.draw(p3[[i]])
  dev.off()
}
  
# manual top markers
presMark
i = "Hcn1"
p3<-plotBrowserTrack(
  ArchRProj = proj_rna, 
  groupBy = "Annotations", 
  geneSymbol = i,
  upstream = 90000,
  downstream = 490000
)

fName <- paste0('./archControl/atacRnaFilt_', i, '_', curDate, '.png')
png(file= fName, width=18, height=14, units = 'in', res = 300)
grid::grid.draw(p3[[i]])
dev.off()