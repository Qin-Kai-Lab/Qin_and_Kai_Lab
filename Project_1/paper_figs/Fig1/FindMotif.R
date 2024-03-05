#BiocManager::install(version = "3.18")


library(Signac)
library(JASPAR2022)
library(Seurat)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
#library(JASPAR2020)
#library(ggseqlogo)
library(ggplot2)

markers = read.csv("atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_AllClusters_500K_2023-11-06_FiltClust.csv")

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")
DefaultAssay(atacFilt) <- "Combined_peaks"
#atacFilt[["PredictActivity"]] = NULL
gc()
# atacFilt <- FindTopFeatures(atacFilt, min.cutoff = 5)
# atacFilt <- RunTFIDF(atacFilt)
# 
# atacFilt<- RunSVD(atacFilt)

targDir = "Paper_figs/TasksList/P19/"
dir.create(targDir, recursive = T, showWarnings = F)

getMotif = function(atacFilt, clusters) {
  combDat = data.frame()
  for (curCluster in clusters) {
    try({
      curMarkers = unique(markers$peak[markers$Cluster == curCluster])
      enriched.motifs <- FindMotifs(
        object = atacFilt,
        features = curMarkers, assay= "Combined_peaks")
      enriched.motifs$Cluster = curCluster
      combDat = rbind(combDat, enriched.motifs)
    })
  }
  return(combDat)
}

getTopMotif = function(curDf, clusters, topN) {
  combDat = data.frame()
  for (cluster in clusters) {
    dfSel = curDf[curDf$Cluster == cluster,]
    dfSel = dfSel[order(dfSel$pvalue, decreasing = F),]
    topDf = head(dfSel, topN)
    combDat = rbind(combDat, topDf)
  }
  return(combDat)
}

clusters = unique(atacFilt$Annotations)
clusters = clusters[!clusters == "C-R"]

# vertebrate motifs
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


atacFilt <- AddMotifs(
  object = atacFilt,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


all(unique(markers$peak) %in% rownames(atacFilt))


# save with motiffs
#saveRDS(atacFilt, "atacIntegrated_macs2_2_RNA")

#atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)


#write.csv(enriched.motifs, "atacIntegrated_macs2_2_RNA_Motiffs_Jaspar2022.csv")


motifVert = getMotif(atacFilt = atacFilt, clusters = clusters )


write.csv(motifVert, paste0(targDir, "peakGeneCor_AllClusters_500K_VertMotifs_2023-12-18.csv"), row.names = F)

# try mouse specific motifs

pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = 10090, all_versions = FALSE)
)


atacFilt <- AddMotifs(
  object = atacFilt,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

#pfm <- getMatrixByName(JASPAR2022, species=10090)
motifVert = getMotif(atacFilt = atacFilt, clusters = clusters )
write.csv(motifVert, paste0(targDir, "peakGeneCor_AllClusters_500K_MouseMotifs_2023-12-19.csv"), row.names = F)

# select top motifs

vMotif = read.csv(paste0(targDir, "peakGeneCor_AllClusters_500K_VertMotifs_2023-12-18.csv"))
mMotif = read.csv(paste0(targDir, "peakGeneCor_AllClusters_500K_MouseMotifs_2023-12-19.csv"))

vTop = getTopMotif(curDf=vMotif, clusters=clusters, topN=10)
mTop = getTopMotif(curDf=mMotif, clusters=clusters, topN=10)

# save
write.csv(vTop, paste0(targDir, "peakGeneCor_AllClusters_500K_Top10_VertMotifs_2023-12-19.csv"), row.names = F)
write.csv(mTop, paste0(targDir, "peakGeneCor_AllClusters_500K_Top10_MouseMotifs_2023-12-19.csv"), row.names = F)


# make motifs for human, rat, and mouse

pfm1 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(10090), all_versions = FALSE))

pfm2 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(9606), all_versions = FALSE))

pfm3 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(10116), all_versions = FALSE))

curNames = unique(c(names(pfm1), names(pfm2), names(pfm3)))

opts = list()
opts[["ID"]] <- curNames

PFMatrixList <- getMatrixSet(JASPAR2022, opts)
length(PFMatrixList) == length(curNames)

atacFilt <- AddMotifs(
  object = atacFilt,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = PFMatrixList
)

motifComb = getMotif(atacFilt = atacFilt, clusters = clusters )

#write.csv(motifComb, paste0(targDir, "peakGeneCor_AllClusters_500K_CombMotifs_2023-12-22.csv"), row.names = F)

cTop = getTopMotif(curDf=motifComb, clusters=clusters, topN=10)

# save
#write.csv(cTop, paste0(targDir, "peakGeneCor_AllClusters_500K_Top10_CombinedMotifs_2023-12-22.csv"), row.names = F)

for (cluster in clusters) {
  curMotifs = cTop[cTop$Cluster == cluster,]
  curMotifs = curMotifs[order(curMotifs$p.adjust, decreasing = F),]
  plotMotif =  curMotifs$motif
  curPlot = MotifPlot(
    object = atacFilt,
    motifs = plotMotif, nrow = 2) &
    theme(text = element_text(size = 28)) +
    theme(axis.text.x = element_text(color = "black", size = 12))
  
  png(filename = paste0(targDir, "peakGeneCor_500K_Top10_CombinedMotifs_", cluster,  "_2023-12-22.png"), width = 32, height = 16, units = "in", res = 300)
  print(curPlot )
  dev.off()
  
}

