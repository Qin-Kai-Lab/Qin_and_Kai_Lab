library(Seurat)
library(Signac)
#library(future)
#plan("multicore", workers = 5)

#options(future.globals.maxSize = 50 * 1024 ^ 3)

#options(future.globals.seed=TRUE)

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/spearman/"
dir.create(targetDir, recursive = T, showWarnings = F)

getPeakRna = function(atacFilt, test, targetDir) {
  try({
    peakRna = LinkPeaks(
      object = atacFilt,
      peak.assay = 'Combined_peaks',
      expression.assay = "RNA",
      peak.slot = "counts",
      expression.slot = "data",
      method = test,
      distance = 5e+05,
      min.cells = 10,
      pvalue_cutoff = 0.05,
      score_cutoff = 0.05,
      verbose = TRUE)
    p2g = data.frame(Links(peakRna))
    write.csv(p2g, file = paste0(targetDir, "peakGeneCor_AllClusters_500K_2023-11-06.csv"), row.names = F)
    curPath = paste0(targetDir, "peakGeneCor_AllClusters_500K_2023-11-06.csv")
    print(curPath)
  })
}

getPeakRna(atacFilt= atacFilt, test = "spearman", targetDir = targetDir)


# edit peaks

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allMarkers = allMarkers[allMarkers$avg_log2FC > 0,]
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]


clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'MG', 'ODC', 'OPC', 'SUB')

corDf = read.csv("atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_AllClusters_500K_2023-11-06.csv")
corDf = corDf[corDf$pvalue < 0.01,]
corDf = corDf[order(corDf$pvalue, decreasing = F),]

combDf = data.frame()
for (cluster in clusters ) {
  rnaMarkDf = allMarkersPct[allMarkersPct$cluster == cluster,]
  rnaGenes = unique(rnaMarkDf$gene[rnaMarkDf$p_val_adj < 0.05])
  corSub = corDf[corDf$gene%in%rnaGenes,]
  corSub$Cluster = cluster
  combDf = rbind(combDf, corSub)
}

#
# combDf[combDf$peak == "chr7-115666776-115667244",]
# combDf[combDf$peak == "chr15-54935460-54935838",]
#

editTable = function(curDf) {
    curDf$position = gsub("chr.{1,2}-", "", curDf$peak)
    selDf = curDf[, c('peak', 'seqnames', 'position', 'gene', 'pvalue', 'Cluster')]
    colnames(selDf)[2] = "Chromosome"
    return(selDf)
}

combEdit = editTable(combDf)

write.csv(combEdit, "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_AllClusters_500K_2023-11-06_FiltClust.csv", row.names = F)
write.csv(corDf, "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_AllClusters_500K_2023-11-06_Filtered.csv", row.names = F)