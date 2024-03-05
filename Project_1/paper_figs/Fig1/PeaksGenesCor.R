library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)
library(Seurat)
library(Signac)
library(foreach)
library(doParallel)
library(GenomicRanges)

#clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'MG', 'ODC', 'OPC', 'SUB')

clusters<-c("ODC", "OPC")

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

# targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/"
# dir.create(targetDir, recursive = T, showWarnings = F)


targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/spearman/"
dir.create(targetDir, recursive = T, showWarnings = F)

getPeakRna = function(atacFilt, test, cluster, targetDir) {
  try({
    atacSub = subset(x = atacFilt, subset = Annotations == cluster)
    peakRna = LinkPeaks(
      object = atacSub,
      peak.assay = 'Combined_peaks',
      expression.assay = "RNA",
      peak.slot = "counts",
      expression.slot = "counts",
      method = test,
      distance = 5e+05,
      min.cells = 10,
      pvalue_cutoff = 0.05,
      score_cutoff = 0.05,
      verbose = TRUE)
    p2g = data.frame(Links(peakRna))
    p2g$Cluster = cluster
    write.csv(p2g, file = paste0(targetDir, "peakGeneCor_", cluster, "_500K_2023-10-31.csv"), row.names = F)
    pritn(cluster)
  })
  return(p2g)
}


useCores<-2

cl <- makeCluster(useCores, type = "FORK")
registerDoParallel(cl)

p2gComb<-data.frame(foreach(i=clusters, .combine=rbind, .packages=c('Seurat', 'Signac') ) %dopar%{
  getPeakRna(atacFilt= atacFilt, test = "spearman", cluster = i, targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/spearman/")
})

parallel::stopCluster(cl = cl)


#write.csv(p2gComb, file = paste0(targetDir, "peakGenes_All_20K_2023-10-27.csv"), row.names = F)


#curFiles = list.files(targetDir, pattern = ".csv", full.names = T)

#curDf = read.csv(curFiles[1])

editTables = function(curFiles, targetDir) {
  targDf = paste0(targetDir, "adjsted_p/")
  dir.create(targDf, recursive = T, showWarnings = F)
  for (curFile in curFiles) {
    sampleName = basename(curFile)
    curDf = read.csv(curFile)
    curDf$fdr_p = p.adjust(curDf$pvalue, method = "fdr")
    curDf$bonferroni_p = p.adjust(curDf$pvalue, method = "bonferroni")
    curDf$position = gsub("chr.{1,2}-", "", curDf$peak)
    selDf = curDf[, c('peak', 'seqnames', 'position', 'gene', 'pvalue', 'Cluster', 'fdr_p', 'bonferroni_p')]
    colnames(selDf)[2] = "Chromosome"
    outFile = paste0(targDf, sampleName)
    write.csv(selDf, outFile, row.names = F)
  }
}

#editTables(curFiles, targetDir)
