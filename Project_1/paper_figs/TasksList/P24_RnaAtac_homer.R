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

clusters = c("OPC", "ODC")

targDir = "Paper_figs/TasksList/P24/RnaAtacSignPeaks/"
dir.create(targDir, recursive = T, showWarnings = F)

for (cluster in clusters) {
  curMarkers = markers[markers$Cluster == cluster,]
  curPeaks = unique(curMarkers$peak)
  #objSub = subset(atacFilt, Annotations==cluster)
  #curExpression = data.frame(AverageExpression(object=objSub,assays = "Combined_peaks", features = curPeaks, group.by = "ident", layer = "counts"))
  
  # check to see if there are peaks with 0 fragments, yes there are
  #curExpression = data.frame(AggregateExpression(object=objSub,assays = "Combined_peaks", features = curPeaks, group.by = "ident"))
  
  combDat = data.frame() 
  for ( i in curPeaks ) {
    curList = strsplit(i, "-")[[1]]
    curDf = t(data.frame(curList))
    colnames(curDf) = c("chromosome", "start", "end")
    combDat = rbind(combDat, curDf)
  }
  
  rownames(combDat) = NULL
  combDat$ID = row.names(combDat)
  combDat$V5 = NA
  combDat$Strand = "+"
  outBed = paste0(targDir, "Peaks_", cluster, ".bed")
  write.table(combDat, outBed, col.names = T, row.names = F, quote = F, sep = "\t")
  
  # need active homer conda environment
  outDir = paste0(targDir, "Motif_", cluster, "/")
  cmd_str = paste0("findMotifsGenome.pl ", outBed, " mm10 ", outDir, " -size 200 -mask")
  system(cmd_str)
}


