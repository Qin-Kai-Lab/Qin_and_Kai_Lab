# another way to adjust p values would be to filter by logfc and then adjust p-values

library(ggvenn)
library(ggpubr)

setwd("/home/flyhunter/Wang/output")

curDate = Sys.Date()

list.files(pattern = ".csv")

rnagenes = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

atacGenes = read.csv("/home/flyhunter/Wang/output/archContrStress/atacRna/StressVsControl/StressVsControl_GeneScoremat_allClusters2023-03-20.csv")

nrow(atacGenes[(atacGenes$FDR_def < 0.05),])

atacGenesFilt = atacGenes[!is.na(atacGenes$FDR_def),]

nrow(atacGenesFilt[(atacGenesFilt$FDR_def < 0.05),])

nrow(atacGenes[(atacGenes$fdr_p < 0.05),])

# basic comparison

clusters = unique(rnagenes$Cell_Type)
vennGenes = function(x, y, p.colx, p.coly, p.limit) {
  for (cluster in clusters) {
    rnaSub = x[(x[['Cell_Type']] == cluster) & (x[[p.colx]] < p.limit), ]
    atacSub = y[(y[['Cluster']] == cluster) & (y[[p.coly]] < p.limit), ]
    if ((nrow(rnaSub) > 0) & (nrow(atacSub) > 0)) {
      genes = list("RNA_Genes" = unique(rnaSub$Genes), "ATAC_Genes" = unique(atacSub$Genes))
      vennPlot <-ggvenn(genes, text_size=5, set_name_size=8)+
        theme(text = element_text(size = 26)) +
        theme(plot.title = element_text(hjust = 0.5, vjust = -6))
      outFile = paste0(targDir, "Venn_", cluster, "_Plimit_", p.limit, "_",  curDate, ".jpeg")
      print(outFile)
      ggsave(file = outFile, plot = vennPlot, width = 9, height = 9, units = 'in', dpi = 300)
    }
  }
}

targDir = "atacRna/StressVsControl/RnaAdjP_AtacFdrDef/"
dir.create(targDir, recursive = T)
vennGenes(x = rnagenes, y = atacGenesFilt, p.colx = "p_val_adj", p.coly = "FDR_def", p.limit = 0.05)

# compare exact same things

edtiDf = function(x, clustCol,  log.limit, log.col, p.col) {
  results = data.frame(matrix(nrow = 0, ncol = 0))
  for (cluster in clusters) {
    xSub = x[(x[[clustCol]] == cluster),]
    xSub$AbsLog = abs(xSub[, log.col])
    xFilt  = xSub[(xSub$AbsLog > 0.25),]
    if (nrow(xFilt) > 0) {
      xFilt$fdr_p_cust = p.adjust(xFilt[, p.col], method = "fdr")
      results = rbind(results, xFilt) 
    }
  }
  return(results)
}

editRna = edtiDf(x = rnagenes, clustCol = "Cell_Type",  log.limit = 0.25, log.col = "avg_log2FC", p.col = "p_val" )

editAtac = edtiDf(x = atacGenesFilt, clustCol = "Cluster",  log.limit = 0.25, log.col = "Log2FC", p.col = "Pval" )

targDir = "atacRna/StressVsControl/RnaFdrCustP_AtacFdrCustP/"
dir.create(targDir, recursive = T)
vennGenes(x = editRna, y = editAtac, p.colx = "fdr_p_cust", p.coly = "fdr_p_cust", p.limit = 0.1)

targDir = "atacRna/StressVsControl/"
write.csv(editAtac, file = paste0(targDir, 'AtacGeneExprAllClust_editFdr_', curDate, '.csv'), row.names = F)
write.csv(editRna, file = paste0(targDir, 'RnaAllClust_editFdr_', curDate, '.csv'), row.names = F)

atacPeaks = data.frame(matrix(nrow =0, ncol =0))
filesList = list.files('atacMergeDE/adj_p/', pattern = ".csv", full.names = T)
for ( i in filesList) {
  df = read.csv(i)
  atacPeaks = rbind(atacPeaks, df)
}

atacPeaks = read.csv('archContrStress/atacRna/StressVsControl/Macs2_0.05_PerCluster/StressVsControl_Macs2PerClust_allClusters2023-03-28.csv')


colnames(atacPeaks)
editPeaks = edtiDf(x = atacPeaks, clustCol = "Cell_Type",  log.limit = 0.25, log.col = "avg_log2FC", p.col = "p_val" )
colnames(editPeaks)[9] <- "Genes"

targDir = "atacRna/StressVsControl/RnaFdrCustP_AtacFdrCustP_SignacMergeMacs/"
dir.create(targDir, recursive = T)
vennGenes(x = editRna, y = editPeaks, p.colx = "fdr_p_cust", p.coly = "fdr_p_cust", p.limit = 0.05)
