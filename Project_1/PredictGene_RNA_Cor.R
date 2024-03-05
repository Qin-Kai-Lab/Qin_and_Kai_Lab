library(Seurat)
library(Signac)
curDate = Sys.Date()

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]

getCor = function(atacFilt, cluster, topN, allMarkersPct, corTest) {
  Estimate = numeric()
  Pval = numeric()
  Gene = character()
  markSub = allMarkersPct[allMarkersPct$cluster == cluster,]
  markSub = markSub[markSub$p_val_adj < 0.05 & markSub$avg_log2FC > 0,]
  if (topN == "all") {
    genesList = markSub$gene
  } else {
    genesList = head(markSub$gene, topN)
  }
  subset_seurat_object <- subset(atacFilt, Annotations == cluster)
  for (curGene in genesList ) {
    try({
      rna_expr = subset_seurat_object@assays$RNA@data[curGene ,]
      atac_expr = subset_seurat_object@assays$PredictActivity@data[curGene ,]
      curCor = cor.test(rna_expr, atac_expr,  method = corTest)
      curCoef = curCor$estimate
      curPval = curCor$p.value
      # combine
      Gene = c(Gene, curGene)
      Estimate = c(Estimate, curCoef)
      Pval = c(Pval, curPval)
    })

  }
  
  combDf = data.frame(Gene, Estimate, Pval)
  combDf$Cluster = cluster
  return(combDf)

}

clusters = unique(atacFilt$Annotations)

targDir = "atacRna/PredictGeneActCor/"
dir.create(targDir, recursive = T, showWarnings = F)
combDf = data.frame()
for (cluster in clusters ) {
  try({
    curDf = getCor(atacFilt=atacFilt, cluster=cluster, topN="all", allMarkersPct=allMarkersPct, corTest="spearman")
    curDf$fdr_p = p.adjust(curDf$Pval, method = "fdr")
    outFile = paste0(targDir, "PredictGene_RNA_spearman_", cluster, "_2023-10-10.csv")
    write.csv(curDf, outFile, row.names = F)
    combDf = rbind(combDf, curDf)
  })

}

write.csv(combDf, paste0(targDir, "PredictGene_RNA_spearman_AllClusters_2023-10-10.csv"), row.names = F)

