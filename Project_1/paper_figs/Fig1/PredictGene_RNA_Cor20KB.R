library(Seurat)
library(Signac)
curDate = Sys.Date()

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]

predictActivity = readRDS("PredGeneActivities_atacIntegrated_macs2_2_20kbUp_20kbDown")
atacFilt[["PredictActivity"]] = NULL

head(colnames(atacFilt))
head(colnames(predictActivity))

length(colnames(atacFilt))
length(colnames(predictActivity))

atacFilt[['PredictActivity']] <- CreateAssayObject(counts = predictActivity)
atacFilt <- NormalizeData(
  object = atacFilt,
  assay = 'PredictActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacFilt$nCount_RNA)
)

DefaultAssay(atacFilt) = "PredictActivity"
atacFilt = ScaleData(atacFilt)

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
      rna_expr = subset_seurat_object@assays$RNA@counts[curGene ,]
      atac_expr = subset_seurat_object@assays$PredictActivity@counts[curGene ,]
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

targDir = "atacRna/PredictGeneActCor/20kbUp_20kbDown/"
dir.create(targDir, recursive = T, showWarnings = F)
combDf = data.frame()
for (cluster in clusters ) {
  try({
    curDf = getCor(atacFilt=atacFilt, cluster=cluster, topN="all", allMarkersPct=allMarkersPct, corTest="spearman")
    curDf$fdr_p = p.adjust(curDf$Pval, method = "fdr")
    outFile = paste0(targDir, "PredictGene_RNA_spearman_", cluster, "_2023-10-23.csv")
    write.csv(curDf, outFile, row.names = F)
    combDf = rbind(combDf, curDf)
  })
  
}

write.csv(combDf, paste0(targDir, "PredictGene_RNA_spearman_AllClusters_2023-10-23.csv"), row.names = F)

## checks

opcDf = combDf[combDf$Cluster == "OPC",]
opcDf = opcDf[opcDf$Estimate > 0,]
opcDf = opcDf[order(opcDf$Estimate, decreasing = T),]
rownames(opcDf) = NULL
Sox6 = opcDf[opcDf$Gene == "Sox6",]
rownames(Sox6)

odcDf = combDf[combDf$Cluster == "ODC",]
odcDf = odcDf[odcDf$Estimate > 0,]
odcDf = odcDf[order(odcDf$Estimate, decreasing = T),]
rownames(odcDf) = NULL
Enpp2 = odcDf[odcDf$Gene == "Enpp2",]
rownames(Enpp2)



