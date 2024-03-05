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

getCor = function(atacFilt, group, topN, allMarkersPct, corTest) {
  Estimate = numeric()
  Pval = numeric()
  Gene = character()
  markSub = allMarkersPct
  markSub = markSub[markSub$p_val_adj < 0.05 & markSub$avg_log2FC > 0,]
  if (topN == "all") {
    genesList = unique(markSub$gene)
  } else {
    genesList = unique(head(markSub$gene, topN))
  }
  subset_seurat_object <- subset(atacFilt, dataset == group)
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
  combDf$Group = group
  return(combDf)
  
}

groups = unique(atacFilt$dataset)

targDir = "Paper_figs/TasksList/P26/"
dir.create(targDir, recursive = T, showWarnings = F)
combDf = data.frame()
for (group in groups ) {
  try({
    curDf = getCor(atacFilt=atacFilt, group=group, topN="all", allMarkersPct=allMarkersPct, corTest="spearman")
    curDf$fdr_p = p.adjust(curDf$Pval, method = "fdr")
    outFile = paste0(targDir, "PredictGene_20KB_RNA_spearman_", group, "_2024-01-03.csv")
    write.csv(curDf, outFile, row.names = F)
    combDf = rbind(combDf, curDf)
  })
  
}

write.csv(combDf, paste0(targDir, "PredictGene_20KB_RNA_spearman_AllClusters_SepGroups_2024-01-03.csv"), row.names = F)