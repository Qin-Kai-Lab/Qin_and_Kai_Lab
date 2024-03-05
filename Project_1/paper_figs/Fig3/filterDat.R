setwd("~/Wang/output")

curDate<-Sys.Date()

targDir = './Paper_figs/Fig3/'

allMarkers = read.csv("Paper_figs/Fig3/All_RNA_Markers_2023-05-08.csv")
colnames(allMarkers)

filterTable = function(df) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  dfFilt = df[(df$avg_log2FC > 0) & (df$p_val < 0.05),]
  clusters = unique(df$cluster)
  for (cluster in clusters) {
    dfSub = dfFilt[(dfFilt$cluster == cluster) & (dfFilt$p_val_adj < 0.05),]
    dfRest = dfFilt[!(dfFilt$cluster == cluster),]
    restGenes = unique(dfRest$gene)
    dfSubFilt = dfSub[!(dfSub$gene %in% restGenes),]
    combDat  = rbind(combDat, dfSubFilt)
  }
  return(combDat)
}

filtGenes = filterTable(df = allMarkers)

write.csv(filtGenes, paste0(targDir, "RNA_Markers_FiltPosUnique_", curDate, ".csv"), row.names = F)

# make sure current genes are not expressed in other clusters
allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allFilt = allMarkers[(allMarkers$avg_log2FC > 0 ) & (allMarkers$p_val < 0.05),]
clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'MG')
allFilt = allFilt[(allFilt$cluster %in% clusters),]

newGenes = filtGenes[!(filtGenes$gene %in% allFilt$gene),]

unique(newGenes$cluster)

write.csv(newGenes, paste0(targDir, "RNA_Markers_FiltPosUniqueAllClust_", curDate, ".csv"), row.names = F)