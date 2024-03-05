library(EnhancedVolcano)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig3/Volcano/'

dir.create(targDir, recursive = T)

# RNA markers

groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

#EnhancedVolcano(groupMarkers,lab = groupMarkers$Genes, x = 'avg_log2FC', y = 'p_val', cutoffLineType = 'blank')

editMarkers = function(groupMarkers) {
  groupMarkers$Cluster[groupMarkers$Cell_Type == 4] = 3
  groupMarkers$Cluster[groupMarkers$Cell_Type == 1] = "ODC"
  groupMarkers$Cluster[groupMarkers$Cell_Type == 2] = "OPC"
  groupMarkers$Cluster[groupMarkers$Cell_Type == 3] = "Intermideate"
  return(groupMarkers)
}

groupMarkers= editMarkers(groupMarkers)

clusters = unique(groupMarkers$Cluster)

makePlot = function(df, logf,cluster, datType, pcut) {
    dfSub = df[(df$Cluster == cluster),]
    volcPlot = EnhancedVolcano(dfSub ,lab =  dfSub$Genes, x = 'avg_log2FC', y = 'p_val', labSize = 8, pCutoffCol = "p_val_adj", pCutoff = pcut, title = cluster, subtitle = NULL, FCcutoff = logf) + 
      theme(plot.title = element_text(size = 25))
    
    png(paste0(targDir,'Volcano_',cluster, "_", datType, "_", curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
    print(volcPlot)
    dev.off()
}

for (cluster in clusters) {
  makePlot(df = groupMarkers, logf = 1, cluster = cluster, datType = "RNA", pcut= 0.05)
}


# AUcell
rm(groupMarkers)
groupMarkers = read.csv("OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
groupMarkers= editMarkers(groupMarkers)
clusters = unique(groupMarkers$Cluster)

for (cluster in clusters) {
  print(cluster)
  makePlot(df = groupMarkers, logf = 0, cluster = cluster, datType = "AuCell_KeggClustProf", pcut = 0.001)
}

