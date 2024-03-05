library(ggplot2)
library(gridExtra)
library(viridis)

allMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
allMarkers = allMarkers[allMarkers$p_val_adj < 0.05,]
clusters = unique(allMarkers$Cell_Type)

makeHeatMaps = function(cluster, allMarkers) {
  #curMarkersDf = allMarkers[allMarkers$Cell_Type == cluster,]
  curMarkersDf = allMarkers[allMarkers$Cell_Type == cluster,]
  curMarkersDf$Genes <- factor(curMarkersDf$Genes, levels = curMarkersDf$Genes[order(curMarkersDf$avg_log2FC, decreasing = TRUE)])
  curHeatMap = ggplot(curMarkersDf, aes(x=Genes, y=Cell_Type, fill=avg_log2FC))+
    geom_tile() +
    theme_classic() +
    scale_fill_viridis(option="magma", discrete=F) +
    theme(text = element_text(size = 24), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 20)) +
    theme(axis.text.y = element_text(size = 20)) + 
    theme(axis.text.x = element_text(size = 3)) +
    #labs(title = cluster) +
    scale_y_discrete(expand=c(0, 0)) 
  
  return(curHeatMap)
}

clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'C-R', 'SUB')
clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'SUB')

combPlots = list()
for (i in 1:length(clusters)) {
  cluster = clusters[i]
  curPlot = makeHeatMaps(cluster=cluster, allMarkers=allMarkers)
  combPlots[[i]] = curPlot
}

grid.arrange(grobs = combPlots, ncol = 5) 

targDir = './Paper_figs/Fig2/HeatMap/'

jpeg(filename = paste0(targDir, "AllClusters_Heatmap_SignGenes_log2Fc_2023-12-05.jpeg"), width = 46, height = 20, units = "in", res = 300)
print(grid.arrange(grobs = combPlots, ncol = 5))
dev.off()

### all clusters together

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Genes)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  allPAdj = numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
    for (i in keggList){
      if (!(i %in% dfSel$Genes)){
        curDesc<-i
        curClust<-cluster
        curP<-1
        curCount<-0
        curPAdj = 1
        allDesc<-c(allDesc, curDesc)
        allClust<-c(allClust, curClust)
        allP<-c(allP, curP)
        allCount<-c(allCount, curCount)
        allPAdj = c(allPAdj, curPAdj)
        # combine all vectors and make a dataframe from them
      }
    }
  }
  finalTable<-data.frame(Genes=allDesc, p_val=allP, p_val_adj=allPAdj, avg_log2FC=allCount, Cell_Type=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab = allMarkers, clusters = clusters)
keggComplete<-rbind.fill(allMarkers, missDat)


keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =  c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'C-R', 'SUB'))
curHeatMap = ggplot(keggComplete, aes(x=Genes, y=Cell_Type, fill=avg_log2FC))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option="magma", discrete=F) +
  theme(text = element_text(size = 24), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) + 
  theme(axis.text.x = element_text(size = 1)) +
  #labs(title = cluster) +
  scale_y_discrete(expand=c(0, 0)) 

jpeg(filename = paste0(targDir, "AllClustersComb_Heatmap_SignifGenes_log2Fc_2023-12-05.jpeg"), width = 46, height = 20, units = "in", res = 300)
print(curHeatMap)
dev.off()
