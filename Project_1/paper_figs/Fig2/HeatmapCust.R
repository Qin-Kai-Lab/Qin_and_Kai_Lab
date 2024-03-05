source('../programs/renameClusters.R')

library(ggplot2)
library(gridExtra)

allMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
allMarkers = allMarkers[allMarkers$p_val_adj < 0.05,]

makeHeatMaps = function(RNA.combined.norm, cluster, allMarkers) {
  curMarkers = unique(allMarkers$Genes[allMarkers$Cell_Type == cluster])
  curMarkersDf = allMarkers[allMarkers$Cell_Type == cluster,]
  objSub = subset(RNA.combined.norm, Annotations == cluster)
  Idents(objSub) = objSub$group
  curExpression = data.frame(AverageExpression(object=objSub,assays = "RNA", features = curMarkers, group.by = "ident", layer = "scale.data"))
  curExpression$Genes = rownames(curExpression)
  exprDat = data.frame(curExpression$RNA.Control)
  exprDat$Group = "Control"
  stressDat = data.frame(curExpression$RNA.Stress)
  stressDat$Group = "Stress"
  exprDat$Genes = rownames(curExpression)
  stressDat$Genes = rownames(curExpression)
  colnames(exprDat)[1] = "Expression"
  colnames(stressDat)[1] = "Expression"
  combDat = rbind(exprDat, stressDat)
  combDat$Genes = as.factor(combDat$Genes)
  combDat$Genes <- factor(combDat$Genes, levels = curMarkersDf$Genes[order(curMarkersDf$avg_log2FC, decreasing = TRUE)])
  curHeatMap = ggplot(combDat, aes(y=Group, x=Genes, fill=Expression))+
    geom_tile() +
    theme_classic() +
    scale_fill_viridis_c(option = "plasma") +
    theme(text = element_text(size = 24), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 20)) +
    theme(axis.text.y = element_text(size = 20)) + 
    theme(axis.text.x = element_text(size = 3)) +
    labs(title = cluster) +
    scale_y_discrete(expand=c(0, 0)) 
  
  return(curHeatMap)
}

#clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'C-R', 'SUB')
clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'SUB')

combPlots = list()
for (i in 1:length(clusters)) {
  cluster = clusters[i]
  curPlot = makeHeatMaps(RNA.combined.norm=RNA.combined.norm, cluster=cluster, allMarkers=allMarkers)
  combPlots[[i]] = curPlot
}

grid.arrange(grobs = combPlots, ncol = 5) 

targDir = './Paper_figs/Fig2/HeatMap/'

jpeg(filename = paste0(targDir, "AllClusters_Heatmap_SignGenes_2023-12-05.jpeg"), width = 46, height = 20, units = "in", res = 300)
print(grid.arrange(grobs = combPlots, ncol = 5))
dev.off()


# make expression based on the difference with total expression

allMarkers = allMarkers[allMarkers$p_val_adj < 0.05,]

makeHeatMaps = function(RNA.combined.norm, cluster, allMarkers) {
  curMarkers = unique(allMarkers$Genes[allMarkers$Cell_Type == cluster])
  curMarkersDf = allMarkers[allMarkers$Cell_Type == cluster,]
  objSub = subset(RNA.combined.norm, Annotations == cluster)
  combExpr = data.frame(AverageExpression(object=objSub,assays = "RNA", features = curMarkers, group.by = "ident", layer = "scale.data"))
  combExpr$Genes = rownames(combExpr)
  combExpr$Group = "All"
  colnames(combExpr)[1] = "Expression"
  Idents(objSub) = objSub$group
  curExpression = data.frame(AverageExpression(object=objSub,assays = "RNA", features = curMarkers, group.by = "ident", layer = "scale.data"))
  curExpression$Genes = rownames(curExpression)
  exprDat = data.frame(curExpression$RNA.Control)
  exprDat$Group = "Control"
  stressDat = data.frame(curExpression$RNA.Stress)
  stressDat$Group = "Stress"
  exprDat$Genes = rownames(curExpression)
  stressDat$Genes = rownames(curExpression)
  colnames(exprDat)[1] = "Expression"
  colnames(stressDat)[1] = "Expression"
  identical(combExpr$Genes, exprDat$Genes)
  identical(combExpr$Genes, stressDat$Genes)
  # adjust  by total expression
  exprDat$Adj_Expression = exprDat$Expression - combExpr$Expression
  stressDat$Adj_Expression = stressDat$Expression - combExpr$Expression
  combDat = rbind(exprDat, stressDat)
  combDat$Genes = as.factor(combDat$Genes)
  combDat$Genes <- factor(combDat$Genes, levels = curMarkersDf$Genes[order(curMarkersDf$avg_log2FC, decreasing = TRUE)])
  curHeatMap = ggplot(combDat, aes(y=Group, x=Genes, fill=Adj_Expression))+
    geom_tile() +
    theme_classic() +
    scale_fill_viridis_c(option = "plasma") +
    theme(text = element_text(size = 24), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 20)) +
    theme(axis.text.y = element_text(size = 20)) + 
    theme(axis.text.x = element_text(size = 3)) +
    labs(title = cluster) +
    scale_y_discrete(expand=c(0, 0)) 
  
  return(curHeatMap)
}

clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'SUB')

combPlots = list()
for (i in 1:length(clusters)) {
  cluster = clusters[i]
  curPlot = makeHeatMaps(RNA.combined.norm=RNA.combined.norm, cluster=cluster, allMarkers=allMarkers)
  combPlots[[i]] = curPlot
}

grid.arrange(grobs = combPlots, ncol = 5) 

targDir = './Paper_figs/Fig2/HeatMap/'

jpeg(filename = paste0(targDir, "AllClusters_Heatmap_SignGenes_DiffAdj_2023-12-07.jpeg"), width = 46, height = 20, units = "in", res = 300)
print(grid.arrange(grobs = combPlots, ncol = 5))
dev.off()




