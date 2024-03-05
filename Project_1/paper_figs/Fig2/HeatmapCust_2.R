source('../programs/renameClusters.R')

library(ggplot2)
library(gridExtra)

allMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
allMarkers = allMarkers[allMarkers$p_val_adj < 0.05,]

#clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'C-R', 'SUB')
clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'SUB')


# make expression based on the difference with total expression

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
  print(max(combDat$Adj_Expression))
  print(min(combDat$Adj_Expression))
  combDat$Genes = as.factor(combDat$Genes)
  combDat$Genes <- factor(combDat$Genes, levels = curMarkersDf$Genes[order(curMarkersDf$avg_log2FC, decreasing = TRUE)])
  curHeatMap = ggplot(combDat, aes(y=Genes, x=Group, fill=Adj_Expression)) +
    geom_tile() +
    theme_classic() +
    scale_fill_viridis_c(option = "plasma") +  # Set fill range from -1 to 1 with breaks in increments of 0.2
    theme(text = element_text(size = 28), plot.title = element_text(hjust = 0.5))+
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title = element_text(size = 20)) +
    theme(axis.text.y = element_text(size = 3)) + 
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
    scale_x_discrete(expand=c(0, 0)) +
    scale_y_discrete(expand=c(0, 0)) &
    theme(axis.title.x = element_blank()) &
    theme(axis.title.y = element_blank()) &
    theme(
      panel.border = element_blank(),
      panel.spacing = unit(0, "lines"),
      plot.margin = unit(c(0.2, 0, 0.2, 0.1), "lines")
      #plot.margin = unit(c(2, 2, 2, 2), "lines")
    )
  
  #+
  # guides(
  #   fill = guide_legend(override.aes = list(alpha = 0))
  # ) +
  # theme(
  #   legend.background = element_rect(fill = "transparent"),
  #   legend.box.background = element_blank(),
  #   legend.title = element_text(color = "transparent"),
  #   legend.text = element_text(color = "transparent"),
  #   legend.key = element_blank()
  # ) 
  
  nPl = curHeatMap +labs(tag = cluster) +
    theme(plot.tag.position = c(0.85, 0.92))
  return(nPl)
}

#clusters = c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'SUB')
clusters = c('CA1', 'ODC', 'CA2', 'OPC', 'CA3', 'MG', 'DG', 'SUB', 'GABA')

combPlots = list()
for (i in 1:length(clusters)) {
  cluster = clusters[i]
  curPlot = makeHeatMaps(RNA.combined.norm=RNA.combined.norm, cluster=cluster, allMarkers=allMarkers)
  combPlots[[i]] = curPlot
}

grid.arrange(grobs = combPlots, ncol = 2) 

targDir = './Paper_figs/Fig2/HeatMap/'

jpeg(filename =  paste0(targDir, "AllClusters_Heatmap_SignGenes_DiffAdj_2023-12-08.jpeg"), width = 20, height = 40, units = "in", res = 600)
print(grid.arrange(grobs = combPlots, ncol = 2), widths = c(1.1, 1.1))
dev.off()