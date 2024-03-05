library(ggplot2)

df = read.csv("Paper_figs/Fig3/All_AUCell_KeggClustProf_Markers_2023-05-09.csv")
clusters = unique(df$cluster)

colnames(df)

colnames(df)[1] = 'pvalue'
colnames(df)[6:7] = c("Cluster", "Description")

# need to rename the columns to match the program, rename final paths and filter based on positvie negative logf2 change at the beginning of the code
# stopped here

getTopPath = function(groupMarkers, topNumb) {
  combDf = data.frame(matrix(nrow = 0, ncol = 0))
  for (curCluster in clusters) {
    curDf = groupMarkers[(groupMarkers$Cluster == curCluster),]
    curDfOrd = curDf[order(curDf$pvalue, decreasing = F),]
    dfTop = head(curDfOrd, topNumb)
    combDf = rbind(combDf, dfTop)
    topMark = unique(combDf$Description)
  }
  return(topMark)
}

topMark = getTopPath(groupMarkers = df, topNumb = 5)

dfFilt = df[df$Description%in%topMark,]

dfFilt$`-LogPadj` = log10(dfFilt$p_val_adj) * -1


dfFilt$Cluster = factor(dfFilt$Cluster, levels =  c('OPC', 'Intermideate', 'ODC'))

plot1 = ggplot(dfFilt, aes(x = Description, y = `-LogPadj`, fill = avg_log2FC)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(text = element_text(size = 26)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

plot2 = plot1  + facet_grid(cols = vars(Cluster))

targDir = "Paper_figs/Fig3/"
topN = "5"
typeGenes = "All"
png(filename = paste0(targDir,"AuCell_Top_", topN, "_", typeGenes, "Box_KEGG_2023-12-14.png"), width = 28, height = 16, units = "in", res = 300)
print(plot2)
dev.off()


##
df = read.csv("Paper_figs/Fig3/All_AUCell_KeggClustProf_Markers_2023-05-09.csv")
#df = df[df$avg_log2FC > 0,]
df = df[df$avg_log2FC < 0,]
clusters = unique(df$cluster)

colnames(df)[1] = 'pvalue'
colnames(df)[6:7] = c("Cluster", "Description")

getTopPath = function(groupMarkers, topNumb) {
  combDf = data.frame(matrix(nrow = 0, ncol = 0))
  for (curCluster in clusters) {
    curDf = groupMarkers[(groupMarkers$Cluster == curCluster),]
    curDfOrd = curDf[order(curDf$pvalue, decreasing = F),]
    dfTop = head(curDfOrd, topNumb)
    combDf = rbind(combDf, dfTop)
    topMark = unique(combDf$Description)
  }
  return(combDf)
}

dfFilt = getTopPath(df, 5)

dfFilt$`-LogPval` = log10(dfFilt$pvalue) * -1

dfFilt = dfFilt[order(dfFilt$Cluster, dfFilt$`-LogPval`, decreasing = F),]

dfFilt$Cluster = factor(dfFilt$Cluster, levels =  c('OPC', 'Intermideate', 'ODC'))
dfFilt$Description = factor(dfFilt$Description, levels =  unique(dfFilt$Description))


plot1 <- ggplot(dfFilt, aes(x = `-LogPval`, y  = Description, fill = avg_log2FC)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(text = element_text(size = 26)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) &
  facet_grid(rows = vars(Cluster), scales = "free_y")

plot2 <- plot1 + facet_grid(rows = vars(Cluster), scales = "free_y")

targDir = "Paper_figs/Fig3/"
topN = "5"
typeGenes = "Negative"
png(filename = paste0(targDir,"AuCell_Only_", topN, "_", typeGenes, "Bar_KEGG_2023-12-14.png"), width = 28, height = 16, units = "in", res = 300)
print(plot2)
dev.off()


