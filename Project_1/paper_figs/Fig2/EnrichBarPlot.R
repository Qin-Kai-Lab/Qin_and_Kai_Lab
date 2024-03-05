library(ggplot2)

allGenes = read.csv("/home/flyhunter/Wang/output/Paper_figs/Fig2/AllClust_AllGenes_2023-06-14_p_val_adj.csv")

allGenes$KEGG = gsub(" - Mus musculus \\(house mouse\\)", "", allGenes$Description)

clusters = c("DG", "CA1")

allGenes$`-Log10(Pval_Adj)` = log10(allGenes$p.adjust) * -1

targDir = './Paper_figs/Fig2/BarPlot/'
dir.create(targDir, recursive = T)

plotPerClust = function(groupMarkers, topNumb, clusters) {
  for (curCluster in clusters) {
    print(curCluster)
    curDf = groupMarkers[(groupMarkers$Cluster == curCluster),]
    curDfOrd = curDf[order(curDf$p.adjust, decreasing = F),]
    dfTop = head(curDfOrd, topNumb)
    
    curPlot = ggplot(dfTop, aes(x = `-Log10(Pval_Adj)`, y = reorder(KEGG, `-Log10(Pval_Adj)`))) +
      geom_bar(stat = "identity", fill = "#F8766D") +
      labs(x = "-Log10(Pval_Adj)", y = "KEGG Pathway") +
      theme_minimal() +
      theme(text = element_text(colour = "black", size = 18)) + 
      theme(axis.text.x = element_text(size =18, color = "black")) +
      theme(axis.text.y = element_text(size =18, color = "black"))
    curPlot
    outFile = paste0(targDir, "ClustProf_Top", topNumb, "_KEGG_",  curCluster, "_2023-10-06.png")
    print(outFile)
    png(filename = outFile, height = 10, width = 16, units = "in", res = 300)
    print(curPlot)
    dev.off()
    
  }
}

plotPerClust(groupMarkers=allGenes, topNumb=10, clusters=clusters)
