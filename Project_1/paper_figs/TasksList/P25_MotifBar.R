library(ggplot2)

targDir = "Paper_figs/TasksList/P19/"

motDf = read.csv(paste0(targDir, "peakGeneCor_AllClusters_500K_Top10_CombinedMotifs_2023-12-22.csv"))

clusters = c("OPC", "ODC")


targDir = "Paper_figs/TasksList/P25/"
dir.create(targDir, showWarnings = F, recursive = T)
for (cluster in clusters ) {
  selClust = motDf[motDf$Cluster == cluster,]
  selClust$`-Log10Pval` = log10(selClust$pvalue) * -1
  selClust = selClust[order(selClust$`-Log10Pval`, decreasing = T),]
  plot1 = ggplot(selClust, aes(x = `-Log10Pval`, y = reorder(motif.name, `-Log10Pval`), fill = fold.enrichment)) +
    geom_bar(stat = "identity") +
    theme_classic() + 
    ylab("Motif") +
    theme(text = element_text(size = 26)) +
    scale_fill_viridis_c(option = "plasma") 
  
  plot1
  
  png(filename =  paste0(targDir, "peakGeneCor_500K_Top10_CombinedMotifsBar_", cluster,  "_2024-01-02.png"), width = 28, height = 16, units = "in", res = 300)
  print(plot1)
  dev.off()
  
}