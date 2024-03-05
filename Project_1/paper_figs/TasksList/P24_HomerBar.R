library(ggplot2)

targDir = "Paper_figs/TasksList/P19/"

odc = read.delim("~/Wang/output/Paper_figs/TasksList/P24/RnaAtacSignPeaks/Motif_ODC/knownResults.txt")
opc = read.delim("~/Wang/output/Paper_figs/TasksList/P24/RnaAtacSignPeaks/Motif_OPC/knownResults.txt")

targDir = "Paper_figs/TasksList/P24/"

  
  makePlot = function(df, targDir, cluster) {
    df$`-Log10Pval` = log10(df$P.value) * -1
    selClust = head(df, 10)
    selClust = selClust[order(selClust$`-Log10Pval`, decreasing = T),]
    selClust$Motif.Name = gsub("\\(.*", "",  selClust$Motif.Name)
    selClust$PercTarg = as.numeric(as.character(gsub("%", "", selClust$X..of.Target.Sequences.with.Motif)))
    selClust$PercBack = as.numeric(as.character(gsub("%", "", selClust$X..of.Background.Sequences.with.Motif)))
    selClust$fold.percent = selClust$PercTarg /  selClust$PercBack
    
    plot1 = ggplot(selClust, aes(x = `-Log10Pval`, y = reorder(Motif.Name, `-Log10Pval`), fill = fold.percent)) +
      geom_bar(stat = "identity") +
      theme_classic() + 
      ylab("Motif") +
      theme(text = element_text(size = 26)) +
      scale_fill_viridis_c(option = "plasma") 
    
    plot1
    
    png(filename =  paste0(targDir, "peakGeneCor_500K_Top10_HomerMotifsBar_", cluster,  "_2024-01-16.png"), width = 28, height = 16, units = "in", res = 300)
    print(plot1)
    dev.off()
    
  }
  
  makePlot(df = odc, targDir = "Paper_figs/TasksList/P24/", cluster = "ODC")
  makePlot(df = opc, targDir = "Paper_figs/TasksList/P24/", cluster = "OPC")