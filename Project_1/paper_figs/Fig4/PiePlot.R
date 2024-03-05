library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

set.seed(13)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

targDir = './Paper_figs/Fig4/'
dir.create(targDir)

# get number of cells per cluster and pie plots

metadat <- RNA.combined.norm@meta.data

renameClusters = function(df, clustCol) {
  newClust = character()
  for ( i in 1:nrow(df)) {
    if (df[i, clustCol] == 1) {
      curClust = "ODC"
    } else if (df[i, clustCol] == 2) {
      curClust = "OPC"
    } else if ((df[i, clustCol] == 3) | (df[i, clustCol] == 4)) {
      curClust = "Intermideate"
    }
    newClust = c(newClust, curClust)
  }
  df$newMonocClust = newClust
  return(df)
}

metadat = renameClusters(df = metadat, clustCol = "MonocClust")

sumClustGr<-data.frame(table(metadat$group, metadat$newMonocClust))
colnames(sumClustGr)[1:2]<-c('Group', 'monoc_cluster')
sumClustGr$monoc_cluster = factor(sumClustGr$monoc_cluster, levels = c("OPC", "Intermideate", "ODC"))

makePlots = function(df, varCol) {
  mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CC79A7")
  
  vars = unique(df[, varCol])
  for (var in vars) {
    oneGr = df[(df[varCol] == var),]
    oneGr <- oneGr %>%
      arrange(desc(monoc_cluster)) %>%
      mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)
    title0<- var
    oneGr$percent<-round(oneGr$Freq/sum(oneGr$Freq)*100, digits = 1)
    oneGr$Percent<-paste0(oneGr$percent, '%', ' (', oneGr$Freq,')')
    piePlot = ggplot(oneGr, aes(x = "", y = Freq, fill = monoc_cluster)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0)+
      #geom_text(aes(y = lab.ypos, label = Percent), color = "black")+
      scale_fill_manual(values = mycols) +
      guides(fill = guide_legend(title = "Cluster")) +
      theme_void()+
      ggtitle(title0)+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(text = element_text(size = 20))   
    
    png(paste0(targDir, "PiePlot_3Groups_", title0, "_2023-12-13.png"), height = 9, width = 12, units = 'in', res = 300)
    print(piePlot)
    dev.off()
  }
                
}

makePlots(df = sumClustGr, varCol = "Group")