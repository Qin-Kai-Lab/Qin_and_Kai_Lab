
library(Seurat)
library(ggplot2)
library(dplyr)


set.seed(13)

curDate<-Sys.Date()

setwd("~/Wang/output")

targDir = './Paper_figs/Fig1/Cell_Numbers/'
dir.create(targDir, recursive = T, showWarnings = F)

# get RNA seq UMAP
source('../programs/renameClusters.R')

curDate = Sys.Date()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

metadat <- RNA.combined.norm@meta.data

# pie chart

sumClustGr<-data.frame(table(metadat$Annotations))
colnames(sumClustGr)[1]<-c('Cluster')
sumClustGr$Cluster = factor(sumClustGr$Cluster, levels = c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

sumClustGr <- sumClustGr %>%
  arrange(desc(Cluster)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)

sumClustGr$percent<-round(sumClustGr$Freq/sum(sumClustGr$Freq)*100, digits = 1)
sumClustGr$Percent<-paste0(sumClustGr$percent, '%', ' (', sumClustGr$Freq,')')


#write.csv(sumClustGr, "RNA_CellsPerClust.csv", row.names = F)

bPlot = ggplot(sumClustGr, aes(x = Cluster, y= Freq, fill = Cluster)) + 
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Number of Cells") +
  theme_classic() +
  theme(
    text = element_text(size = 16),  # Increase font size of axis labels
    plot.title = element_text(size = 16, face = "bold")  # Customize title font
  ) + 
  geom_text(aes(label = Percent), vjust = -0.5, size = 4)

bPlot

png(file = paste0(targDir, "RNA_CellsPerClust_Bar.png"), height = 10, width = 12, res = 300, units = "in")
print(bPlot)
dev.off()