library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

set.seed(13)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

targDir <- 'OPC_ODC/Monocle3/86PC/cells_genes/'

# get number of cells per cluster and pie plots

metadat <- RNA.combined.norm@meta.data

sumClustGr<-data.frame(table(metadat$group, metadat$MonocClust))

colnames(sumClustGr)[1:2]<-c('Group', 'monoc_cluster')

# pie chart control

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#CC79A7")

oneGr<-sumClustGr[(sumClustGr$Group == 'Control'),]

oneGr <- oneGr %>%
  arrange(desc(monoc_cluster)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)

title0<- unique(oneGr$Group)

oneGr$percent<-round(oneGr$Freq/sum(oneGr$Freq)*100, digits = 1)

oneGr$Percent<-paste0(oneGr$percent, '%', ' (', oneGr$Freq,')')


ggplot(oneGr, aes(x = "", y = Freq, fill = monoc_cluster)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Percent), color = "black")+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle(title0)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 20))   

ggsave(paste0(targDir, 'OPC_ODC_MonocClust86PC_Pie_', title0, '_' , curDate, '.jpeg'), height = 9, width = 12, units = 'in', dpi = 300)

# pie chart stress

oneGr<-sumClustGr[(sumClustGr$Group == 'Stress'),]

oneGr <- oneGr %>%
  arrange(desc(monoc_cluster)) %>%
  mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)

title0<- unique(oneGr$Group)

oneGr$percent<-round(oneGr$Freq/sum(oneGr$Freq)*100, digits = 1)

oneGr$Percent<-paste0(oneGr$percent, '%', ' (', oneGr$Freq,')')


ggplot(oneGr, aes(x = "", y = Freq, fill = monoc_cluster)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = Percent), color = "black")+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle(title0)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 20))

ggsave(paste0(targDir, 'OPC_ODC_MonocClust86PC_Pie_', title0, '_' , curDate, '.jpeg'), height = 9, width = 12, units = 'in', dpi = 300)

# save 
write.csv(sumClustGr, paste0(targDir, 'OpcOdc_monocClust86PC_cellNumb_', curDate, '.csv'), row.names = F)