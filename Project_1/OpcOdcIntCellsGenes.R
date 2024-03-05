library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

set.seed(13)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

targDir <- 'OPC_ODC/Monocle3/cells_genes/'

dir.create(targDir, recursive = T)

cds <- as.cell_data_set(RNA.combined.norm)


cds <- preprocess_cds(cds, num_dim = 85)

cds <- align_cds(cds, num_dim = 85, alignment_group = "group")
cds <- reduce_dimension(cds)
conStr<-plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE, cell_size=1)

cds <- cluster_cells(cds)
plot_cells(cds)

plot_cells(cds, color_cells_by = "cluster")

clustDat<-cds@clusters$UMAP

clustID<-clustDat$clusters

cellNames<-names(clustID)
cellClusters<-as.character(clustID)
cellClustInf<-data.frame(cbind(cellNames, cellClusters))

metadat<-data.frame(cds@colData)

identical(rownames(metadat), cellClustInf$cellNames)

metadat$monocClust <- cellClustInf$cellClusters

# get number of cells per cluster and pie plots

sumClustGr<-data.frame(table(metadat$group, metadat$monocClust))

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

ggsave(paste0(targDir, 'OPC_ODC_MonocClust85PC_Pie_', title0, '_' , curDate, '.jpeg'), height = 9, width = 12, units = 'in', dpi = 300)
 
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

ggsave(paste0(targDir, 'OPC_ODC_MonocClust85PC_Pie_', title0, '_' , curDate, '.jpeg'), height = 9, width = 12, units = 'in', dpi = 300)

# save 
write.csv(sumClustGr, paste0(targDir, 'OpcOdc_monocClust85PC_cellNumb_', curDate, '.csv'), row.names = F)

# violin plots

genes<-c('Plp1', 'Ptgds', 'Trf', 'Calcrl', 'Cspg4', 'Vcan')

identical(rownames(cds@colData), cellClustInf$cellNames)

cds@colData$Monoc_clusters <- cellClustInf$cellClusters

rowData(cds)$gene_short_name<-rownames(cds)

cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]

plot_genes_violin(cds_subset, group_cells_by="Monoc_clusters", ncol=2, label_by_short_name = T) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(text = element_text(size = 20))

ggsave(paste0(targDir, 'OPC_ODC_MonocClust85PC_ViolinP_' , curDate, '.jpeg'), height = 12, width = 16, units = 'in', dpi = 300)