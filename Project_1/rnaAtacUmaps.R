library(ArchR)
library(stringr)
library(gridExtra)

targDir = './RnaAtacPlots/'

dir.create(targDir)

# get RNA seq UMAP
source('../programs/renameClusters.R')

curDate = Sys.Date()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

df<-data.frame(RNA.combined.norm@reductions$umap@cell.embeddings)
df$Cell_id<-row.names(df)

metadata<-data.frame(RNA.combined.norm@meta.data)
metadata$Cell_id<-row.names(metadata)

identical(metadata$Cell_id, df$Cell_id)

df$Group = metadata$group
df$Annotations = as.factor(metadata$Annotations)

# atac

projFilt <- readRDS('atacArchRnaFilt')
archUmap <- getEmbedding(ArchRProj = projFilt, embedding = "UMAP", returnDF = TRUE)
archUmap$Cell_id = rownames(archUmap)
identical(archUmap$Cell_id, projFilt$cellNames)
archUmap$Annotations = projFilt$Annotations
archUmap$Group = str_to_title(projFilt$Group)
colnames(archUmap)[1:2] = c('UMAP_1', 'UMAP_2')

levels(x = archUmap) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
archUmap$Annotations <- factor(archUmap$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
# plots
plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = archUmap )+
  geom_point(size = 1) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))
  

plot1

plot2 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Group), data = archUmap )+
  geom_point(size = 1) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

plot2

p3<-grid.arrange(plot1,plot2, ncol=2, nrow=1)

ggsave(plot = p3, file = paste0(targDir, 'AtacUmapTileMat_allClusers_', curDate, '.png'), width = 22, height = 12, dpi = 300, units = 'in')

# make plots with ATAC peaks expression

archUmap <- getEmbedding(ArchRProj = projFilt, embedding = "UMAP_macs2", returnDF = TRUE)
archUmap$Cell_id = rownames(archUmap)
identical(archUmap$Cell_id, projFilt$cellNames)
archUmap$Annotations = projFilt$Annotations
archUmap$Group = str_to_title(projFilt$Group)
colnames(archUmap)[1:2] = c('UMAP_1', 'UMAP_2')

levels(x = archUmap) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
archUmap$Annotations <- factor(archUmap$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
# plots
plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = archUmap )+
  geom_point(size = 0.04) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))


plot1

plot2 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Group), data = archUmap )+
  geom_point(size = 0.04) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

plot2

p3<-grid.arrange(plot1,plot2, ncol=2, nrow=1)

ggsave(plot = p3, file = paste0(targDir, 'AtacUmapMacs2_allClusers_', curDate, '.png'), width = 22, height = 12, dpi = 300, units = 'in')

# make Umaps with Rna seq 
#levels(x = df) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG') did not work
df$Annotations <- factor(df$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = df) +
  geom_point(size = 0.000) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))


plot1

plot2 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Group), data = df) +
  geom_point(size = 0.001) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

plot2

p3<-grid.arrange(plot1,plot2, ncol=2, nrow=1)

ggsave(plot = p3, file = paste0(targDir, 'RNA_UMAP_allClusers_', curDate, '.png'), width = 22, height = 12, dpi = 300, units = 'in')