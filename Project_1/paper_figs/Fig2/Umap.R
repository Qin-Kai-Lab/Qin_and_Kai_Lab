library(ArchR)
library(stringr)
library(gridExtra)

targDir = './Paper_figs/Fig2/'

dir.create(targDir, recursive = T)

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

df$Annotations <- factor(df$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# split UMAP
plot2 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Group), data = df) +
  geom_point(size = 0.001) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

plot2

ggsave(plot = plot2, file = paste0(targDir, 'RNA_UMAP_StrContr_', curDate, '.png'), width = 18, height = 12, dpi = 300, units = 'in')