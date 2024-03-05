library(ArchR)
library(stringr)
library(gridExtra)

targDir = './Paper_figs/Fig1/'

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

custom_colors <- c("SUB" = "red", "CA1" = "blue", "CA2" = "forestgreen", 'CA3'= "purple", 'DG' = 'magenta1', 'GABA' = "orange", 'C-R' = "cyan1", 'OPC' = "brown", 'ODC' = "black", 'MG' = "gold")

# Plot with custom colors and order
plot1 <- ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = df) +
  geom_point(size = 0.5) +
  scale_color_manual(values = custom_colors, breaks = names(custom_colors)) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

# Display the plot
print(plot1)

ggsave(plot = plot1, file = paste0(targDir, 'RNA_UMAP_allClusers_CustomCol_2023-12-21.png'), width = 18, height = 12, dpi = 300, units = 'in')

write.csv(df, paste0(targDir, 'RNA_UMAP_allClusers.csv'), row.names = F)

# split UMAP
plot2 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Group), data = df) +
  geom_point(size = 0.001) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))

plot2

p3<-grid.arrange(plot1,plot2, ncol=2, nrow=1)

ggsave(plot = p3, file = paste0(targDir, 'RNA_UMAP_allClusers_Split_', curDate, '.png'), width = 22, height = 12, dpi = 300, units = 'in')

# ATAC

projFilt <- readRDS('atacArchRnaFilt')
archUmap <- getEmbedding(ArchRProj = projFilt, embedding = "UMAP_macs2", returnDF = TRUE)
archUmap$Cell_id = rownames(archUmap)
identical(archUmap$Cell_id, projFilt$cellNames)
archUmap$Annotations = projFilt$Annotations
archUmap$Group = str_to_title(projFilt$Group)
colnames(archUmap)[1:2] = c('UMAP_1', 'UMAP_2')

levels(x = archUmap) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
archUmap$Annotations <- factor(archUmap$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

custom_colors <- c("SUB" = "red", "CA1" = "blue", "CA2" = "forestgreen", 'CA3'= "purple", 'DG' = 'magenta1', 'GABA' = "orange", 'C-R' = "cyan1", 'OPC' = "brown", 'ODC' = "black", 'MG' = "gold")


# plots
plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = archUmap )+
  geom_point(size = 0.5) +
  scale_color_manual(values = custom_colors, breaks = names(custom_colors)) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))


plot1

png(filename = paste0(targDir, 'ATAC_UMAP_ArchRMacs2_allClusers_CustomCol_2023-12-21.png'), width = 18, height = 12, res = 300, units = 'in')
print(plot1)
dev.off()

# merged UMAP

intObj = readRDS("atacIntegrated_macs2_2_RNA_WNN")

df<-data.frame(intObj@reductions$umap@cell.embeddings)
df$Cell_id<-row.names(df)

metadata<-data.frame(intObj@meta.data)
metadata$Cell_id<-row.names(metadata)

identical(metadata$Cell_id, df$Cell_id)

df$Group = metadata$group
df$Annotations = as.factor(metadata$Annotations)

df$Annotations <- factor(df$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

colnames(df)[1:2] = c("UMAP_1", "UMAP_2")

custom_colors <- c("SUB" = "red", "CA1" = "blue", "CA2" = "forestgreen", 'CA3'= "purple", 'DG' = 'magenta1', 'GABA' = "orange", 'C-R' = "cyan1", 'OPC' = "brown", 'ODC' = "black", 'MG' = "gold")

plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = df) +
  geom_point(size = 0.5) +
  scale_color_manual(values = custom_colors, breaks = names(custom_colors)) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))


plot1

ggsave("Paper_figs/Fig1/RNA_ATAC_Merged_WNN_UMAP_2023-12-21.png", width = 18, height = 12, dpi = 300, units = 'in')