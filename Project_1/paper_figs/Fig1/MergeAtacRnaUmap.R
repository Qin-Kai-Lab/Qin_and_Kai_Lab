library(Seurat)
library(dplyr)
library(Signac)
library(ggplot2)

intObj = readRDS("atacIntegrated_macs2_2_RNA")


DimPlot(intObj)

intObj[["PredictActivity"]] = NULL

DefaultAssay(intObj) <- "Combined_peaks"
intObj <- FindTopFeatures(intObj, min.cutoff = 5)
intObj <- RunTFIDF(intObj)
intObj <- RunSVD(intObj)

gc()

DefaultAssay(intObj) = "RNA"
intObj <- NormalizeData(intObj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# intObj <- RunUMAP(object = intObj, dims = 1:40)
# DimPlot(intObj, label = TRUE, repel = TRUE)


intObj <- FindMultiModalNeighbors(
  object = intObj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
intObj <- RunUMAP(
  object = intObj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

intObj$Annotations <- factor(intObj$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

DimPlot(intObj, label = FALSE, repel = TRUE, reduction = "umap")


df<-data.frame(intObj@reductions$umap@cell.embeddings)
df$Cell_id<-row.names(df)

metadata<-data.frame(intObj@meta.data)
metadata$Cell_id<-row.names(metadata)

identical(metadata$Cell_id, df$Cell_id)

df$Group = metadata$group
df$Annotations = as.factor(metadata$Annotations)

df$Annotations <- factor(df$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

colnames(df)[1:2] = c("UMAP_1", "UMAP_2")

plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = df) +
  geom_point(size = 0.001) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  guides(shape = guide_legend(override.aes = list(size = 6)))


plot1



#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

ggsave("Paper_figs/Fig1/RNA_ATAC_Merged_WNN_UMAP_2023-12-18.png", width = 18, height = 12, dpi = 300, units = 'in')

#saveRDS(intObj, "atacIntegrated_macs2_2_RNA_WNN")