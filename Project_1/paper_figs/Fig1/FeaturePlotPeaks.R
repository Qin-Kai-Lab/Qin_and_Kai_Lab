library(Signac)
library(Seurat)
library(ggplot2)

targDir = "Paper_figs/Fig1/FeatuePlot/"

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

DefaultAssay(atacFilt) = "PredictActivity"

genes = c("Sox6", "Enpp2")

fPlot<-FeaturePlot(atacFilt, features = genes, min.cutoff = "q10")&
  theme(text = element_text(size = 22))

fPlot

ggsave(paste0(targDir,'AtacInteg_macs2_2_Sox6_Enpp2.jpeg'), plot = fPlot, height = 10, width = 16, units = 'in', dpi = 300)

refPlot = DimPlot(atacFilt) +
  theme(text = element_text(size = 22))

ggsave(paste0(targDir,'AtacInteg_macs2_2_UMAP.jpeg'), plot = refPlot, height = 10, width = 14, units = 'in', dpi = 300)


DefaultAssay(atacFilt) = "Combined_peaks"

genes = c("chr7-115666776-115667244", "chr15-54935460-54935838")

fPlot<-FeaturePlot(atacFilt, features = genes, min.cutoff = "q10")

fPlot

ggsave(paste0(targDir,'AtacInteg_macs2_2_Sox6_Enpp2_Peaks.jpeg'), plot = fPlot, height = 10, width = 16, units = 'in', dpi = 300)