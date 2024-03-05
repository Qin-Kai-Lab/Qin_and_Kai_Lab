library(Seurat)
library(monocle3)
cds<-readRDS('OpcOdcInt_MonocClust_86PC')

source('../programs/renameClusters.R')

targDir = './Paper_figs/Fig1/FeatuePlot/Custom/'

dir.create(targDir, recursive = T, showWarnings = F)


# plot function

x = c('Gpr17', 'Cxcr4', 'Cxcr2', 'Cyp46a1', "Ackr3")


allGenes = rownames(RNA.combined.norm)

x[!x%in%allGenes]
allGenes[(grepl("CXCR7", allGenes, ignore.case = T))]

fPlot<-FeaturePlot(RNA.combined.norm, features = x, min.cutoff = "q9")

ggsave(paste0(targDir,'FeaturePlot_2024-01-17.png'), plot = fPlot, height = 14, width = 16, units = 'in', dpi = 300)

# opc/PDC

opcOdc = subset(RNA.combined.norm, Annotations == "OPC" | Annotations == "ODC")

fPlot<-FeaturePlot(opcOdc, features = x, min.cutoff = "q9")

ggsave(paste0(targDir,'FeaturePlot_OPC_ODC_2024-01-19.png'), plot = fPlot, height = 14, width = 16, units = 'in', dpi = 300)


uMap_obj = cds@reduce_dim_aux@listData[["UMAP"]]@listData[["model"]]@listData[["umap_model"]][["embedding"]]

opcOdc[['UMAP']] <- CreateDimReducObject(embeddings = uMap_obj, key = "UMAP_", global = T, assay = "RNA")
fPlot<-FeaturePlot(opcOdc, features = x, min.cutoff = "q9")


fPlot<-FeaturePlot(cds, features = x, min.cutoff = "q9")
ggsave(paste0(targDir,'FeaturePlot_OPC_ODC_2024-01-19.png'), plot = fPlot, height = 14, width = 16, units = 'in', dpi = 300)