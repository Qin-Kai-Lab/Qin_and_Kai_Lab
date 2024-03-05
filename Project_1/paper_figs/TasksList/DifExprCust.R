library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)

editClusterNames = function(meta) {
  meta$newMonocClust = as.character(meta$newMonocClust)
  meta$Cluster = meta$newMonocClust
  for ( i in 1:nrow(meta) ) {
    if ( meta$Cell_ID[i] %in%opc_inter  ) {
      meta$Cluster[i] = "OPC_Inter"
    } else if (meta$Cell_ID[i]%in%odc) {
      meta$Cluster[i] = "ODC"
    } else {
      meta$Cluster[i] = meta$newMonocClust[i]
    }
  }
  return(meta)
}

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

opcOdc =  readRDS('integ_OPC_ODC')
opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

opcOdc$Cell_ID = rownames(opcOdc@meta.data)

atacFiltSub = subset(atacFilt, Merged_CellName%in%rownames(opcOdc@meta.data))
opcOdcSub = subset(opcOdc, Cell_ID%in%rownames(atacFiltSub@meta.data))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))


opcOdcSub[["PredictActivity"]] = atacFiltSub[["PredictActivity"]]

DefaultAssay(opcOdcSub) = "PredictActivity"

opc_inter = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "OPC" | opcOdcSub@meta.data$newMonocClust == "Intermideate"]
odc = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "ODC" ]

opcOdcSub@meta.data = editClusterNames(opcOdcSub@meta.data)

# Sema3e
curClust = subset(opcOdcSub, newMonocClust == "OPC")
Idents(curClust) = curClust$group
curExpr = FindMarkers(curClust, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "MAST", assay = "PredictActivity")
#curExpr1 = FindMarkers(curClust, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
curExpr$Gene = rownames(curExpr)
Sema3e = curExpr[curExpr$Gene == "Sema3e",]

Gpr17 = curExpr[curExpr$Gene == "Gpr17",]

# Cadm1
curClust = subset(opcOdcSub, Cluster == "OPC_Inter")
Idents(curClust) = curClust$group
curExpr = FindMarkers(curClust, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "MAST", assay = "PredictActivity")
curExpr$Gene = rownames(curExpr)
Cadm1 = curExpr[curExpr$Gene == "Cadm1",]

# plots
atacFiltSub$newMonocClust = opcOdcSub$newMonocClust
atacFiltSub$Cluster = opcOdcSub$Cluster

atacFiltSub<- LinkPeaks(
  object = atacFiltSub ,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = c('Sema3e','Cadm1', 'Gpr17')
)

clusters<-c('OPC', "Intermideate", 'ODC')
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

p1 = CoveragePlot(
  object = atacFiltSub ,
  region = "Sema3e",
  features = "Sema3e",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 20000,
  extend.downstream = 20000, split.by = "Annotations") & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p1)

p2 = CoveragePlot(
  object = atacFiltSub ,
  region = "Cadm1",
  features = "Cadm1",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 20000,
  extend.downstream = 20000, split.by = "Annotations") & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p2)

p3 = CoveragePlot(
  object = atacFiltSub ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 20000,
  extend.downstream = 20000, split.by = "newMonocClust") & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p3)

targDir = './Paper_figs/Fig3/'

ggsave(paste0(targDir, "CovPlot_Sema3e_2024-02-06.png"), plot = p1, height = 12, width = 16, units = 'in', dpi = 300)

ggsave(paste0(targDir, "CovPlot_Cadm1_2024-02-06.png"), plot = p2, height = 12, width = 16, units = 'in', dpi = 300)

ggsave(paste0(targDir, "CovPlot_Gpr17_2024-02-06.png"), plot = p3, height = 12, width = 16, units = 'in', dpi = 300)

##
opc = subset(atacFiltSub, newMonocClust == "OPC")

opc<- LinkPeaks(
  object = opc ,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = c('Sema3e','Cadm1', 'Gpr17')
)

p4 = CoveragePlot(
  object = opc ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 20000,
  extend.downstream = 20000, split.by = "newMonocClust") & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p4)

# pseudo bulk
p3 = CoveragePlot(
  object = atacFiltSub ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 20000,
  extend.downstream = 20000, split.by = "newMonocClust", show.bulk = T) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 


p4 = CoveragePlot(
  object = opc ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  group.by = "dataset",
  extend.upstream = 0,
  extend.downstream = 5000, split.by = "newMonocClust") & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 