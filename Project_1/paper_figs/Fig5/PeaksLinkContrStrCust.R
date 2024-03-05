library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)
# library(future)
# plan("multisession", workers = 4)
# options(future.globals.maxSize = 50 * 1024 ^ 3)

setwd("~/Wang/output/")

targDir = './Paper_figs/Fig5/'


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

opc_inter = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "OPC" | opcOdcSub@meta.data$newMonocClust == "Intermideate"]
odc = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "ODC" ]

opcOdcSub@meta.data = editClusterNames(opcOdcSub@meta.data)

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

atacFiltSub$newMonocClust = opcOdcSub$newMonocClust
atacFiltSub$Cluster = opcOdcSub$Cluster

opcObj = subset(atacFiltSub, newMonocClust == "OPC")

rm(atacFilt, atacFiltSub, opcOdc, opcOdcSub)

gc()

Groups = unique(opcObj$dataset)

for (curGroup in Groups) {
  subObj = subset(opcObj, dataset == curGroup)
  try({
    peakRna = LinkPeaks(
      subObj,
      peak.assay = "Combined_peaks",
      expression.assay = "RNA",
      peak.slot = "counts",
      expression.slot = "data",
      method = "pearson",
      distance = 5000,
      min.cells = 0,
      genes.use = "Gpr17",
      n_sample = 200,
      pvalue_cutoff = 1,
      score_cutoff = 0,
      gene.id = FALSE,
      verbose = TRUE)
    
    
    p2g = data.frame(Links(peakRna))
    p2g$Group = curGroup
    write.csv(p2g, file = paste0(targDir, "Peaks_Gpr17_Cor_5KB_", curGroup, "_2024-02-12.csv"), row.names = F)
      
  })
  
}

control =  read.csv("~/Wang/output/Paper_figs/Fig5/Peaks_Gpr17_Cor_500KB_Control_2024-02-12.csv")
stress = read.csv("~/Wang/output/Paper_figs/Fig5/Peaks_Gpr17_Cor_500KB_Stress_2024-02-12.csv")

control = control[control$pvalue < 0.05,]
stress  = stress[stress$pvalue < 0.05,]

strUniq = stress[!stress$peak %in% control$peak,]

