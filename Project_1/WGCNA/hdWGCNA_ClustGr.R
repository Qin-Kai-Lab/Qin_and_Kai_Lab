library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(irlba)
library(corrplot)
library(igraph)

theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)
# optionally enable multithreading
enableWGCNAThreads(nThreads = 14)

curDate<-Sys.Date()
source('../programs/renameClusters.R')
setwd("/home/flyhunter/Wang/output")

targDir = "hdWGCNA/StrVsContr/"
dir.create(targDir, recursive = T, showWarnings = F)

##
prepDat = function(obj, cluster, cond) {
  seurat_obj = subset(obj, subset = ((Annotations == cluster) & (group == cond)))
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction", 
    fraction = 0.05, 
    wgcna_name = paste(cluster, cond, sep = "_") 
  )
  
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("Annotations"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    k = 25, # nearest-neighbors parameter
    max_shared = 15, # maximum number of shared cells between two metacells
    ident.group = 'Annotations', # set the Idents of the metacell seurat object
    min_cells = 89
  )
  
  seurat_obj <- NormalizeMetacells(seurat_obj)
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  #seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj), fastpath=FALSE)
  #seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='pca', dims=1:15)
  
  return(seurat_obj)
}


getModules = function(obj, cluster, cond, targDir) {
  
  seurat_obj <- SetDatExpr(
    obj,
    assay = 'RNA', # using RNA assay
    slot = 'data' # using normalized data
  )
  
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )
  
  seurat_obj <- ConstructNetwork(seurat_obj, setDatExpr=FALSE,
    tom_name = paste(cluster, cond, sep = "_"),  overwrite_tom = TRUE, tom_outdir = paste0(targDir,"TOM"))
  
  seurat_obj <- ModuleEigengenes(seurat_obj)
  
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  write.csv(MEs, paste0(targDir, "MEs_", cluster, "_", cond, ".csv"), row.names = T)
  
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
  )
  # rename the modules
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste(cluster, cond, sep = "_") 
  )
  
  p <- PlotKMEs(seurat_obj, ncol=5)
  png(paste0(targDir, "ModulesPlots_", cluster, "_", cond, ".png"), height = 14, width = 16, units = "in", res = 300)
  print(p)
  dev.off()
  
  modules <- GetModules(seurat_obj)
  hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
  write.csv(modules,  paste0(targDir, "Modules_", cluster, "_", cond, ".csv"), row.names = T)
  write.csv(hub_df,  paste0(targDir, "Top10HubGenes_", cluster, "_", cond, ".csv"), row.names = T)
  
  ModuleNetworkPlot(seurat_obj, outdir = paste0(targDir,"ModuleNetworks_", cluster, "_", cond), plot_size = c(14, 14))
  
}


# runn All
clusters = c('CA2', 'GABA', 'MG', "ODC")

for (cluster in clusters) {
  try({
    cond = "Stress"
    seurat_obj = prepDat(obj = RNA.combined.norm, cluster = cluster, cond = cond)
    getModules(obj = seurat_obj, cluster = cluster, cond = cond, targDir = targDir)
  })
}