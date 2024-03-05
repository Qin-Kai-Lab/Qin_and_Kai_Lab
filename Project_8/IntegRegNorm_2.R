# current plan is to merge the A5 and B3 before integration split additionally in groups for integration and merform integration


library(Seurat)
library(ggplot2)
library(clustree)
library(MAST)
#library(future)
#plan("multisession", workers = 10)
#options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8/output")



addGroups = function(objSubFilt) {
  a = 1:98
  b = 99:196
  objSubFilt@meta.data$Slide_Group = NA
  for (i in 1:nrow(objSubFilt@meta.data) ) {
    if (objSubFilt@meta.data$fov[i] <= 98) {
      objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$Slide[i], "A", sep = "_")
    } else if (objSubFilt@meta.data$fov[i] >= 99) {
      objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$Slide[i], "B", sep = "_")
    }
  }
  return(objSubFilt)
}

plotUmap<-function(x, curPc){
  for (i in x){
    combObj <- FindClusters(combObj, resolution = i)
    combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
    current_resol<-paste0("Integrated_RegNorm_snn_res.", i)
    #Idents(combObj)<-current_resol
    p2 <- DimPlot(combObj, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
    ggsave(paste0("Integrated_RegNorm_Clustering/Resols/", curPc, "/Resol_",i, ".png"), plot=p2, height = 12, width = 14, units = 'in', dpi = 300)
    
  }
}

# a5Obj = readRDS("A5_selFovs.rds")
# a5Obj$Group = "A5"
# a5Obj$Slide = 1
# 
# a5Fov_list = strsplit(rownames(a5Obj@meta.data), "_")
# a5Fov <- as.numeric(as.character(sapply(a5Fov_list, function(vec) vec[2])))
# a5Obj$fov = a5Fov
# 
# b3Obj = readRDS("B3_selFovs.rds")
# b3Obj$Group = "B3"
# b3Obj$Slide = 2
# b3Fov_list = strsplit(rownames(b3Obj@meta.data), "_")
# b3Fov <- as.numeric(as.character(sapply(b3Fov_list, function(vec) vec[2])))
# b3Obj$fov = b3Fov
# 
# a5Obj = addGroups(a5Obj)
# b3Obj = addGroups(b3Obj)
# 
# combObj = merge(a5Obj, b3Obj )
# combObj[["Nanostring"]] <- JoinLayers(combObj[["Nanostring"]])
# combObj[["Nanostring"]] <- split(combObj[["Nanostring"]], f = combObj$Slide_Group)
# 
# ImageDimPlot(combObj, fov = "a5_fov", axes = TRUE, cols = "glasbey")
# ImageDimPlot(combObj, fov = "b3_fov", axes = TRUE, cols = "glasbey")
# 
# combObj <- NormalizeData(combObj)
# combObj <- FindVariableFeatures(combObj)
# combObj <- ScaleData(combObj)
# combObj <- RunPCA(combObj, npcs = 50)
# 
# DefaultAssay(combObj)
# 
# combObj <- IntegrateLayers(object = combObj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.dr",verbose = FALSE)
# elbPlot = ElbowPlot(combObj, ndims = 50)
# ggsave("RegNorm_2_Elbow.png", plot = elbPlot, height = 10, width = 10, dpi = 300, units = "in")
# 
# 
dir.create("Integrated_RegNorm_Clustering/PCA/", recursive = T)
dir.create("Integrated_RegNorm_Clustering/Resols/", recursive = T)

#saveRDS(combObj, "RegNormInteg_selFov.rds")

combObj = readRDS("RegNormInteg_selFov.rds")

pcs = c(10, 15, 20, 25 , 30 , 35, 40, 45, 50)
resols<-c(0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (curPc in pcs ) {
  # clustering
  combObj <- FindNeighbors(combObj, dims = 1:curPc, reduction = "integrated.dr")
  combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
  
  curPcPlot = DimPlot(combObj, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_RegNorm_Clustering/PCA/PCs_",curPc, ".png"), plot = curPcPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  contrGr = c("1_A", "2_B")
  contr = subset(combObj, Slide_Group %in% contrGr)
  contrPlot = DimPlot(contr, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_RegNorm_Clustering/PCA/PCs_",curPc, "_Controls.png"), plot=contrPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  dir.create(paste0("Integrated_RegNorm_Clustering/Resols/", curPc), recursive = T, showWarnings = F)
  clustree = FindClusters(combObj, resolution = resols)
  p <- clustree(clustree)
  ggsave(paste0("Integrated_RegNorm_Clustering/PCA/ClustTree_PCs_",curPc, ".png"), plot = p, height = 12, width = 14, units = 'in', dpi = 300)
  
  plotUmap(resols, curPc)
  
}