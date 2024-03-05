library(Seurat)
library(ggplot2)
library(clustree)
library(MAST)
#library(future)
#plan("multisession", workers = 10)
#options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8/output")



# addGroups = function(objSubFilt) {
#   a = 1:98
#   b = 99:196
#   objSubFilt@meta.data$Slide_Group = NA
#   for (i in 1:nrow(objSubFilt@meta.data) ) {
#     if (objSubFilt@meta.data$fov[i] <= 98) {
#       objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$Slide[i], "A", sep = "_")
#     } else if (objSubFilt@meta.data$fov[i] >= 99) {
#       objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$Slide[i], "B", sep = "_")
#     }
#   }
#   return(objSubFilt)
# }
# 
# # plotUmap<-function(x, curPc){
# #   for (i in x){
# #     combObj <- FindClusters(combObj, resolution = i)
# #     combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
# #     current_resol<-paste0("Integrated_RegNorm_snn_res.", i)
# #     #Idents(combObj)<-current_resol
# #     p2 <- DimPlot(combObj, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
# #     ggsave(paste0("Integrated_Clustering_SCT/Resols/", curPc, "/Resol_",i, ".png"), plot=p2, height = 12, width = 14, units = 'in', dpi = 300)
# #     
# #   }
# # }
# 
# a5Obj = readRDS("A5_selFovs_QC.rds")
# a5Obj$Group = "A5"
# a5Obj$Slide = 1
# 
# # a5Fov_list = strsplit(rownames(a5Obj@meta.data), "_")
# # a5Fov <- as.numeric(as.character(sapply(a5Fov_list, function(vec) vec[2])))
# # a5Obj$fov = a5Fov
# 
# b3Obj = readRDS("B3_selFovs_QC.rds")
# b3Obj$Group = "B3"
# b3Obj$Slide = 2
# # b3Fov_list = strsplit(rownames(b3Obj@meta.data), "_")
# # b3Fov <- as.numeric(as.character(sapply(b3Fov_list, function(vec) vec[2])))
# # b3Obj$fov = b3Fov
# 
# a5Obj = addGroups(a5Obj)
# b3Obj = addGroups(b3Obj)
# 
# combObj = merge(a5Obj, b3Obj )
# combObj[["Nanostring"]] <- JoinLayers(combObj[["Nanostring"]])
# combObj[["Nanostring"]] <- split(combObj[["Nanostring"]], f = combObj$Slide_Group)
# 
# combObj <- SCTransform(combObj, assay = "Nanostring")
# combObj <- RunPCA(combObj)
# 
# DefaultAssay(combObj)
# 
# combObj <- IntegrateLayers(object = combObj, method = CCAIntegration, normalization.method = "SCT", verbose = F)
# elbPlot = ElbowPlot(combObj, ndims = 50)
# png("Sct_3_Elbow.png", height = 10, width = 10, res = 300, units = "in")
# print(elbPlot)
# dev.off()
# # 
# # 
# dir.create("Integrated_Clustering_SCT/PCA/", recursive = T)
# dir.create("Integrated_Clustering_SCT/Resols/", recursive = T)
# # 
# saveRDS(combObj, "SctInteg_3_selFov.rds")

combObj = readRDS("SctInteg_3_selFov.rds")

pcs = c(10, 15, 20, 25 , 30 , 35, 40, 45, 50)
resols<-c(0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

for (curPc in pcs ) {
  # clustering
  combObj <- FindNeighbors(combObj, dims = 1:curPc, reduction = "integrated.dr")
  combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
  
  curPcPlot = DimPlot(combObj, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_Clustering_SCT/PCA/PCs_",curPc, ".png"), plot = curPcPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  contrGr = c("1_A", "2_B")
  contr = subset(combObj, Slide_Group %in% contrGr)
  contrPlot = DimPlot(contr, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_Clustering_SCT/PCA/PCs_",curPc, "_Controls.png"), plot=contrPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  combObj= FindClusters(combObj, resolution = resols)
  p <- clustree(combObj)
  ggsave(paste0("Integrated_Clustering_SCT/PCA/ClustTree_PCs_",curPc, ".png"), plot = p, height = 12, width = 14, units = 'in', dpi = 300)
  
  dir.create(paste0("Integrated_Clustering_SCT/Images/", curPc), recursive = T, showWarnings = F)
  dir.create(paste0("Integrated_Clustering_SCT/Resols/", curPc), recursive = T, showWarnings = F)
  #combObj <- FindClusters(combObj, resolution = resols)
  for (curResol in resols) {
    #combObj <- FindClusters(combObj, resolution = curResol)
    #combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
    current_resol<-paste0("SCT_snn_res.", curResol)
    Idents(combObj) = combObj@meta.data[, current_resol]
    
    p1 = ImageDimPlot(combObj, fov = "b3_fov", axes = TRUE, group.by = "ident") +
      ggtitle("B3")
    
    p2 = ImageDimPlot(combObj, fov = "a5_fov", axes = TRUE, group.by = "ident")+
      ggtitle("A5")
    
    p3 = p2 + p1
    ggsave(paste0("Integrated_Clustering_SCT/Images/", curPc, "/Resol_", curResol, ".png"), plot=p3, height = 14, width = 20, units = 'in', dpi = 300)
    
    p4 <- DimPlot(combObj, reduction = "umap", label = T, repel=T, raster=FALSE)+ggtitle(current_resol)
    ggsave(paste0("Integrated_Clustering_SCT/Resols/", curPc, "/Resol_",curResol, ".png"), plot=p4, height = 12, width = 14, units = 'in', dpi = 300)
    
  }
  
}