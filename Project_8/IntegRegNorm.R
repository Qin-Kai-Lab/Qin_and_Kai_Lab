library(Seurat)
library(ggplot2)
library(clustree)
library(MAST)

# dat = readRDS("/home/flyhunter/Projects/Project8/Seurat_Object_combined_slides.RDS")
# q1 = unique(data.frame(dat@meta.data$fov, dat@meta.data$Run_Tissue_name))
# table(q1$dat.meta.data.Run_Tissue_name)
# DimPlot(dat)

setwd("/home/flyhunter/Projects/Project8/output")

dat2 = readRDS("/home/flyhunter/Projects/Project8/Seurat_object_combined_slides_quality_control.RDS")
# DimPlot(dat2)
# ImageDimPlot(dat, fov = '1', axes = TRUE, cols = "glasbey")

#dat2 = SCTransform(dat2)

addGroups = function(objSubFilt) {
  a = 1:98
  b = 99:196
  objSubFilt@meta.data$Slide_Group = NA
  for (i in 1:nrow(objSubFilt@meta.data) ) {
    if (objSubFilt@meta.data$fov[i] <= 98) {
      objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$slide_ID_numeric[i], "A", sep = "_")
    } else if (objSubFilt@meta.data$fov[i] >= 99) {
      objSubFilt@meta.data$Slide_Group[i] = paste(objSubFilt@meta.data$slide_ID_numeric[i], "B", sep = "_")
    }
  }
  return(objSubFilt)
}

plotUmap<-function(x, curPc){
  for (i in x){
    objSubFilt <- FindClusters(objSubFilt, resolution = i)
    objSubFilt <- RunUMAP(objSubFilt, dims = 1:curPc, reduction = "integrated.dr")
    current_resol<-paste0("Integrated_RegNorm_snn_res.", i)
    #Idents(objSubFilt)<-current_resol
    p2 <- DimPlot(objSubFilt, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
    ggsave(paste0("Integrated_RegNorm_Clustering/Resols/", curPc, "/Resol_",i, ".png"), plot=p2, height = 12, width = 14, units = 'in', dpi = 300)
    
  }
}

a5_fovs = c(9, 10, 16, 17, 23, 24, 25, 30, 31, 32, 33, 34, 38, 39, 40, 41, 
            69, 75, 76, 79, 80, 81, 82, 83, 86, 87, 88, 89,
            156, 157, 163, 164, 170, 171, 177, 178, 179, 180, 184, 185, 
            186, 187, 188, 192, 193, 194, 195,
            108, 109, 110, 115, 116, 117, 118, 123, 124, 125, 129, 130, 131, 132, 135, 136, 137, 138, 144)

b3_fovs = c(9, 10, 15, 16, 17, 22, 23, 24, 25, 26, 30, 31, 32, 33, 34, 38, 39, 40,
            54, 55, 61, 62, 63, 68, 69, 70, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 88,
            113, 114, 120, 121, 122, 123, 128, 129, 130, 131, 132, 136, 137, 138, 139,
            167, 168, 173, 174, 175, 177, 178, 179, 180, 181, 182, 184, 185, 186, 187)


q1 = unique(dat2@meta.data[, c('slide_ID_numeric', 'Run_Tissue_name')])

a5_fovs = paste0("1_", a5_fovs)
b3_fovs = paste0("2_", b3_fovs)

selFovs = c(a5_fovs, b3_fovs)

dat2$Slide_fov = paste(dat2$slide_ID_numeric, dat2$fov, sep = "_")

objSub = subset(dat2, Slide_fov%in%selFovs)

colnames(objSub@meta.data)

q1 = unique(objSub@meta.data[, c("qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg", "qcFlagsCellComplex", "qcFlagsCellArea")])

objSubFilt = subset(objSub, qcFlagsRNACounts == "Pass" & qcFlagsCellCounts == "Pass" & qcFlagsCellPropNeg == "Pass" & qcFlagsCellComplex == "Pass" & qcFlagsCellArea == "Pass")

objSubFilt = addGroups(objSubFilt)

objSubFilt[["RNA"]] <- split(objSubFilt[["RNA"]], f = objSubFilt$Slide_Group)

objSubFilt <- NormalizeData(objSubFilt)
objSubFilt <- FindVariableFeatures(objSubFilt)
objSubFilt <- ScaleData(objSubFilt)
objSubFilt <- RunPCA(objSubFilt, npcs = 50)

DefaultAssay(objSubFilt)

objSubFilt <- IntegrateLayers(object = objSubFilt, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.dr",verbose = FALSE)
ElbowPlot(objSubFilt, ndims = 50)



dir.create("Integrated_RegNorm_Clustering/PCA/", recursive = T)
dir.create("Integrated_RegNorm_Clustering/Resols/", recursive = T)

saveRDS(objSubFilt, "RegNormInteg")

pcs = c(10, 15, 20, 25 , 30 , 35, 40, 45, 50)
resols<-c(0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (curPc in pcs ) {
  # clustering
  objSubFilt <- FindNeighbors(objSubFilt, dims = 1:curPc, reduction = "integrated.dr")
  objSubFilt <- RunUMAP(objSubFilt, dims = 1:curPc, reduction = "integrated.dr")
  
  curPcPlot = DimPlot(objSubFilt, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_RegNorm_Clustering/PCA/PCs_",curPc, ".png"), plot = curPcPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  contrGr = c("1_A", "2_B")
  contr = subset(objSubFilt, Slide_Group %in% contrGr)
  contrPlot = DimPlot(contr, reduction = "umap", group.by = "Slide_Group", label = F)
  ggsave(paste0("Integrated_RegNorm_Clustering/PCA/PCs_",curPc, "_Controls.png"), plot=contrPlot, height = 12, width = 14, units = 'in', dpi = 300)
  
  dir.create(paste0("Integrated_RegNorm_Clustering/Resols/", curPc), recursive = T, showWarnings = F)
  
  plotUmap(resols, curPc)
  
}