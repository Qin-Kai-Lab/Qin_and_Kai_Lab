library(ggplot2)
library(Seurat)

combObj = readRDS("SctInteg_selFov.rds")

pcs = c(10, 15, 20, 25 , 30 , 35, 40, 45, 50)
resols<-c(0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

for (curPc in pcs ) {
  # clustering
  combObj <- FindNeighbors(combObj, dims = 1:curPc, reduction = "integrated.dr")
  dir.create(paste0("Integrated_Clustering/Images/", curPc), recursive = T, showWarnings = F)
  for (curResol in resols) {
    combObj <- FindClusters(combObj, resolution = curResol)
    combObj <- RunUMAP(combObj, dims = 1:curPc, reduction = "integrated.dr")
    p1 = ImageDimPlot(combObj, fov = "b3_fov", axes = TRUE, group.by = "ident") +
      ggtitle("B3")
    
    p2 = ImageDimPlot(combObj, fov = "a5_fov", axes = TRUE, group.by = "ident")+
      ggtitle("A5")
    
    p3 = p2 + p1
    
    ggsave(paste0("Integrated_Clustering/Images/", curPc, "/Resol_", curResol, ".png"), plot=p3, height = 14, width = 20, units = 'in', dpi = 300)
    
  }
}


ImageDimPlot(combObj, fov = "a5_fov", axes = TRUE, group.by = "Slide_Group")

ImageDimPlot(combObj, fov = "b3_fov", axes = TRUE, group.by = "Slide_Group")