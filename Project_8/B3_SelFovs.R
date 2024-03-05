library(Seurat)
library(future)
plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8")

data.dir = "flatfile_brain_new/SlideB3"
fov = "b3_fov"
assay = "Nanostring"

sel_fovs = c(9, 10, 15, 16, 17, 22, 23, 24, 25, 26, 30, 31, 32, 33, 34, 38, 39, 40,
             54, 55, 61, 62, 63, 68, 69, 70, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 88,
             113, 114, 120, 121, 122, 123, 128, 129, 130, 131, 132, 136, 137, 138, 139,
             167, 168, 173, 174, 175, 177, 178, 179, 180, 181, 182, 184, 185, 186, 187)

data <- ReadNanostring(data.dir = data.dir, type = c("centroids", "segmentations"), fov.filter = sel_fovs)


segs <- CreateSegmentation(data$segmentations)
cents <- CreateCentroids(data$centroids)
segmentations.data <- list(centroids = cents, segmentation = segs)
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                          "centroids"), molecules = data$pixels, assay = assay)
obj <- CreateSeuratObject(counts = data$matrix, assay = assay)
cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                   Cells(x = coords, boundary = "centroids"))
cells <- intersect(Cells(obj), cells)
coords <- subset(x = coords, cells = cells)
obj[[fov]] <- coords

ImageDimPlot(obj, fov = "b3_fov", axes = TRUE, cols = "glasbey")

obj$NewCellId = as.character(paste(rownames(obj@meta.data), "B3", sep = "_"))
obj = RenameCells(obj, new.names = obj@meta.data$NewCellId)

setwd("output")

saveRDS(obj, "B3_selFovs.rds")