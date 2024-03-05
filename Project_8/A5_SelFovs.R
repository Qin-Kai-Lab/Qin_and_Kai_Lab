library(Seurat)
library(future)
plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8")

data.dir = "flatfile_brain_new/SlideA5"
fov = "a5_fov"
assay = "Nanostring"

sel_fovs = c(9, 10, 16, 17, 23, 24, 25, 30, 31, 32, 33, 34, 38, 39, 40, 41, 
            69, 75, 76, 79, 80, 81, 82, 83, 86, 87, 88, 89,
            156, 157, 163, 164, 170, 171, 177, 178, 179, 180, 184, 185, 
            186, 187, 188, 192, 193, 194, 195,
            108, 109, 110, 115, 116, 117, 118, 123, 124, 125, 129, 130, 131, 132, 135, 136, 137, 138, 144)

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

ImageDimPlot(obj, fov = "a5_fov", axes = TRUE, cols = "glasbey")

obj$NewCellId = as.character(paste(rownames(obj@meta.data), "A5", sep = "_"))
obj = RenameCells(obj, new.names = obj@meta.data$NewCellId)

setwd("output")

saveRDS(obj, "A5_selFovs.rds")