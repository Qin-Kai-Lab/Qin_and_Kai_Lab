library(Seurat)
#library(future)
#plan("multisession", workers = 10)
#options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8")

data.dir = "flatfile_brain_new/SlideA5"
fov = "a5_fov"
assay = "Nanostring"
metadata = read.csv("/home/flyhunter13/Projects/Project_8/flatfile_brain_new/SlideA5/SlideA5_metadata_file.csv.gz")

sel_fovs = c(9, 10, 16, 17, 23, 24, 25, 30, 31, 32, 33, 34, 38, 39, 40, 41, 
             8, 15, 18, 22, 26, 37, 46, 47, 48,
             69, 75, 76, 79, 80, 81, 82, 83, 86, 87, 88, 89, 74, 77, 90,
             156, 157, 163, 164, 170, 171, 177, 178, 179, 180, 184, 185, 186, 187, 188, 192, 193, 194, 195,
             162, 169, 172, 176, 191,
             108, 109, 110, 115, 116, 117, 118, 123, 124, 125, 129, 130, 131, 132, 135, 136, 137, 138, 144, 101,
             102, 111, 122, 139, 142, 143, 145, 
             68, 70, 84, 149, 150, 155, 183, 100, 107)

data <- ReadNanostring(data.dir = data.dir, type = c("centroids", "segmentations"), fov.filter = sel_fovs)


segs <- CreateSegmentation(data$segmentations)
cents <- CreateCentroids(data$centroids)
segmentations.data <- list(centroids = cents, segmentation = segs)
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                          "centroids"), molecules = data$pixels, assay = assay)

# filter metadata by fovs
metadataFilt = metadata[metadata$fov%in%sel_fovs,]
metadataFilt$CurId = paste0(metadataFilt$cell_ID, '_', metadataFilt$fov)
rownames(metadataFilt) = metadataFilt$CurId

identical(rownames(metadataFilt), colnames(data$matrix))

obj <- CreateSeuratObject(counts = data$matrix, assay = assay, meta.data = metadataFilt)


VlnPlot(
  object = obj,
  features = c("nCount_Nanostring", "nFeature_Nanostring"),
  ncol = 2,
  pt.size = 1
)

low_feature <- quantile(obj$nFeature_Nanostring, probs = 0.01)

low_count <- quantile(obj$nCount_Nanostring, probs = 0.01)

objFilt = subset(obj, nFeature_Nanostring > low_feature & nCount_Nanostring > low_count)
objFiltSub = subset(objFilt, qcFlagsRNACounts == "Pass" & qcFlagsCellCounts == "Pass" & qcFlagsCellPropNeg == "Pass" & qcFlagsCellComplex == "Pass" & qcFlagsCellArea == "Pass")


cells <- intersect(Cells(x = coords, boundary = "segmentation"), 
                   Cells(x = coords, boundary = "centroids"))
cells <- intersect(Cells(objFiltSub), cells)
coords <- subset(x = coords, cells = cells)
objFiltSub[[fov]] <- coords

ImageDimPlot(objFiltSub, fov = "a5_fov", axes = TRUE)

objFiltSub$NewCellId = as.character(paste(rownames(objFiltSub@meta.data), "A5", sep = "_"))
objFiltSub = RenameCells(objFiltSub, new.names = objFiltSub@meta.data$NewCellId)

setwd("output")

saveRDS(objFiltSub, "A5_selFovs_QC.rds")