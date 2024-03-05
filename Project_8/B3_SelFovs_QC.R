library(Seurat)
#library(future)
#plan("multisession", workers = 10)
#options(future.globals.maxSize = 8000 * 1024^2)
setwd("~/Projects/Project_8")

data.dir = "flatfile_brain_new/SlideB3"
fov = "b3_fov"
assay = "Nanostring"
metadata = read.csv("/home/flyhunter13/Projects/Project_8/flatfile_brain_new/SlideB3/SlideB3_metadata_file.csv.gz")

sel_fovs = c(9, 10, 15, 16, 17, 22, 23, 24, 25, 26, 30, 31, 32, 33, 34, 38, 39, 40, 2, 3, 8, 18, 27, 29, 37,
             54, 55, 61, 62, 63, 68, 69, 70, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 88, 67, 84, 86, 87, 89, 90,
             113, 114, 120, 121, 122, 123, 128, 129, 130, 131, 132, 136, 137, 138, 139, 115, 124, 127, 135,
             167, 168, 173, 174, 175, 177, 178, 179, 180, 181, 182, 184, 185, 186, 187, 172, 188,
             72, 106, 107, 160, 161)

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

ImageDimPlot(objFiltSub, fov = "b3_fov", axes = TRUE)

objFiltSub$NewCellId = as.character(paste(rownames(objFiltSub@meta.data), "B3", sep = "_"))
objFiltSub = RenameCells(objFiltSub, new.names = objFiltSub@meta.data$NewCellId)

setwd("output")

saveRDS(objFiltSub, "B3_selFovs_QC.rds")