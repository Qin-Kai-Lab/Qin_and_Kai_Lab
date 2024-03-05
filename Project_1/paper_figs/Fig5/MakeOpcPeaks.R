library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

setwd("~/Wang/output")

opcOdc =  readRDS('integ_OPC_ODC')
opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")
opcOdc$Cell_ID = rownames(opcOdc@meta.data)
atacFiltSub = subset(atacFilt, Merged_CellName%in%rownames(opcOdc@meta.data))
opcOdcSub = subset(opcOdc, Cell_ID%in%rownames(atacFiltSub@meta.data))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

opcContr = opcOdcSub@meta.data$CellName[opcOdcSub@meta.data$newMonocClust == "OPC" & opcOdcSub@meta.data$group == "Control"]
opcStress = opcOdcSub@meta.data$CellName[opcOdcSub@meta.data$newMonocClust == "OPC" & opcOdcSub@meta.data$group == "Stress"]

rm(atacFilt, atacFiltSub, opcOdc)
gc()

# control
counts <- Read10X_h5("/home/flyhunter/Wang/data/control/outs/extracted/filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "/home/flyhunter/Wang/data/control/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

control_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz',
  min.cells = 1
)

control_obj <- CreateSeuratObject(
  counts = control_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)


rm(control_assay, counts, metadata)

control_obj$Cell_ID = rownames(control_obj@meta.data)

control_obj = subset(control_obj, Cell_ID%in%opcContr)

contrPeaks <- CallPeaks(
  object = control_obj,
  effective.genome.size = 1.87e+09,
  macs2.path= '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
)

contrPeaks <- subsetByOverlaps(x = contrPeaks, ranges = blacklist_mm10, invert = TRUE)
contrPeaks <- keepStandardChromosomes(contrPeaks, pruning.mode = "coarse", species="Mus_musculus")



# stress
counts <- Read10X_h5("/home/flyhunter/Wang/data/stress/outs/extracted/filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "/home/flyhunter/Wang/data/stress/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

stress_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz',
  min.cells = 1
)

stress_obj <- CreateSeuratObject(
  counts = stress_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)


rm(stress_assay, counts, metadata)

stress_obj$Cell_ID = rownames(stress_obj@meta.data)

stress_obj = subset(stress_obj, Cell_ID%in%opcStress)

stressPeaks <- CallPeaks(
  object = stress_obj,
  effective.genome.size = 1.87e+09,
  macs2.path= '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
)

stressPeaks <- subsetByOverlaps(x = stressPeaks, ranges = blacklist_mm10, invert = TRUE)
stressPeaks <- keepStandardChromosomes(stressPeaks, pruning.mode = "coarse", species="Mus_musculus")

identical(stressPeaks, contrPeaks)

# merge

combined.peaks <- reduce(x = c(contrPeaks, stressPeaks))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

control_frag<-CreateFragmentObject('/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz', cells = opcContr)

stress_frag<-CreateFragmentObject('/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz', cells = opcStress)

#identical(control_frag, stress_frag)

control_count<- FeatureMatrix(fragments = control_frag, features = combined.peaks)
#
stress_count<- FeatureMatrix(fragments = stress_frag, features = combined.peaks)

identical(control_count, stress_count)

###
###
control_chrom<-CreateChromatinAssay(control_count, fragments = control_frag)

stress_chrom<-CreateChromatinAssay(stress_count, fragments = stress_frag)

contrObj <- CreateSeuratObject(
  counts = control_chrom,
  assay = 'combined_peaks',
  project = 'ATAC'
)

stressObj <- CreateSeuratObject(
  counts = stress_chrom,
  assay = 'combined_peaks',
  project = 'ATAC'
)

contrObj$Group = "Control"
stressObj$Group = "Stress"
contrObj$Orig_ID = rownames(contrObj@meta.data)
stressObj$Orig_ID = rownames(stressObj@meta.data)
contrObj$Merge_CellID =paste0(contrObj$Orig_ID, "_", 1)
stressObj$Merge_CellID =paste0(stressObj$Orig_ID, "_", 2)

combObj = merge(contrObj, stressObj)

combObj = NormalizeData(combObj)
combObj = ScaleData(combObj)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(combObj) <- annotations

opcSubRna= subset(opcOdcSub, newMonocClust == "OPC")

identical(opcSubRna@meta.data$Cell_ID, combObj@meta.data$Merge_CellID)

combObj = RenameCells(combObj, new.names = combObj@meta.data$Merge_CellID)

combObj[["RNA"]] = opcSubRna[["RNA"]]


saveRDS(combObj, "OPC_CombPeaks_macs2_2.rds")