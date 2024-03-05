library(Seurat)
library(Signac)
#library(hdf5r)
library(EnsDb.Mmusculus.v79)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

# import QCed cell identities
controlCellNames<-read.table('cellNames_doubletFinderFiltered_control.txt', F)
controlCells<-controlCellNames$V1

# import seurat Control data
control = Read10X(data.dir = "../data/control")
RNA_control = CreateSeuratObject(counts = control$`Gene Expression`)
RNA_control$group = "Control"
RNA_control[["CellName"]] <- colnames(RNA_control)

# atac
#hcounts <- Read10X_h5('./control/extracted/filtered_feature_bc_matrix.h5')
#hcounts peaks and conrol peaks are identical


RNA_control[['ATAC']] <- CreateChromatinAssay(counts = control$Peaks,  sep=c(":", "-"), genome = "mm10",
                                              fragments = './control/atac_fragments.tsv.gz', min.cells=1)

DefaultAssay(RNA_control) <- 'ATAC'

granges(RNA_control)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(RNA_control) <- annotations

# quality control based on RNAseq
controlCellNames<-read.table('cellNames_doubletFinderFiltered_control.txt', F)
controlCells<-controlCellNames$V1
controlFiltRna <- subset(RNA_control, subset = CellName %in% controlCells )

rm(RNA_control, control)
rm(annotations, controlCellNames, controlCells)
gc()

# quality control with ATACseq
controlFiltRna <- NucleosomeSignal(object = controlFiltRna)

controlFiltRna$nucleosome_group <- ifelse(controlFiltRna$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = controlFiltRna, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

controlFiltRna <- TSSEnrichment(controlFiltRna, fast = FALSE)

controlFiltRna$high.tss <- ifelse(controlFiltRna$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(controlFiltRna, group.by = 'high.tss') + NoLegend()

#controlFiltRna$pct_reads_in_peaks <- controlFiltRna$peak_region_fragments / controlFiltRna$passed_filters * 100
#controlFiltRna$blacklist_ratio <- controlFiltRna$blacklist_region_fragments / controlFiltRna$peak_region_fragments

VlnPlot(
  object = controlFiltRna,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


low_atac<-quantile(controlFiltRna$nCount_ATAC, probs = 0.02)

low_tss<-quantile(controlFiltRna$TSS.enrichment, probs = 0.01)

controlAllFilt <- subset(
  x = controlFiltRna,
  subset = nCount_ATAC < 110000 &
    nCount_ATAC > low_atac &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

saveRDS(controlAllFilt, file = 'atacRnaControl')
