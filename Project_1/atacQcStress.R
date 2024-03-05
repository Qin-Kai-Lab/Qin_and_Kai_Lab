library(Seurat)
library(Signac)
#library(hdf5r)
library(EnsDb.Mmusculus.v79)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")


stressCellNames<-read.table('cellNames_doubletFinderFiltered_stress.txt', F)
stressCells<-stressCellNames$V1

stress = Read10X(data.dir = "../data/stress")
RNA_stress = CreateSeuratObject(counts = stress$`Gene Expression`)
RNA_stress$group = "Stress"
RNA_stress[["CellName"]] <- colnames(RNA_stress)


###
RNA_stress[['ATAC']] <- CreateChromatinAssay(counts = stress$Peaks,  sep=c(":", "-"), genome = "mm10",
                                              fragments = './stress/atac_fragments.tsv.gz', min.cells=1)

DefaultAssay(RNA_stress) <- 'ATAC'

granges(RNA_stress)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(RNA_stress) <- annotations

# quality control based on RNAseq

stressFiltRna <- subset(RNA_stress, subset = CellName %in% stressCells )

rm(RNA_stress, stress)
rm(annotations, stressCellNames, stressCells)
gc()

# quality control with ATACseq
stressFiltRna <- NucleosomeSignal(object = stressFiltRna)

stressFiltRna$nucleosome_group <- ifelse(stressFiltRna$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = stressFiltRna, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

stressFiltRna <- TSSEnrichment(stressFiltRna, fast = FALSE)

stressFiltRna$high.tss <- ifelse(stressFiltRna$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(stressFiltRna, group.by = 'high.tss') + NoLegend()

#controlFiltRna$pct_reads_in_peaks <- controlFiltRna$peak_region_fragments / controlFiltRna$passed_filters * 100
#controlFiltRna$blacklist_ratio <- controlFiltRna$blacklist_region_fragments / controlFiltRna$peak_region_fragments

VlnPlot(
  object = stressFiltRna,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


low_atac<-quantile(stressFiltRna$nCount_ATAC, probs = 0.02)
low_atac

low_tss<-quantile(stressFiltRna$TSS.enrichment, probs = 0.02)
low_tss


stressAllFilt <- subset(
  x = stressFiltRna,
  subset = nCount_ATAC < 150000 &
    nCount_ATAC > low_atac &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)


saveRDS(stressAllFilt, file = 'atacRnaStress')