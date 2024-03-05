library(Seurat)
library(Signac)
#library(hdf5r)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

counts <- Read10X_h5("/home/flyhunter/Wang/data/stress/outs/extracted/filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "/home/flyhunter/Wang/data/stress/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz',
  min.cells = 1
)

brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

rm(brain_assay, counts, metadata, per_barcode_metrics)

gc()


# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(brain) <- annotations

# quality control based on RNAseq
controlCellNames<-read.table('/home/flyhunter/Wang/data/cellNames_doubletFinderFiltered_stress.txt', F)
controlCells<-controlCellNames$V1
brainFiltRna <- subset(brain, subset = gex_barcode %in% controlCells )

# quality control with ATACseq
brainFiltRna <- NucleosomeSignal(object = brainFiltRna)

metaCols<-colnames(brainFiltRna@meta.data)

brainFiltRna$nucleosome_group <- ifelse(brainFiltRna$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = brainFiltRna, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

brainFiltRna <- TSSEnrichment(brainFiltRna, fast = FALSE)

brainFiltRna$high.tss <- ifelse(brainFiltRna$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(brainFiltRna, group.by = 'high.tss') + NoLegend()

brainFiltRna$pct_reads_in_peaks <- brainFiltRna$atac_peak_region_fragments / brainFiltRna$atac_fragments * 100

VlnPlot(
  object = brainFiltRna,
  features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"),
  ncol = 4,
  pt.size = 0
)

print(quantile(brainFiltRna$pct_reads_in_peaks, probs = 0.01))

print(quantile(brainFiltRna$pct_reads_in_peaks, probs = 0.98))

low_atac<-quantile(brainFiltRna$nCount_peaks, probs = 0.02)

low_tss<-quantile(brainFiltRna$TSS.enrichment, probs = 0.01)

quantile(brainFiltRna$nCount_peaks, probs = 0.98)


brainAllFilt <- subset(
  x = brainFiltRna,
  subset = nCount_peaks < 150000 &
    nCount_peaks > low_atac &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    pct_reads_in_peaks > 5 &
    pct_reads_in_peaks < 52
)


peaks <- CallPeaks(
  object = brainAllFilt,
  effective.genome.size = 2.6e+09,
  macs2.path= '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
  
)


peaks1 <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
peaks1 <- keepStandardChromosomes(peaks1, pruning.mode = "coarse", species="Mus_musculus")

saveRDS(peaks1, file = 'macs2_peaks_stress')

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(brainAllFilt),
  features = peaks1,
  cells = colnames(brainAllFilt)
)

saveRDS(macs2_counts, file = 'macs2_counts_stress')

fragpath <- "/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz"

brainAllFilt[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotations
)

DefaultAssay(brainAllFilt)<-'macs2'

saveRDS(brainAllFilt, file = 'atacAllFiltMacs2_stress')