controlCells = read.table('OPC_ODC_contrCells.txt')
controlCells = controlCells$V1


library(Seurat)
library(Signac)
#library(hdf5r)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)


setwd('/home/flyhunter/Wang/data')

counts <- Read10X_h5("/home/flyhunter/Wang/data/control/outs/extracted/filtered_feature_bc_matrix.h5")
metadata <- read.csv(
  file = "/home/flyhunter/Wang/data/control/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz',
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


brainFiltRna <- subset(brain, subset = gex_barcode %in% controlCells )

# call peaks with macs2

peaks <- CallPeaks(
  object = brainFiltRna,
  effective.genome.size = 1.87e+09,
  macs2.path= '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
  
)

peaks1 <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
peaks1 <- keepStandardChromosomes(peaks1, pruning.mode = "coarse", species="Mus_musculus")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(brainFiltRna),
  features = peaks1,
  cells = colnames(brainFiltRna)
)

#fragpath <- "/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz"

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(brain) <- annotations


brainFiltRna[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(brainFiltRna),
  annotation = annotations
)

DefaultAssay(brainFiltRna)<-'macs2'

brainFiltRna[["peaks"]]<-NULL

# filter
brainFiltRna <- NucleosomeSignal(object = brainFiltRna)

metaCols<-colnames(brainFiltRna@meta.data)

brainFiltRna$nucleosome_group <- ifelse(brainFiltRna$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = brainFiltRna, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

brainFiltRna <- TSSEnrichment(brainFiltRna, fast = FALSE)

brainFiltRna$high.tss <- ifelse(brainFiltRna$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(brainFiltRna, group.by = 'high.tss') + NoLegend()

#brainFiltRna$pct_reads_in_peaks <- brainFiltRna$atac_peak_region_fragments / brainFiltRna$atac_fragments * 100
metaCols<-colnames(brainFiltRna@meta.data)

VlnPlot(
  object = brainFiltRna,
  features = c("nCount_macs2", "TSS.enrichment", "nucleosome_signal", "nFeature_macs2"),
  ncol = 4,
  pt.size = 1
)

low_feature <- quantile(brainFiltRna$nFeature_macs2, probs = 0.01)

high_feature <- quantile(brainFiltRna$nFeature_macs2, probs = 0.99)

low_atac<-quantile(brainFiltRna$nCount_macs2, probs = 0.01)
high_atac<-quantile(brainFiltRna$nCount_macs2, probs = 0.99)

low_tss<-quantile(brainFiltRna$TSS.enrichment, probs = 0.01)

high_nss <- quantile(brainFiltRna$nucleosome_signal, probs = 0.99, na.rm = T)
low_nss <- quantile(brainFiltRna$nucleosome_signal, probs = 0.01, na.rm = T)

brainAllFilt <- subset(
  x = brainFiltRna,
  subset = nCount_macs2 < high_atac &
    nCount_peaks > low_atac &
    nucleosome_signal < high_nss &
    nucleosome_signal > low_nss &
    TSS.enrichment > low_tss & 
    nFeature_macs2 > low_feature &
    nFeature_macs2 < high_feature
)

brainAllFilt
brainFiltRna

saveRDS(brainAllFilt, file = 'OpcOdc_control_macs2')
