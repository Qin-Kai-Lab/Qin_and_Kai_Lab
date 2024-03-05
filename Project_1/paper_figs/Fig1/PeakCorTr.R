library(GenomicRanges)
library(Signac)

findOverlappingRanges = function(refStart, refEnd, rangesDf) {
  selRanges = data.frame()
  for ( i in 1:nrow(rangesDf) ) {
    testStart = rangesDf$V2[i]
    testEnd = rangesDf$V3[i]
    if ((refStart <= testEnd) && (testStart <= refEnd)) {
      curRow = rangesDf[i,]
      selRanges = rbind(selRanges,  curRow)
    }
  }
  return(selRanges)
}



df = read.csv("Enpp2_Sox6_genes_peaks_corel.csv")

rangesDf = read.csv("/home/flyhunter/Wang/output/atacRna/peaksGenesCor/atacIntegMacs2_2/pearson/peakGeneCor_ODC_2023-04-13.csv")
rangesDf = rangesDf[rangesDf$seqnames == "chr15",]


rangesDf = read.csv("/home/flyhunter/Wang/output/atacRna/peaksGenesCor/atacIntegMacs2_2/pearson/peakGeneCor_OPC_2023-04-13.csv")
rangesDf = rangesDf[rangesDf$seqnames == "chr7",]


q1 = rangesDf[rangesDf$peak == "chr7-115666776-115667244",]



curRange = StringToGRanges("chr15-54951964-54953128", sep = c("-", "-"))
refRange = StringToGRanges("chr15-54935460-54935838", sep = c("-", "-"))

overlaps <- findOverlaps(curRange, refRange)

overlapping_ranges <- curRange[which(overlapsAny(curRange, refRange))]

(refStart <= testEnd) && (testStart <= refEnd)

(refStart <= testEnd) && (testStart <= refEnd)


rangesDf = read.csv("/home/flyhunter/Wang/output/atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_ODC_50K_2023-10-31.csv")
rangesDf[rangesDf$peak == "chr15-54935460-54935838",]
rangesDf = rangesDf[rangesDf$seqnames == "chr15",]
rangesDf[rangesDf$gene == "Enpp2",]


rangesDf = read.csv("/home/flyhunter/Wang/output/atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_OPC_50K_2023-10-31.csv")
rangesDf[rangesDf$peak == "chr7-115666776-115667244",]
rangesDf = rangesDf[rangesDf$seqnames == "chr7",]
rangesDf[rangesDf$gene == "Sox6",]



#
atacFilt = readRDS("atacIntegrated_macs2_2_RNA")
atacSub = subset(x = atacFilt, subset = Annotations == "ODC")

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC', 'SUB', 'MG', 'C-R')

atacSub<- LinkPeaks(
  object = atacSub ,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = c("Enpp2"), distance = 100000, n_sample = 20
)

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC')
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

p1 = CoveragePlot(
  object = atacSub ,
  region = "Enpp2",
  features = "Enpp2",
  expression.assay = "RNA",
  idents = "ODC",
  extend.upstream = 0,
  extend.downstream = 50000
) 

print(p1)


enpDf = data.frame(Links(atacFilt))

allRnages = row.names(atacSub)
curRanges = StringToGRanges(allRnages, sep = c("-", "-"))
overlapping_ranges <- curRanges[which(overlapsAny(curRanges, refRange))]

subset_seurat_object <- subset(atacSub, features = "chr15-54935460-54935838")

gene1_expr = atacSub@assays$RNA@counts["Enpp2" ,]

atac1_expr = atacSub@assays$Combined_peaks@counts["chr15-54935460-54935838" ,]

cor.test(gene1_expr, atac1_expr, method = "spearman")

cor.test(gene1_expr, atac1_expr, method = "pearson")


non_zero_atac <- which(atac1_expr != 0)
atacNames = names(non_zero_atac)

non_zero_indices <- which(gene1_expr != 0)

matchRna = gene1_expr[atacNames]

cor.test(matchRna, non_zero_atac, method = "spearman")

cor.test(matchRna, non_zero_atac, method = "pearson")


peaks <- CallPeaks(
  object = atacSub,
  effective.genome.size = 1.87e+09,
  macs2.path= '/home/flyhunter/miniconda3/envs/macs2/bin/macs2'
  
)

peaksDf = data.frame(peaks)

chr15 = peaksDf[peaksDf$seqnames == "chr15",]






peaks1 <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
peaks1 <- keepStandardChromosomes(peaks1, pruning.mode = "coarse", species="Mus_musculus")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(atacSub),
  features = peaks1,
  cells = colnames(atacSub)
)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(brainFiltRna) <- annotations


atacSub[["macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(atacSub),
  annotation = annotations
)

# find gene coordinates
subset_gene <- annotations[annotations$gene_name == "Enpp2"]
gene1 = data.frame(subset_gene@ranges)
startPos = gene1$start[1] - 20000
endPos = gene1$end[nrow(gene1)] + 20000
chr = as.character(subset_gene@seqnames@values)

DefaultAssay(atacSub) = "macs2"

target_range <- GRanges(seqnames = chr, ranges = IRanges(start = startPos, end = endPos))
curRanges = StringToGRanges(rownames(atacSub), sep = c("-", "-"))
refRange = StringToGRanges("chr15-54935460-54935838", sep = c("-", "-"))

overlapping_ranges <- curRanges[which(overlapsAny(curRanges, refRange))]

chr15_peaks = paste(chr15$seqnames, chr15$start, chr15$end, sep = "-")

chr15_ranges = StringToGRanges(chr15_peaks, sep = c("-", "-"))



overlapping_ranges <- chr15_ranges[which(overlapsAny(chr15_ranges, target_range))]


findOverlaps(chr15_ranges, target_range)



gene1_expr = atacSub@assays$RNA@counts["Enpp2" ,]

atac1_expr = atacSub@assays$macs2@counts["chr15-54935494-54935831" ,]

cor.test(gene1_expr, atac1_expr, method = "spearman")

cor.test(gene1_expr, atac1_expr, method = "pearson")




gene1_expr = atacFilt@assays$RNA@counts["Enpp2" ,]

atac1_expr = atacFilt@assays$Combined_peaks@counts["chr15-54935460-54935838" ,]

cor.test(gene1_expr, atac1_expr, method = "spearman")

cor.test(gene1_expr, atac1_expr, method = "pearson")





# all clusters 
df = read.csv("Enpp2_Sox6_genes_peaks_corel.csv")
corDf = read.csv("atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/spearman/peakGeneCor_AllClusters_500K_2023-11-06.csv")

corDf[corDf$peak == "chr7-115666776-115667244",]
corDf[corDf$peak == "chr15-54935460-54935838",]