library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
setwd("~/Wang/output")

CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

opcOdc =  readRDS('integ_OPC_ODC')
opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

opcOdc$Cell_ID = rownames(opcOdc@meta.data)

atacFiltSub = subset(atacFilt, Merged_CellName%in%rownames(opcOdc@meta.data))
opcOdcSub = subset(opcOdc, Cell_ID%in%rownames(atacFiltSub@meta.data))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))


atacFiltSub$MonocClust = opcOdc$newMonocClust

# 1
peak = "chr18-31788658-31789364"
curSlot = "counts"
curGene = "Gpr17"
curName =  "Gpr17"

transcripts <- StringToGRanges(peak)

frags <- Fragments(object = atacFiltSub[["Combined_peaks"]])
cells <- colnames(x = atacFiltSub[["Combined_peaks"]])
counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                        process_n = 2000, cells = cells, verbose = T)

atacDat = counts[rownames(x = counts) != "", ]
rnaDat = atacFiltSub@assays$RNA[curSlot][curGene, ]

identical(names(atacDat), names(rnaDat))

curTest = cor.test(rnaDat, atacDat, method = "spearman")
curP = curTest$p.value
curP
# pearson test
curSlot = "data"
rnaDat = atacFiltSub@assays$RNA[curSlot][curGene, ]
cor.test(rnaDat, atacDat, method = "pearson")
cor.test(rnaDat, atacDat, method = "spearman")

# 2
curName =  "Gpr17"
extUp = 2500
extDown = 500

annotation <- Annotation(object = atacFiltSub[["Combined_peaks"]])
curAnnot = annotation[annotation$gene_name == curName,]
transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
curStrand = as.character(transcripts@strand)
if ( curStrand == "-") {
  gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = end(transcripts), end = end(transcripts)), strand = strand("-"), gene_name = curName)
} else {
  gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = start(transcripts), end = start(transcripts)), strand = strand("+"), gene_name = curName)
}

transcripts <- Extend(x = gr, upstream = extUp, 
                      downstream = extDown)

frags <- Fragments(object = atacFiltSub[["Combined_peaks"]])
cells <- colnames(x = atacFiltSub[["Combined_peaks"]])
counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                        process_n = 2000, cells = cells, verbose = T)
gene.key <- transcripts$gene_name
names(x = gene.key) <- GRangesToString(grange = transcripts)
rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
Tsscounts <- counts[rownames(x = counts) != "", ]

identical(names(Tsscounts), names(atacDat))
tssTest = cor.test(Tsscounts,atacDat,  method = "spearman")
tssP = tssTest$p.value
tssP 

# 3 
opcObj = readRDS("OPC_CombPeaks_macs2_2.rds")
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

p5 = CoveragePlot(
  object = opcObj ,
  region = "chr18-32147000-32148243",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  4500,
  extend.downstream = 4500, group.by = "Group", window = 100, annotation = F, assay.scale = "common", show.bulk = T, ymax = 25) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 


p5

# 4
curName = "chr18-32147000-32148243"
transcripts <- StringToGRanges(curName)

frags <- Fragments(object = opcObj[["combined_peaks"]])
cells <- colnames(x = opcObj[["combined_peaks"]])
counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                        process_n = 2000, cells = cells, verbose = T)
counts <- counts[rownames(x = counts) != "", ]
contrGr = opcObj@meta.data$Merge_CellID[opcObj@meta.data$Group == "Control"]
strGr = opcObj@meta.data$Merge_CellID[opcObj@meta.data$Group == "Stress"]
contrDat = counts[names(counts)%in%contrGr]
strDat = counts[names(counts)%in%strGr]
wilcT = wilcox.test(strDat, contrDat, alternative = "greater")
wilcP = wilcT$p.value
wilcP

# 5
peak = "chr18-32147000-32148243"
curSlot = "counts"
curGene = "Gpr17"
curName =  "Gpr17"

transcripts <- StringToGRanges(peak)

frags <- Fragments(object = atacFiltSub[["Combined_peaks"]])
cells <- colnames(x = atacFiltSub[["Combined_peaks"]])
counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                        process_n = 2000, cells = cells, verbose = T)

atacDat = counts[rownames(x = counts) != "", ]
rnaDat = atacFiltSub@assays$RNA[curSlot][curGene, ]

identical(names(atacDat), names(rnaDat))

curTest_1 = cor.test(rnaDat, atacDat, method = "spearman")
curP_1 = curTest_1$p.value
curP_1

# 6  depends on #5 for peak calculations atacDat
curName =  "Gpr17"
extUp = 2500
extDown = 500

annotation <- Annotation(object = atacFiltSub[["Combined_peaks"]])
curAnnot = annotation[annotation$gene_name == curName,]
transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
curStrand = as.character(transcripts@strand)
if ( curStrand == "-") {
  gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = end(transcripts), end = end(transcripts)), strand = strand("-"), gene_name = curName)
} else {
  gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = start(transcripts), end = start(transcripts)), strand = strand("+"), gene_name = curName)
}

transcripts <- Extend(x = gr, upstream = extUp, 
                      downstream = extDown)

frags <- Fragments(object = atacFiltSub[["Combined_peaks"]])
cells <- colnames(x = atacFiltSub[["Combined_peaks"]])
counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                        process_n = 2000, cells = cells, verbose = T)
gene.key <- transcripts$gene_name
names(x = gene.key) <- GRangesToString(grange = transcripts)
rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
Tsscounts <- counts[rownames(x = counts) != "", ]

identical(names(Tsscounts), names(atacDat))
tssTest_1 = cor.test(Tsscounts,atacDat,  method = "spearman")
tssP_1  = tssTest_1 $p.value
tssP_1  

#7
opcObj = readRDS("OPC_CombPeaks_macs2_2.rds")
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")
opcObj <- RegionStats(opcObj, genome = BSgenome.Mmusculus.UCSC.mm10)

opcObj = LinkPeaks(
  opcObj,
  peak.assay = "combined_peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  distance = 500000,
  min.cells = 10,
  genes.use = "Gpr17",
  n_sample = 200,
  pvalue_cutoff = 1,
  score_cutoff = 0,
  gene.id = FALSE,
  verbose = TRUE)

p5 = CoveragePlot(
  object = opcObj,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  0,
  extend.downstream = 2500, group.by = "Group", ymax = 25) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

p5