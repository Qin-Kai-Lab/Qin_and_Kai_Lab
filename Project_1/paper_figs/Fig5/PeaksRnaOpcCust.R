library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(MAST)
library(ggplot2)
library(data.table)

setwd("/home/flyhunter/Wang/output")
opcObj = readRDS("OPC_CombPeaks_macs2_2.rds")
opcObj <- RegionStats(opcObj, genome = BSgenome.Mmusculus.UCSC.mm10)

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

exprMean = AverageExpression(
  opcObj,
  assays = "combined_peaks",
  features = NULL,
  return.seurat = FALSE,
  group.by = "Group",
  add.ident = NULL,
  layer = "counts",
  verbose = TRUE,
)

exprMean = data.frame(exprMean)
strPeaks = exprMean[(exprMean$combined_peaks.Control < exprMean$combined_peaks.Stress),]
strPeaks$Peak = rownames(strPeaks)

# get the reference range
curName = "Gpr17"
annotation <- Annotation(object = opcObj[["combined_peaks"]])
curAnnot = annotation[annotation$gene_name == curName,]
transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = end(transcripts), end = end(transcripts)), strand = strand("-"), gene_name = curName)
transcripts = Extend(gr, upstream = 500000, downstream = 500000)

# select peaks for testing
curPeaks = StringToGRanges(strPeaks$Peak)
overlapping_regions <- subsetByOverlaps(curPeaks,transcripts)

peaksList = GRangesToString(overlapping_regions)

makeGroupComp = function(opcObj, peaksList, curTest, curSlot, minCells, curGroup, curGene) {
  combP = numeric()
  combPeaks = character()
  combEstimate = numeric()
  curObj = subset(opcObj, Group == curGroup)
  for (peak in  peaksList ) {
    atacDat = curObj@assays$combined_peaks["counts"][peak, ]
    rnaDat = curObj@assays$RNA[curSlot][curGene, ]
    cellNumb = length(atacDat[atacDat > 0])
    if (cellNumb >= minCells) {
      try({
        if (curTest=="spearman2") {
          testRes = cor.test(rnaDat, atacDat, method = "spearman")
          pval =  testRes$p.value
          est = testRes$estimate
        } else if (curTest=="spearman1") {
          testRes = cor.test(rnaDat, atacDat, method = "spearman", alternative = "greater")
          pval =  testRes$p.value
          est = testRes$estimate
        } else if (curTest=="pearson2") {
          testRes = cor.test(rnaDat, atacDat, method = "pearson")
          pval =  testRes$p.value
          est = testRes$estimate
        } else if (curTest=="pearson1") {
          testRes = cor.test(rnaDat, atacDat, method = "pearson", alternative = "greater")
          pval =  testRes$p.value
          est = testRes$estimate
        }
        combPeaks = c(combPeaks, peak)
        combP = c(combP, pval)
        combEstimate = c(combEstimate, est)
      })
    }
  }
  combDf = data.frame(Peaks = combPeaks, Estimate = combEstimate, Pval = combP)
  combDf$fdr_p = p.adjust(combDf$Pval, method = "fdr")
  combDf$bonferroni_p =  p.adjust(combDf$Pval, method = "bonferroni")
  combDf = combDf[order(combDf$Pval, decreasing = F),]
  return(combDf)
}

sp2ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman2", curSlot="counts", minCells=3, curGroup="Stress", curGene="Gpr17")
sp2ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman2", curSlot="counts", minCells=3, curGroup="Control", curGene="Gpr17")
sp1ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman1", curSlot="counts", minCells=3, curGroup="Stress", curGene="Gpr17")
sp1ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman1", curSlot="counts", minCells=3, curGroup="Control", curGene="Gpr17")

pr2ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="pearson2", curSlot="data", minCells=3, curGroup="Stress", curGene="Gpr17")
pr2ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="pearson2", curSlot="data", minCells=3, curGroup="Control", curGene="Gpr17")
pr1ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="pearson1", curSlot="data", minCells=3, curGroup="Stress", curGene="Gpr17")
pr1ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="pearson1", curSlot="data", minCells=3, curGroup="Control", curGene="Gpr17")

sp2ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman2", curSlot="data", minCells=3, curGroup="Stress", curGene="Gpr17")
sp2ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman2", curSlot="data", minCells=3, curGroup="Control", curGene="Gpr17")
sp1ResStr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman1", curSlot="data", minCells=3, curGroup="Stress", curGene="Gpr17")
sp1ResContr = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="spearman1", curSlot="data", minCells=3, curGroup="Control", curGene="Gpr17")

strFilt = sp1ResStr[(sp1ResStr$Pval < 0.1 & sp1ResStr$Estimate > 0) ,]
contrFilt = sp1ResContr[(sp1ResContr$Pval < 0.1 & sp1ResContr$Estimate > 0) ,]

finalPeaks = strFilt$Peaks[!strFilt$Peaks%in%contrFilt$Peaks]


opcObj <- RegionStats(opcObj, genome = BSgenome.Mmusculus.UCSC.mm10)
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")
p5 = CoveragePlot(
  object = opcObj ,
  region = "chr18-32237263-32237559",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  50,
  extend.downstream = 50, group.by = "Group", window = 100, annotation = F, assay.scale = "separate", show.bulk = T) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

p5





strFilt = pr1ResStr[(pr1ResStr$Pval < 0.1 & pr1ResStr$Estimate > 0) ,]
contrFilt = pr1ResContr[(pr1ResContr$Pval < 0.1 & pr1ResContr$Estimate > 0) ,]

finalPeaks = strFilt$Peaks[!strFilt$Peaks%in%contrFilt$Peaks]


custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")
p5 = CoveragePlot(
  object = opcObj ,
  region = "chr18-32237263-32237559",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  50,
  extend.downstream = 50, group.by = "Group", window = 100, annotation = F, assay.scale = "separate", show.bulk = T) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

p5



dir.create("Paper_figs/Fig5/DiffExprPeaksTss")


endPost = as.numeric(as.character(gsub("^.*\\-","", t1Res$Peaks)))
endPos1 = 31949636- endPost
t1Res$Distance_TSS = endPos1 
