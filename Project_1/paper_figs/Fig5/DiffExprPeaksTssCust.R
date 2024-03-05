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
transcripts = Extend(gr, upstream = 200000, downstream = 200000)

# select peaks for testing
curPeaks = StringToGRanges(strPeaks$Peak)
overlapping_regions <- subsetByOverlaps(curPeaks,transcripts)

peaksList = GRangesToString(overlapping_regions)

makeGroupComp = function(opcObj, peaksList, curTest, curSlot, minCells) {
  combP = numeric()
  combPeaks = character()
  for (peak in  peaksList ) {
    contrGroup = opcObj@assays$combined_peaks[curSlot][peak, ][opcObj@meta.data$Group == "Control"]
    strGroup = opcObj@assays$combined_peaks[curSlot][peak, ][opcObj@meta.data$Group == "Stress"]
    contrNumb = length(contrGroup[contrGroup > 0])
    strNumb = length(strGroup[strGroup  > 0])
    if (contrNumb >=minCells | strNumb >= minCells) {
      try({
        if (curTest=="wilcox") {
          testRes = wilcox.test(contrGroup, strGroup)
          pval =  testRes$p.value
        } else if (curTest=="wilcox1") {
          testRes = wilcox.test(contrGroup, strGroup, alternative = "less")
          pval =  testRes$p.value
          } else if (curTest=="ttest2") {
          testRes = t.test(contrGroup, strGroup)
          pval =  testRes$p.value
        } else if (curTest=="ttest1") {
          testRes = t.test(contrGroup, strGroup, alternative = "less")
          pval =  testRes$p.value
        }
        combPeaks = c(combPeaks, peak)
        combP = c(combP, pval)
      })
    }
  }
  combDf = data.frame(Peaks = combPeaks, Pval = combP)
  combDf$fdr_p = p.adjust(combDf$Pval, method = "fdr")
  combDf$bonferroni_p =  p.adjust(combDf$Pval, method = "bonferroni")
  combDf = combDf[order(combDf$Pval, decreasing = F),]
  return(combDf)
}

wilcRes = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="wilcox", curSlot="counts", minCells = 20)
wilc1Res = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="wilcox1", curSlot="counts", minCells = 20)
t2Res = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="ttest2", curSlot="data", minCells = 20)
t1Res = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="ttest1", curSlot="data", minCells = 20)

dir.create("Paper_figs/Fig5/DiffExprPeaksTss")


endPost = as.numeric(as.character(gsub("^.*\\-","", t1Res$Peaks)))
endPos1 = 31949636- endPost
t1Res$Distance_TSS = endPos1 

write.csv(t1Res, "Paper_figs/Fig5/DiffExprPeaksTss/StrHigherPeaks_Gpr17_Tss200KB_MinC3_2024-02-22.csv", row.names = F)

t1Res = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="ttest1", curSlot="data", minCells = 0)
endPost = as.numeric(as.character(gsub("^.*\\-","", t1Res$Peaks)))
endPos1 = 31949636- endPost
t1Res$Distance_TSS = endPos1 

#write.csv(t1Res, "Paper_figs/Fig5/DiffExprPeaksTss/StrHigherPeaks_Gpr17_Tss100KB_MinC0_2024-02-21.csv", row.names = F)

# make test with all peaks
peaksList = strPeaks$Peak
wilcResAll = makeGroupComp(opcObj=opcObj, peaksList=peaksList, curTest="wilcox", curSlot="counts", minCells = 5)
write.csv(wilcResAll, "Paper_figs/Fig5/DiffExprPeaksTss/StrHigherPeaks_All_MinC5_2024-02-21.csv", row.names = F)

## check that my subsampling is correct

# mySample = opcObj@assays$RNA[curSlot][curName, ][opcObj@meta.data$Group == "Stress"]
# strObj = subset(opcObj, Group == "Stress")
# defMethod = strObj@assays$RNA@counts[curName, ]
# identical(mySample, defMethod)
# 
# mySample = opcObj@assays$RNA[curSlot][curName, ][opcObj@meta.data$Group == "Stress"]
# strObj = subset(opcObj, Group == "Stress")
# defMethod = strObj@assays$RNA@data[curName, ]
# identical(mySample, defMethod)