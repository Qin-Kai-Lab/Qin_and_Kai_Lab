library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)
library(rtracklayer)

setwd("~/Wang/output")

source('../programs/renameClusters.R')

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv")

curDate<-Sys.Date()

atacInt<-readRDS('atacIntegrated_macs2_2')
DefaultAssay(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacInt) <- annotations
#saveRDS(annotations, "EnsDb_Granges")

#combine rna and atac
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)
rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )
atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

identical(colnames(rnaFilt), colnames(atacFilt))
identical(rownames(rnaFilt@meta.data), rownames(atacFilt@meta.data))
identical(rnaFilt$Merged_CellName, atacFilt$Merged_CellName)

atacFilt[['RNA']]<-rnaFilt@assays$RNA
DefaultAssay(atacInt)<-'Combined_peaks'
atacFilt$Annotations = rnaFilt$Annotations

rm(RNA.combined.norm, rnaFilt, atacInt)
gc()

Idents(atacFilt) = atacFilt$Annotations

#
cond<-'Control_Stress'
targetDir<-paste0('Atac_Rna_coverPlots/macs2_Integ/', cond, '/')
dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(atacFilt))

atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)

table(atacFilt$Annotations)

# if needed ex
clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC', 'SUB', 'MG', 'C-R')

atacFilt<- LinkPeaks(
  object = atacFilt ,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = "Kcnmb2"
)


CoveragePlot(
    object = atacFilt ,
    region = "Kcnmb2",
    features = "Kcnmb2",
    expression.assay = "RNA",
    idents = clusters,
    extend.upstream = 100000,
    extend.downstream = 100000
  ) & theme(text = element_text(size = 22)) 
###

# find gene coordinates
subset_gene <- annotations[annotations$gene_name == "Kcnmb2"]
gene1 = data.frame(subset_gene@ranges)
startPos = gene1$start[1] - 50000
endPos = gene1$end[nrow(gene1)] + 50000
chr = as.character(subset_gene@seqnames@values)

# workw with atacMarix
subset_seurat_object <- subset(atacFilt, Annotations == "GABA")
atacCount = atacFilt@assays$Combined_peaks@counts
atacCount <- atacCount[grep(chr, rownames(atacCount)), ]
allAtacRanges = unique(rownames(atacCount))

target_range <- GRanges(seqnames = chr, ranges = IRanges(start = startPos, end = endPos))
curRanges = StringToGRanges(allAtacRanges, sep = c("-", "-"))

overlapping_ranges <- curRanges[which(overlapsAny(curRanges, target_range))]

## try with fragments table since 19 ranges look suspiciously small
stressFragments = read.table("/home/flyhunter/Wang/data/stress/atac_frag_chr3.tsv", F, sep = "\t" )
stressFragments$V4 = paste0(stressFragments$V4, "_2")
targCells = colnames(subset_seurat_object)
stressFragmentsSub = stressFragments[stressFragments$V4%in%targCells,]
head(stressFragmentsSub)
stress_range_string = paste(stressFragmentsSub$V1, stressFragmentsSub$V2, stressFragmentsSub$V3, sep = "-")
curStressRanges = StringToGRanges(stress_range_string, sep = c("-", "-"))
overlapping_ranges <- curStressRanges[which(overlapsAny(curStressRanges, target_range))]

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

stressOverlapRangeDf = findOverlappingRanges(refStart = startPos, refEnd = endPos, rangesDf=stressFragmentsSub)
nrow(stressOverlapRangeDf)
#non_unique_rows <- stressOverlapRangeDf[duplicated(stressOverlapRangeDf$V2) | duplicated(stressOverlapRangeDf$V2, fromLast = TRUE), ]

fragmentsSum = aggregate(V5~V4, data = stressOverlapRangeDf, FUN = sum)

## gene expression
gene1_expr = subset_seurat_object@assays$RNA@counts["Kcnmb2" ,]
