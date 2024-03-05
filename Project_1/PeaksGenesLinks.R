library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(karyoploteR)


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
  genes.use = c("Enpp2","Sox6")
)

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC')
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

p1 = CoveragePlot(
  object = atacFilt ,
  region = "Sox6",
  features = "Sox6",
  expression.assay = "RNA",
  idents = clusters,
  extend.upstream = 20000,
  extend.downstream = 20000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p1)

# outFile<-paste0("CovPlot_AtacRna_OPC_Sox6.png")
# ggsave(outFile, plot = p1, height = 12, width = 16, units = 'in', dpi = 300)



p2 = CoveragePlot(
  object = atacFilt ,
  region = "Enpp2",
  features = "Enpp2",
  expression.assay = "RNA",
  idents = clusters,
  extend.upstream = 20000,
  extend.downstream = 20000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p2)

outFile<-paste0("CovPlot_AtacRna_ODC_Enpp2.jpeg")
ggsave(outFile, plot = p2, height = 12, width = 16, units = 'in', dpi = 300)


enpDf = data.frame(Links(atacFilt))
sox6 = data.frame(Links(atacFilt))

combPeaks = data.frame(Links(atacFilt))
selpeaks = combPeaks[c(1,4),]

#write.csv(selpeaks, "Enpp2_Sox6_genes_peaks_corel.csv", row.names = F)

# get sequence 
getSequence = function(selpeaks, i, extSize) {
  curRange = selpeaks$peak[i]
  rangeStr = strsplit(curRange, "-")[[1]]
  curStart = as.numeric(rangeStr[2]) - extSize
  curEnd = as.numeric(rangeStr[3]) + extSize
  newRange =  paste0(rangeStr[1], ":", curStart, "-", curEnd)
  gr <- toGRanges(newRange)
  seq <- as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)[[1]])
  return(seq)
}

sox6Str = getSequence(selpeaks=selpeaks, i=1, extSize=100)
Enpp2Str = getSequence(selpeaks=selpeaks, i=2, extSize=100)


write.table(sox6Str, "Sox6_peak_sequence_2023-10-13.txt", row.names = F, col.names = F, quote = F)
write.table(Enpp2Str, "Enpp2_peak_sequence_2023-10-13.txt", row.names = F, col.names = F, quote = F)

