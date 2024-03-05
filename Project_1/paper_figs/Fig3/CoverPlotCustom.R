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

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")


clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC', 'SUB', 'MG', 'C-R')

atacFilt<- LinkPeaks(
  object = atacFilt ,
  peak.assay = "Combined_peaks",
  expression.assay = "RNA",
  genes.use = c('Sema3e','Cadm1', 'Gpr17')
)

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC')
custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

p1 = CoveragePlot(
  object = atacFilt ,
  region = "Sema3e",
  features = "Sema3e",
  expression.assay = "RNA",
  idents = clusters,
  extend.upstream = 20000,
  extend.downstream = 20000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p1)

p2 = CoveragePlot(
  object = atacFilt ,
  region = "Cadm1",
  features = "Cadm1",
  expression.assay = "RNA",
  idents = clusters,
  extend.upstream = 20000,
  extend.downstream = 20000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p2)

p3 = CoveragePlot(
  object = atacFilt ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  idents = clusters,
  extend.upstream = 20000,
  extend.downstream = 20000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p3)

targDir = './Paper_figs/Fig3/'

ggsave(paste0(targDir, "CovPlot_Sema3e_2024-02-06.png"), plot = p1, height = 12, width = 16, units = 'in', dpi = 300)

ggsave(paste0(targDir, "CovPlot_Cadm1_2024-02-06.png"), plot = p2, height = 12, width = 16, units = 'in', dpi = 300)

ggsave(paste0(targDir, "CovPlot_Gpr17_2024-02-08.png"), plot = p2, height = 12, width = 16, units = 'in', dpi = 300)