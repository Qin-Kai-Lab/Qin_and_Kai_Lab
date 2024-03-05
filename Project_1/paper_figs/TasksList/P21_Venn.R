#library(ggplot2)
library(ggvenn)
library(ggpubr)

targDir <- 'OPC_ODC/Monocle3/86PC/'

atac = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_TimeCor_ATAC_2023-12-15.csv'))

rna = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv'))

# rna = rna[rna$q_value < 0.01,]
# atac = atac[atac$q_value < 0.01,]

# rna = rna[rna$p_value < 0.05,]
# atac = atac[atac$p_value < 0.05,]

# Positive genes

rnaFilt = rna[rna$estimate > 0,]
atacFilt = atac[atac$estimate > 0,]

genes = list("RNA" = unique(rnaFilt$gene_id), "ATAC" = unique(atacFilt$gene_id))

vennPlot <-ggvenn(genes, text_size=5, set_name_size=8)+
  theme(text = element_text(size = 26)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -6))

vennPlot 

length(atacFilt$gene_id[atacFilt$gene_id%in%rnaFilt$gene_id])

targDir = "Paper_figs/TasksList/P21/"
dir.create(targDir, recursive = T, showWarnings = F)

png(filename = paste0(targDir,"TimeExpression_RnaAtac_Venn_Positive_Pval0.05_2024-01-02.png"), width = 9, height = 9, units = "in", res = 300)
print(vennPlot)
dev.off()

# negative genes

rnaFilt = rna[rna$estimate < 0,]
atacFilt = atac[atac$estimate < 0,]

genes = list("RNA" = unique(rnaFilt$gene_id), "ATAC" = unique(atacFilt$gene_id))

vennPlot <-ggvenn(genes, text_size=5, set_name_size=8)+
  theme(text = element_text(size = 26)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -6))

vennPlot 

length(atacFilt$gene_id[atacFilt$gene_id%in%rnaFilt$gene_id])

targDir = "Paper_figs/TasksList/P21/"
dir.create(targDir, recursive = T, showWarnings = F)

png(filename = paste0(targDir,"TimeExpression_RnaAtac_Venn_Negative_Pval0.05_2024-01-02.png"), width = 9, height = 9, units = "in", res = 300)
print(vennPlot)
dev.off()