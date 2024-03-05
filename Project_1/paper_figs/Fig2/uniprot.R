library(biomaRt)

targDir = './Paper_figs/Fig2/'
curDate<-Sys.Date()

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

# Connect to the appropriate biomart database (e.g., Ensembl)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Convert gene symbols to UniProt protein IDs
gene_symbols <- unique(groupMarkers$Genes)  # Replace with your gene symbols of interest


convert <- getBM(attributes = c("mgi_symbol", "uniprotswissprot", "uniprotsptrembl"), 
                 filters = "mgi_symbol", 
                 values = gene_symbols, 
                 mart = ensembl)


colnames(convert)[1] = "Genes"

annot = plyr::join(groupMarkers, convert, type = "left", by = "Genes", match = "all")
write.csv(annot, paste0(targDir, 'allDiffExprLogfc0.25_ContrVsStress_Uniprot_', curDate, ".csv"), row.names = F)

## uniprot for all clusters together
source('../programs/renameClusters.R')

library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)
library(biomaRt)
library(MAST)

targDir = './Paper_figs/Fig2/'
dir.create(targDir, recursive = T, showWarnings = F)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

levels(RNA.combined.norm)

Idents(RNA.combined.norm) = RNA.combined.norm$group

groupMarkers = FindMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", ident.1 = "Stress", ident.2 = "Control")
groupMarkers$Genes = rownames(groupMarkers)

gene_symbols <- unique(groupMarkers$Genes)  # Replace with your gene symbols of interest

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
convert <- getBM(attributes = c("mgi_symbol", "uniprotswissprot", "uniprotsptrembl"), 
                 filters = "mgi_symbol", 
                 values = gene_symbols, 
                 mart = ensembl)


colnames(convert)[1] = "Genes"

annot = plyr::join(groupMarkers, convert, type = "left", by = "Genes", match = "all")
write.csv(annot, paste0(targDir, 'DiffExprClustComb_ContrVsStress_Uniprot_', curDate, ".csv"), row.names = F)
