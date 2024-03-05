library(MAST)
library(dplyr)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

source('../programs/renameClusters.R')

targDir = './Paper_figs/Fig2/'

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
findMarkersGr<-function(dat, clust, pos){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$group
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0, logfc.threshold = 0, test.use = "MAST", ident.1 = "Stress", ident.2 = "Control")
    allMarkers$Genes<-rownames(allMarkers)
    #write.csv(allMarkers, paste0('DiffExprAll_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F)

sumTab = data.frame(table(groupMarkers$Cell_Type))

groupMarkers$AbsLog = abs(groupMarkers$avg_log2FC)

write.csv(groupMarkers, paste0(targDir, "allDiffExpr_ContrVsStress_.csv"), row.names = F)

avgCl =  groupMarkers %>%
  group_by(Cell_Type) %>%
  summarize(avg_log2FC = mean(AbsLog))

write.csv(avgCl, paste0(targDir, "avgLogfDiff_ContrVsStress_AllGenes_.csv"), row.names = F)

# filtered markers 
groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
groupMarkers$AbsLog = abs(groupMarkers$avg_log2FC)
groupMarkers = groupMarkers[(groupMarkers$p_val_adj < 0.05),]

avgCl =  groupMarkers %>%
  group_by(Cell_Type) %>%
  summarize(avg_log2FC = mean(AbsLog))

write.csv(avgCl, paste0(targDir, "avgLogfDiff_ContrVsStress_SignGenes_.csv"), row.names = F)

## alternative 
gene_names <- rownames(RNA.combined.norm@assays$RNA)
gene_index <- which(gene_names == "Malat2")
gene_expression <- RNA.combined.norm@assays$RNA@data[gene_index, ]

# calculate the distance
gene_expression <- as.matrix(RNA.combined.norm@assays$RNA@data)

gene_expression[1,1]
colnames(gene_expression)[1]
rownames(gene_expression)[1]
