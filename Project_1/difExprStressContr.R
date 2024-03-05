library(MAST)
library(limma)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

levels(RNA.combined.norm)

DefaultAssay(RNA.combined.norm)

# counts per gene
#r1<-data.frame(rowSums(RNA.combined.norm@assays$RNA@counts))

#
# function to claculate differential gene expression
findMarkersGr<-function(dat, clust, pos){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$group
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", ident.1 = "Stress", ident.2 = "Control")
    allMarkers$Genes<-rownames(allMarkers)
    write.csv(allMarkers, paste0('DiffExprLogfc0.25_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F)
write.csv(groupMarkers, paste0('allDiffExprLogfc0.25_ContrVsStress_',curDate, '.csv'), row.names = F)

rm(groupMarkers)

#  minimum cells per cluster filtering provides the same results

# use wilcoxon test
findMarkersGr<-function(dat, clust, pos){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$group
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox", 
                              min.cells.feature = 20, ident.1 = "Stress", ident.2 = "Control")
    allMarkers$Genes<-rownames(allMarkers)
    write.csv(allMarkers, paste0('DiffExprLogfc0.25_Wilcox_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F)
write.csv(groupMarkers, paste0('allDiffExprLogfc0.25_Wilcox_ContrVsStress_',curDate, '.csv'), row.names = F)

# summarizw number of genes that per log fold change group
addLogGr<-function(x){
  deGenes<-read.csv(x)
  deGenesFilt<-deGenes[(deGenes$p_val<0.05),]
  deGenesFilt$AbsLog<-abs(deGenesFilt$avg_log2FC)
  deGenesFilt$FoldCh<-(2^deGenesFilt$AbsLog)
  deGenesFilt$FoldChGr<-NA
  for (i in 1:nrow(deGenesFilt)) {
    if (deGenesFilt$FoldCh[i] >= 1.25 & deGenesFilt$FoldCh[i] < 1.5) {
      deGenesFilt$FoldChGr[i] = 1.25
    } else if (deGenesFilt$FoldCh[i] >= 1.5 & deGenesFilt$FoldCh[i] < 2) {
      deGenesFilt$FoldChGr[i] = 1.5
    } else if ( deGenesFilt$FoldCh[i] >= 2 ) {
      deGenesFilt$FoldChGr[i] = '>=2'
    } else if ( deGenesFilt$FoldCh[i] < 1.25 ) {
      deGenesFilt$FoldChGr[i] = '<1.25' 
    }
    
  }
  return(deGenesFilt)
}

# make a summary
sumLogF<-function(x){
  sumLog<-data.frame(table(x$Cell_Type, x$FoldChGr))
  colnames(sumLog)<-c("Cell_type", "Fold_change_group", "Freq")
  return(sumLog)
}

# default
deGenesDef<-addLogGr(x='DE_StressVsControl/default/allDiffExpr0.2ResGroups_2022-10-24.csv')
deGenesDefSum<-sumLogF(deGenesDef)

deGenesDefSumTot<-data.frame(table(deGenesDef$FoldChGr))
colnames(deGenesDefSumTot)[1]<-'Cell_type'

curDate<-Sys.Date()

write.csv(deGenesDefSum, paste0('genesPerCoefGroup_cellType_', curDate, '.csv'), row.names = F )

write.csv(deGenesDefSumTot, paste0('genesPerCoefGroup_', curDate, '.csv'), row.names = F )

rm(deGenesDef, deGenesDefSum, deGenesDefSumTot)
# filtered

deGenesDef<-addLogGr(x='DE_StressVsControl/min20CellsGroup/allDiffExpr0.2ResGroups_2022-10-24_min20Cells.csv')
deGenesDefSum<-sumLogF(deGenesDef)

deGenesDefSumTot<-data.frame(table(deGenesDef$FoldChGr))
colnames(deGenesDefSumTot)[1]<-'Cell_type'

curDate<-Sys.Date()

write.csv(deGenesDefSum, paste0('genesPerCoefGroup_cellType_min20Cells_', curDate, '.csv'), row.names = F )

write.csv(deGenesDefSumTot, paste0('genesPerCoefGroup_min20Cells_', curDate, '.csv'), row.names = F )

rm(deGenesDef, deGenesDefSum, deGenesDefSumTot)