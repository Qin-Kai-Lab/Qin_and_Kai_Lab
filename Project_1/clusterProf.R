library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

curDate<-Sys.Date()

df<-read.csv('./DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25_Ens.csv')

dfFilt<-df[(df$p_val_adj < 0.05),]

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

targetDir<-'./Enrichment/clusterProfiler/GO/Ontologies/'

dir.create(file.path(targetDir), recursive = T)


# function

allEnriched<-function(dataTab, clusters, ont, pool){
  for ( cluster in clusters){
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
    if (nrow(dfSel) > 0) {
      genes<-unique(dfSel$Ens)
      
      ego2 <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05, keyType = "ENSEMBL", ont = ont, pool = pool)
      
      ego2Res<-ego2@result
      
      newdir<-paste0(targetDir, ont, '/')
      dir.create(newdir)
      
      tablePath<-paste0(newdir, cluster, '_', curDate, '.csv')
      write.csv(ego2Res, tablePath, row.names = F)
      
      plot1<-barplot(ego2, showCategory=20)
      plotPath<-paste0(newdir, cluster, '_20Terms_', curDate, '.jpeg')
      ggsave(plotPath, plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
      
    }
    
  }
  
}

allEnriched(dataTab=dfFilt, clusters = clusters, ont='ALL',  pool = T)

allEnriched(dataTab=dfFilt, clusters = clusters, ont='CC', pool = F)

allEnriched(dataTab=dfFilt, clusters = clusters, ont='MF',  pool = F)

allEnriched(dataTab=dfFilt, clusters = clusters, ont='BP',  pool = F)

