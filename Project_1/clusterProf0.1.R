library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

curDate<-Sys.Date()

df<-read.csv('./DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25_Ens.csv')

dfFilt<-df[(df$p_val_adj < 0.05),]

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

entrez <- as.data.frame(mapIds(org.Mm.eg.db, keys = dfFilt$Genes, keytype = "SYMBOL", column="ENTREZID"))

colnames(entrez)<-'entrez'

dfEntr<-cbind(dfFilt, entrez)

# function kegg enrichment

targetDir<-'./Enrichment/clusterProfiler/'

dir.create(file.path(targetDir), recursive = T, showWarnings = FALSE)


allEnrichedKegg<-function(dataTab, clusters){
  for ( cluster in clusters){
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
    if (nrow(dfSel) > 0) {
      genes<-unique(dfSel$entrez)
      
      ego2 <- enrichKEGG(gene = genes, organism = 'mmu',
                         pvalueCutoff = 0.05)
      
      ego2Res<-ego2@result
      
      newdir<-paste0(targetDir, 'kegg', '/')
      dir.create(newdir, showWarnings = FALSE)
      
      tablePath<-paste0(newdir, cluster, '_', curDate, '.csv')
      write.csv(ego2Res, tablePath, row.names = F)
      
       if (nrow(ego2Res[(ego2Res$p.adjust < 0.05),]) > 0){
         plot1<-barplot(ego2, showCategory=20)
         plotPath<-paste0(newdir, cluster, '_20Terms_', curDate, '.jpeg')
         ggsave(plotPath, plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
         
       }
      
      
    }
    
  }
  
}

allEnrichedKegg(dataTab=dfEntr, clusters = clusters)

# go analysis
targetDir<-'./Enrichment/clusterProfiler/GO_entrez/Ontologies/'

dir.create(file.path(targetDir), recursive = T, showWarnings = FALSE)

allEnriched<-function(dataTab, clusters, ont, pool){
  for ( cluster in clusters){
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
    if (nrow(dfSel) > 0) {
      genes<-unique(dfSel$entrez)
      
      ego2 <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05, keyType = "ENTREZID", ont = ont, pool = pool)
      
      ego2Res<-ego2@result
      
      newdir<-paste0(targetDir, ont, '/')
      dir.create(newdir, showWarnings = FALSE)
      
      tablePath<-paste0(newdir, cluster, '_', curDate, '.csv')
      write.csv(ego2Res, tablePath, row.names = F)
      
      
      if (nrow(ego2Res[(ego2Res$p.adjust < 0.05),]) > 0){
        plot1<-barplot(ego2, showCategory=20)
        plotPath<-paste0(newdir, cluster, '_20Terms_', curDate, '.jpeg')
        ggsave(plotPath, plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
      }

      
    }
    
  }
  
}

allEnriched(dataTab=dfEntr, clusters = clusters, ont='ALL',  pool = T)

allEnriched(dataTab=dfEntr, clusters = clusters, ont='CC', pool = F)

allEnriched(dataTab=dfEntr, clusters = clusters, ont='MF',  pool = F)

allEnriched(dataTab=dfEntr, clusters = clusters, ont='BP',  pool = F)


