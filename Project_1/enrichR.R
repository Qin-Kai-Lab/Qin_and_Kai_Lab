library(enrichR)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

curDate<-Sys.Date()

df<-read.csv('./DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv')

dfFilt<-df[(df$p_val_adj < 0.05),]

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

databases<-read.table('../data/enrichR_databases.txt', F)
databases<-as.character(databases$V1)

#databases<-"KEGG_2019_Mouse"

# make directories

makeDirs<-function() {
  #dir.create('./Enrichment/enrichR')
  for (i in databases){
  dirPath<-paste0('./Enrichment/enrichR/', i)
  dir.create((dirPath))
  }
}

makeDirs()


# function for all databases and clusters
totalEnriched<-function(x){
  for (cluster in clusters){
    dfSel<-x[(x$Cell_Type==cluster),]
    if (nrow(dfSel) > 0 ) {
      genes<-unique(dfSel$Genes)
      for (db in databases){
        # enrichment
        enriched <- enrichr(genes = genes, databases = db)[[1]]
        # make a table
        tableFile=paste0('./Enrichment/enrichR/', db, '/', cluster, '_', curDate, ".csv")
        print(tableFile)
        write.csv(enriched, tableFile, row.names = F)
        # filter for plot
        enrichedSign<-enriched[(enriched$P.value < 0.05),]
        if (nrow(enrichedSign) > 0 ) {
          enrichedSort<-enrichedSign[order(enrichedSign$P.value, decreasing = F),]
          enrichedSort$"LogPval"<- (-log10(enrichedSort$P.value))
          topGenes<-head(enrichedSort, 30)
          topGenes$Term <- factor(topGenes$Term, levels = topGenes$Term)
        # make a plot
          plot1<- ggplot(data=topGenes, aes(x= LogPval, y=Term)) +
            geom_bar(stat="identity", fill="#00BFC4")+
            geom_text(aes(label=Term), x=0.1, hjust="left",  color="black", size=6)+
            theme_classic()+
            theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
            labs(x ="-log10(p)")+
            theme(text = element_text(size = 24))+
            scale_y_discrete(limits=rev)
      # save plot
          plotFile=paste0('./Enrichment/enrichR/', db, '/', cluster, '_', curDate, ".jpeg")
          ggsave(plotFile, plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
        }
      }
      
    }
  }
}

totalEnriched(dfFilt)
