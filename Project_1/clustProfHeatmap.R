library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output/Enrichment/clusterProfiler/kegg")

curDate<-Sys.Date()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

filesList<-list.files(pattern = '.csv')

getTopPval<-function(x){
  top10Kegg<-data.frame(matrix(ncol=0, nrow=0))
  for ( i in x){
    clusterID<-gsub('\\_.*', '', i)
    df<-read.csv(i)
    dfOrder<-df[order(df$p.adjust, decreasing = F),]
    topP<-head(dfOrder, 10)
    topP$Cluster<-clusterID
    top10Kegg<-rbind(top10Kegg, topP)
  }
  return(top10Kegg)
}

keggTopP<-getTopPval(filesList)

# add variable present pathways and if not found in cluster give 0



findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Description)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    for (i in keggList){
      if (i %nin% dfSel$Description){
        curDesc<-i
        curClust<-cluster
        curP<-1
        curCount<-0
        allDesc<-c(allDesc, curDesc)
        allClust<-c(allClust, curClust)
        allP<-c(allP, curP)
        allCount<-c(allCount, curCount)
        # combine all vectors and make a dataframe from them
      }
    }
  }
  finalTable<-data.frame(Description=allDesc, pvalue=allP, Count=allCount, Cluster=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab =  keggTopP, clusters = clusters)

keggComplete<-rbind.fill(keggTopP, missDat)


addPCat<-function(dataTab){
  dataTab$P.value<-NA
  for (i in 1:nrow(dataTab)){
    if (dataTab$pvalue[i] > 0.05){
      dataTab$P.value[i]<- " > 0.05"
    } else if (dataTab$pvalue[i] < 0.05 & dataTab$pvalue[i] >= 0.01  ) {
      dataTab$P.value[i]<- " < 0.05"
    } else if (dataTab$pvalue[i] < 0.01 & dataTab$pvalue[i] >= 0.001) {
      dataTab$P.value[i]<- " < 0.01"
    } else if (dataTab$pvalue[i] < 0.001) {
      dataTab$P.value[i]<- " < 0.001"
    }
  }
  return(dataTab)
}

keggComplete<-addPCat(keggComplete)
keggComplete$Cluster <- factor(keggComplete$Cluster, levels=clusters)

# heatmap with p vlaues
ggplot(keggComplete, aes(y=Description, x=Cluster, fill=P.value))+
  geom_tile() +
  scale_fill_viridis(discrete=T ) +
  theme_classic()+
  theme(text = element_text(size = 16))

ggsave(paste0('HeatMap_Top10_', curDate, '.jpeg'), height = 16, width = 20, units = 'in', dpi = 300)

# heatmap with Counts
ggplot(keggComplete, aes(y=Description, x=Cluster, fill=Count))+
  geom_tile() +
  scale_fill_viridis(discrete=F ) +
  theme_classic()+
  theme(text = element_text(size = 16))

ggsave(paste0('HeatMap_Top10Count_', curDate, '.jpeg'), height = 16, width = 20, units = 'in', dpi = 300)


rm(keggTopP, missDat, keggComplete)

# all kegg Pathways

getAllPval<-function(x){
  allKegg<-data.frame(matrix(ncol=0, nrow=0))
  for ( i in x){
    clusterID<-gsub('\\_.*', '', i)
    df<-read.csv(i)
    dfOrder<-df[order(df$p.adjust, decreasing = F),]
    topP<-dfOrder
    topP$Cluster<-clusterID
    allKegg<-rbind(allKegg, topP)
  }
  return(allKegg)
}

keggAllP<-getAllPval(filesList)

missDatAll<-findMissing(dataTab =  keggAllP, clusters = clusters)

keggCompleteAll<-rbind.fill(keggAllP, missDatAll)

keggCompleteAll<-addPCat(keggCompleteAll)
keggCompleteAll$Cluster <- factor(keggCompleteAll$Cluster, levels=clusters)

# heatmap P-values

ggplot(keggCompleteAll, aes(y=Description, x=Cluster, fill=P.value))+
  geom_tile() +
  scale_fill_viridis(discrete=T ) +
  theme_classic()+
  theme(text = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 8)) 

ggsave(paste0('HeatMap_All_', curDate, '.jpeg'), height = 26, width = 20, units = 'in', dpi = 300)

# heatmap Counts
ggplot(keggCompleteAll, aes(y=Description, x=Cluster, fill=Count))+
  geom_tile() +
  scale_fill_viridis(discrete=F ) +
  theme_classic()+
  theme(text = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 8)) 

ggsave(paste0('HeatMap_AllCount_', curDate, '.jpeg'), height = 26, width = 20, units = 'in', dpi = 300)
