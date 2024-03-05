library(ggplot2)
library(viridisLite)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output/Enrichment/enrichR/KEGG_2019_Mouse")

curDate<-Sys.Date()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

filesList<-list.files(pattern = '.csv')



getTopPval<-function(x){
  top10Kegg<-data.frame(matrix(ncol=0, nrow=0))
  for ( i in x){
    clusterID<-gsub('\\_.*', '', i)
    df<-read.csv(i)
    dfOrder<-df[order(df$P.value, decreasing = F),]
    topP<-head(dfOrder, 10)
    topP$Cluster<-clusterID
    top10Kegg<-rbind(top10Kegg, topP)
  }
  return(top10Kegg)
}

keggTopP<-getTopPval(filesList)

# add variable present pathways and if not found in cluster give 0


ggplot(keggTopP, aes(y=Term, x=Cluster, fill=P.value))+
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) 