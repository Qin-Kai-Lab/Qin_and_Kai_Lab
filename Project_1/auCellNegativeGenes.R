
inDir<-'./Enrichment/AUcell/Kegg_ClustProfiler/pos_and_neg/'

curDate<-Sys.Date()

targetDir<-'./Enrichment/AUcell/Kegg_ClustProfiler/negative/'

dir.create(targetDir, recursive = T)



files.list<-list.files(inDir, pattern = '.csv')

for ( i in files.list){
  targFile<-paste0(inDir, i)
  df<-read.csv(targFile)
  dfSel<-df[(df$avg_log2FC < 0),]
  outPath<-paste0(targetDir, i)
  write.csv(dfSel, outPath, row.names = F)
}

groupMarkers<-read.csv("./Enrichment/AUcell/Kegg_ClustProfiler/negative/All_ContrVsStress_2022-11-28.csv")


# make heatmaps

getTopPval<-function(dataTab, clusters){
  top10Kegg<-data.frame(matrix(ncol=0, nrow=0))
  df<-dataTab[(dataTab$p_val_adj < 0.05),]
  for ( i in clusters){
    dfClust<-df[(df$Cell_Type==i),]
    dfClust$absLog<-abs(dfClust$avg_log2FC)
    dfOrder<-dfClust[order(dfClust$absLog, decreasing = T),]
    topP<-head(dfOrder, 10)
    top10Kegg<-rbind(top10Kegg, topP)
  }
  return(top10Kegg)
}

keggTopP<-getTopPval(dataTab = groupMarkers, clusters = clusters)

###

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Genes)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
    for (i in keggList){
      if (i %nin% dfSel$Genes){
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
  finalTable<-data.frame(Genes=allDesc, p_val=allP, p_val_adj=allP, avg_log2FC=allCount, Cell_Type=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab =  keggTopP, clusters = clusters)

keggComplete<-rbind.fill(keggTopP, missDat)

###

addPCat<-function(dataTab){
  dataTab$P.value<-NA
  for (i in 1:nrow(dataTab)){
    if (dataTab$p_val_adj[i] > 0.05){
      dataTab$P.value[i]<- " > 0.05"
    } else if (dataTab$p_val_adj[i] < 0.05 & dataTab$p_val_adj[i] >= 0.01  ) {
      dataTab$P.value[i]<- " < 0.05"
    } else if (dataTab$p_val_adj[i] < 0.01 & dataTab$p_val_adj[i] >= 0.001) {
      dataTab$P.value[i]<- " < 0.01"
    } else if (dataTab$p_val_adj[i] < 0.001) {
      dataTab$P.value[i]<- " < 0.001"
    }
  }
  return(dataTab)
}

keggComplete<-addPCat(keggComplete)
keggComplete$Cluster <- factor(keggComplete$Cell_Type, levels=clusters)



# heatmap with p vlaues
ggplot(keggComplete, aes(y=Genes, x=Cell_Type, fill=P.value))+
  geom_tile() +
  scale_fill_viridis(discrete=T ) +
  theme_classic()+
  theme(text = element_text(size = 16))

ggsave(paste0(targetDir, 'HeatMap_Top10Pval_StressVsContr', curDate, '.jpeg'), height = 16, width = 20, units = 'in', dpi = 300)

# heatmap with Counts
ggplot(keggComplete, aes(y=Genes, x=Cell_Type, fill=avg_log2FC))+
  geom_tile() +
  scale_fill_viridis(discrete=F ) +
  theme_classic()+
  theme(text = element_text(size = 16))

ggsave(paste0(targetDir, 'HeatMap_TopLog2F_StressVsContr_', curDate, '.jpeg'), height = 16, width = 20, units = 'in', dpi = 300)