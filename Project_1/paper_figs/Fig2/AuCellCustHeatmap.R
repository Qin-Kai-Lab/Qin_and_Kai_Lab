library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)

df = read.csv("AUCell_KeggClustProF_All_ContrVsStress_2022-11-28.csv") # all
#df = df[df$avg_log2FC > 0,] # positive
df = df[df$avg_log2FC < 0,] # negative

colnames(df)[1:2] = c("Cluster", "pvalue") 
colnames(df)[7] = c("Description") 

clusters = c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC', "C-R", 'MG')

df = df[df$Cluster%in%clusters,]


# need to rename the columns to match the program, rename final paths and filter based on positvie negative logf2 change at the beginning of the code
# stopped here

getTopPath = function(groupMarkers, topNumb) {
  combDf = data.frame(matrix(nrow = 0, ncol = 0))
  for (curCluster in clusters) {
    curDf = groupMarkers[(groupMarkers$Cluster == curCluster),]
    curDfOrd = curDf[order(curDf$pvalue, decreasing = F),]
    dfTop = head(curDfOrd, topNumb)
    combDf = rbind(combDf, dfTop)
    topMark = unique(combDf$Description)
  }
  return(topMark)
}

topMark = getTopPath(groupMarkers = df, topNumb = 5)

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Description)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  allPAdj = numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    for (i in keggList){
      if (i %nin% dfSel$Description){
        curDesc<-i
        curClust<-cluster
        curP<-1
        curCount<-0
        curPAdj = 1
        allDesc<-c(allDesc, curDesc)
        allClust<-c(allClust, curClust)
        allP<-c(allP, curP)
        allCount<-c(allCount, curCount)
        allPAdj = c(allPAdj, curPAdj)
        # combine all vectors and make a dataframe from them
      }
    }
  }
  finalTable<-data.frame(Description=allDesc, pvalue=allP, p_val_adj=allPAdj, avg_log2FC=allCount, Cluster=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab = df, clusters = clusters)
keggComplete<-rbind.fill(df, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]


addPCat<-function(dataTab, pcol, newCol){
  dataTab[newCol]<-NA
  for (i in 1:nrow(dataTab)){
    if (dataTab[[pcol]][i] > 0.05){
      dataTab[[newCol]][i]<- " > 0.05"
    } else if (dataTab[[pcol]][i] < 0.05 & dataTab[[pcol]][i] >= 0.01  ) {
      dataTab[[newCol]][i]<- " < 0.05"
    } else if (dataTab[[pcol]][i] < 0.01 & dataTab[[pcol]][i] >= 0.001) {
      dataTab[[newCol]][i]<- " < 0.01"
    } else if (dataTab[[pcol]][i] < 0.001) {
      dataTab[[newCol]][i]<- " < 0.001"
    }
  }
  return(dataTab)
}

keggComplete<-addPCat(dataTab=keggComplete, pcol = "p_val_adj", newCol = "P_Adjust")

keggComplete$Cluster = factor(keggComplete$Cluster, levels =  c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG', 'C-R', 'SUB'))

keggComplete$Description = gsub(" - Mus musculus \\(house mouse\\)", "", keggComplete$Description)

curHeatMap = ggplot(keggComplete, aes(y=Description, x=Cluster, fill=P_Adjust))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(discrete=T ) +
  theme(text = element_text(size = 28))+
  theme(legend.text = element_text(size = 26)) +
  theme(legend.title = element_text(size = 26)) +
  theme(axis.text.y = element_text(size = 26))

print(curHeatMap)

targDir = "Paper_figs/Fig2/Heatmap/"
#dir.create(targDir, recursive = T, showWarnings = F)

topN = "5"
typeGenes = "Negative"
png(filename = paste0(targDir,"AuCell_Top_", topN, "_", typeGenes, "_KEGG_2023-11-27.png"), width = 22, height = 16, units = "in", res = 300)
print(curHeatMap)
dev.off()
