source('../programs/renameClusters.R')

library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)

targDir = './Paper_figs/Fig2/DotPlot/'
dir.create(targDir, recursive = T, showWarnings = F)
clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
curDate<-Sys.Date()
levels(RNA.combined.norm)
DefaultAssay(RNA.combined.norm)

getTopMarkers = function(df, topNumb) {
  clusters = unique(df$Cell_Type)
  markers = character()
  dfPct = df[(df$pct.1 >0.25) | (df$pct.2 >0.25), ]
  for (cluster in clusters) {
    dfSub =  dfPct[(dfPct[['Cell_Type']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$AbsLog, decreasing = T),]
    topMark = head(dfOrd$Genes, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}


# find missing
findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Description)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
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
  colnames(finalTable) = c("Description", "p_val", "avg_log2FC", "Cell_Type")
  return(finalTable)
}

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
topMark = getTopMarkers(df = groupMarkers, topNumb = 10)

colnames(groupMarkers)[7] = "Description"
#topMarkDf = groupMarkers[(groupMarkers$Description%in%topMark),]

clusters = unique(groupMarkers$Cell_Type)
missDat<-findMissing(dataTab =  groupMarkers, clusters = clusters)
keggComplete<-rbind.fill(groupMarkers, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]
keggComplete$p_val_adj[is.na(keggComplete$p_val_adj)] = 1

keggComplete$PvalMod = keggComplete$p_val_adj
keggComplete$PvalMod[keggComplete$PvalMod == 0] <-  1e-300
keggComplete$`-Log10(Pval_Adj)` = log10(keggComplete$PvalMod) * -1
keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cell_Type, size = `-Log10(Pval_Adj)`, color = avg_log2FC)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 10))+ 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  xlab("Genes") +
  theme(text = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =18, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top10Genes_RNA_', curDate, '.png'), height = 16, width = 26, units = 'in', res = 300)
print(custDotPlot)
dev.off()


### KEGG
groupMarkers = read.csv("AllClust_ClustProF_KEGG_2023-05-04.csv")
nrow(groupMarkers)
groupMarkers = groupMarkers[!(groupMarkers$Cluster == "AllClust_ClustProF_KEGG_2023-05-04.csv"),]
nrow(groupMarkers)

getTopPath = function(groupMarkers, topNumb) {
  combDf = data.frame(matrix(nrow = 0, ncol = 0))
  for (curCluster in clusters) {
    curDf = groupMarkers[(groupMarkers$Cluster == curCluster),]
    curDfOrd = curDf[order(curDf$p.adjust, decreasing = F),]
    dfTop = head(curDfOrd, topNumb)
    combDf = rbind(combDf, dfTop)
    topMark = unique(combDf$Description)
  }
  return(topMark)
}

topMark = getTopPath(groupMarkers = groupMarkers, topNumb = 10)

#topMarkDf = groupMarkers[(groupMarkers$Description%in%topMark),]

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Description)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  allGeneRatio <- character()
  allPAdj = numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    for (i in keggList){
      if (i %nin% dfSel$Description){
        curDesc<-i
        curClust<-cluster
        curP<-1
        curCount<-0
        curGeneRatio = "0/0"
        curPAdj = 1
        allDesc<-c(allDesc, curDesc)
        allClust<-c(allClust, curClust)
        allP<-c(allP, curP)
        allCount<-c(allCount, curCount)
        allGeneRatio = c(allGeneRatio, curGeneRatio)
        allPAdj = c(allPAdj, curPAdj)
        # combine all vectors and make a dataframe from them
      }
    }
  }
  finalTable<-data.frame(Description=allDesc, pvalue=allP, Count=allCount, Cluster=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab = groupMarkers, clusters = clusters)

keggComplete<-rbind.fill(groupMarkers, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]

addPCat<-function(dataTab){
  dataTab$Pval_Adj<-NA
  for (i in 1:nrow(dataTab)){
    if (dataTab$p.adjust[i] > 0.05){
      dataTab$Pval_Adj[i]<- " > 0.05"
    } else if (dataTab$p.adjust[i] < 0.05 & dataTab$p.adjust[i] >= 0.01  ) {
      dataTab$Pval_Adj[i]<- " < 0.05"
    } else if (dataTab$p.adjust[i] < 0.01 & dataTab$p.adjust[i] >= 0.001) {
      dataTab$Pval_Adj[i]<- " < 0.01"
    } else if (dataTab$p.adjust[i] < 0.001) {
      dataTab$Pval_Adj[i]<- " < 0.001"
    }
  }
  return(dataTab)
}

keggComplete$p.adjust[is.na(keggComplete$p.adjust)] <-  1

keggComplete<-addPCat(keggComplete)

keggComplete$Cluster = factor(keggComplete$Cluster , levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cluster, size = Count, color = Pval_Adj)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 12), breaks = c(0,5,10,20))+ 
  scale_color_viridis(discrete=T ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  xlab("KEGG Pathways") +
  theme(text = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =14, color = "black")) +
  theme(axis.text.y = element_text(size =18, color = "black")) &
  guides(color = guide_legend(override.aes = list(size = 8)))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top10Genes_KeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()



# AUcell

deDat<-read.csv('AUCell_KeggClustProF_All_ContrVsStress_2022-11-28.csv')

topMark = getTopMarkers(df = deDat, topNumb = 10)
colnames(deDat)[7] = "Description"
#topMarkDf = deDat[(deDat$Description%in%topMark),]

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Description)
  allDesc<-character()
  allClust<-character()
  allP<-numeric()
  allCount<-numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cell_Type==cluster),]
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
  colnames(finalTable) = c("Description", "p_val", "avg_log2FC", "Cell_Type")
  return(finalTable)
}

missDat<-findMissing(dataTab =  deDat, clusters = clusters)

keggComplete<-rbind.fill(deDat, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]

keggComplete$PvalMod = keggComplete$p_val_adj

keggComplete$PvalMod[keggComplete$PvalMod == 0] <-  1e-300

keggComplete$`-Log10(Pval_Adj)` = log10(keggComplete$PvalMod) * -1

keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cell_Type, size = `-Log10(Pval_Adj)`, color = avg_log2FC)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 12))+ 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  xlab("KEGG Pathways") + 
  theme(text = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =18, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top10Genes_AUCellKeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()