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

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

# find top genes
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

topMark = getTopMarkers(df = groupMarkers, topNumb = 3)

RNA.combined.norm$Cluster_Group = paste(RNA.combined.norm$Annotations, RNA.combined.norm$group, sep = "_")

Group = unique(RNA.combined.norm$group)

Cluster_Group <- lapply(clusters, function(x) {
  lapply(Group, function(y) {
    paste(x, y, sep = "_")
  })
})

Cluster_Group = unlist(Cluster_Group)

RNA.combined.norm$Cluster_Group <- factor(RNA.combined.norm$Cluster_Group , levels = Cluster_Group)

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 14, group.by = 'Cluster_Group')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size

dotPlot

png(paste0(targDir,'DotPlot_Top3Genes_ClusterGroup_', curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
print(dotPlot)
dev.off()

# custom dot plot
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
topMark = getTopMarkers(df = groupMarkers, topNumb = 8)

colnames(groupMarkers)[7] = "Description"
#topMarkDf = groupMarkers[(groupMarkers$Description%in%topMark),]

clusters = unique(groupMarkers$Cell_Type)
missDat<-findMissing(dataTab =  groupMarkers, clusters = clusters)
keggComplete<-rbind.fill(groupMarkers, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]
keggComplete$p_val_adj[is.na(keggComplete$p_val_adj)] = 1

keggComplete$PvalMod = keggComplete$p_val
keggComplete$PvalMod[keggComplete$PvalMod == 0] <-  1e-300
keggComplete$`-Log10Pval` = log10(keggComplete$PvalMod) * -1
keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cell_Type, size = `-Log10Pval`, color = avg_log2FC)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 14))+ 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =18, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top8Genes_RNA_', curDate, '.png'), height = 16, width = 26, units = 'in', res = 300)
print(custDotPlot)
dev.off()

# kegg 

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

topMark = getTopPath(groupMarkers = groupMarkers, topNumb = 8)

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

keggComplete$Cluster = factor(keggComplete$Cluster, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cluster, size = Count, color = P.value)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 16), breaks = c(0,5,10,20))+ 
  scale_color_viridis(discrete=T ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =14, color = "black")) +
  theme(axis.text.y = element_text(size =14, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top8Genes_KeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()

# AUcell

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AUCmat)
  DefaultAssay(RNA.combined.norm) <- 'AUC'
  RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))
  return(RNA.combined.norm)
}

# add AUcell data
RNA.combined.norm<-addEnrich(x='./cellsAUC_keggClustProf.RData')

DefaultAssay(RNA.combined.norm)

# add differential expression data

deDat<-read.csv('AUCell_KeggClustProF_All_ContrVsStress_2022-11-28.csv')

topMark = getTopMarkers(df = deDat, topNumb = 5)

RNA.combined.norm$Cluster_Group = paste(RNA.combined.norm$Annotations, RNA.combined.norm$group, sep = "_")

Group = unique(RNA.combined.norm$group)

Cluster_Group <- lapply(clusters, function(x) {
  lapply(Group, function(y) {
    paste(x, y, sep = "_")
  })
})

Cluster_Group = unlist(Cluster_Group)

RNA.combined.norm$Cluster_Group <- factor(RNA.combined.norm$Cluster_Group , levels = Cluster_Group)

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 14, group.by = 'Cluster_Group')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+ # all text elements size
  theme(axis.text = element_text(size = 14))  # axes text size

dotPlot

png(paste0(targDir,'DotPlot_Top5Genes_AUCellKeggClustProf_ClusterGroup_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(dotPlot)
dev.off()

# custom dot plot

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

keggComplete$PvalMod = keggComplete$p_val

keggComplete$PvalMod[keggComplete$PvalMod == 0] <-  1e-300

keggComplete$`-Log10Pval` = log10(keggComplete$PvalMod) * -1

keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cell_Type, size = `-Log10Pval`, color = avg_log2FC)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 16))+ 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =14, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top10Genes_AUCellKeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()