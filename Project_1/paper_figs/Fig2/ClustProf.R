library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

targDir = './Paper_figs/Fig2/'
dir.create(targDir, recursive = T, showWarnings = F)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

dfFilt<-groupMarkers[(groupMarkers$p_val_adj < 0.05),]
entrez <- as.data.frame(mapIds(org.Mm.eg.db, keys = dfFilt$Genes, keytype = "SYMBOL", column="ENTREZID"))
colnames(entrez)<-'entrez'
dfEntr<-cbind(dfFilt, entrez)
colnames(dfEntr)[1] = "Cluster"

table(dfEntr$Cluster)

allEnrichedKegg<-function(dataTab, clusters, targDir, curDate){
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( cluster in clusters){
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    if (nrow(dfSel) > 0) {
      genes<-unique(dfSel$entrez)
      
      ego2 <- enrichKEGG(gene = genes, organism = 'mmu',
                         pvalueCutoff = 0.05)
      
      ego2Res<-ego2@result
      ego2Res$Cluster = cluster
      
      dir.create(targDir, showWarnings = FALSE, recursive = T)
      
      tablePath<-paste0(targDir, cluster, '_', curDate, '.csv')
      write.csv(ego2Res, tablePath, row.names = F)
      combDat = rbind(combDat, ego2Res)
      
      # if (nrow(ego2Res[(ego2Res$p.adjust < 0.05),]) > 0){
      #   plot1<-barplot(ego2, showCategory=20)
      #   plotPath<-paste0(targDir, cluster, '_20Terms_', curDate, '.jpeg')
      #   
      #   png(filename = plotPath, height = 16, width = 20, units = 'in', res = 300)
      #   print(plot1)
      #   dev.off()
        
      #}
    }
    
  }
  return(combDat)
}

# p_val_adj in the name means that only genes that are significant based on seurat default pval adjustment are used

dfEntr<-cbind(dfFilt, entrez)
colnames(dfEntr)[1] = "Cluster"
keggAll = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig2/ClustProf/Kegg_all_p_val_adj/', curDate = curDate)
write.csv(keggAll, paste0(targDir, "AllClust_AllGenes_2023-11-15_p_val_adj.csv"), row.names = F)

# only positive 
dfEntr<-cbind(dfFilt, entrez)
colnames(dfEntr)[1] = "Cluster"
dfEntr = dfEntr[(dfEntr$avg_log2FC > 0),]
keggPos = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig2/ClustProf/Kegg_positive_p_val_adj/', curDate = curDate)
write.csv(keggPos, paste0(targDir, "AllClust_PosGenes_2023-11-15_p_val_adj.csv"), row.names = F)

# only negative 
dfEntr<-cbind(dfFilt, entrez)
colnames(dfEntr)[1] = "Cluster"
dfEntr = dfEntr[(dfEntr$avg_log2FC < 0),]
keggNeg = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig2/ClustProf/Kegg_negative_p_val_adj/', curDate = curDate)
write.csv(keggNeg, paste0(targDir, "AllClust_NegGenes_2023-11-15_p_val_adj.csv"), row.names = F)

# make DotPlots
# custom dot plot

targDir = './Paper_figs/Fig2/DotPlot/'


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

topMark = getTopPath(groupMarkers = keggAll, topNumb = 8)

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

missDat<-findMissing(dataTab = keggAll, clusters = clusters)
keggComplete<-rbind.fill(keggAll, missDat)
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

keggComplete$Description = gsub(" - Mus musculus \\(house mouse\\)", "", keggComplete$Description)


# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cluster, size = Count, color = P.value)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 16), breaks = c(0,5,10,20))+ 
  scale_color_viridis(discrete=T ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =14, color = "black")) +
  theme(axis.text.y = element_text(size =14, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top8Genes_KeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()

# positive kegg genes
topMark = getTopPath(groupMarkers = keggPos, topNumb = 8)
missDat<-findMissing(dataTab = keggPos, clusters = clusters)
keggComplete<-rbind.fill(keggPos, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]
keggComplete<-addPCat(keggComplete)
keggComplete$Cluster = factor(keggComplete$Cluster, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
keggComplete$Description = gsub(" - Mus musculus \\(house mouse\\)", "", keggComplete$Description)
# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cluster, size = Count, color = P.value)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 16), breaks = c(0,5,10,20))+ 
  scale_color_viridis(discrete=T ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =14, color = "black")) +
  theme(axis.text.y = element_text(size =14, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top8PositiveGenes_KeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()

# negative kegg genes
topMark = getTopPath(groupMarkers = keggNeg, topNumb = 8)
missDat<-findMissing(dataTab = keggNeg, clusters = clusters)
keggComplete<-rbind.fill(keggNeg, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]
keggComplete<-addPCat(keggComplete)
keggComplete$Cluster = factor(keggComplete$Cluster, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
keggComplete$Description = gsub(" - Mus musculus \\(house mouse\\)", "", keggComplete$Description)
# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cluster, size = Count, color = P.value)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 14), breaks = c(0,5,10,20))+ 
  scale_color_viridis(discrete=T ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 16, color = "black")) +
  theme(axis.text.x = element_text(size =14, color = "black")) +
  theme(axis.text.y = element_text(size =14, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top8NegativeGenes_KeggClustProf_', curDate, '.png'), height = 18, width = 24, units = 'in', res = 300)
print(custDotPlot)
dev.off()