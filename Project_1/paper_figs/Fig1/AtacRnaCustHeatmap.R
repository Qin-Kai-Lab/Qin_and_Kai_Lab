library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)

corDf = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/PredictGene_RNA_spearman_AllClusters_2023-10-10.csv")
clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

findMissing<-function(dataTab, clusters){
  keggList<-unique(dataTab$Gene)
  allGenes<-character()
  allClust<-character()
  allP<-numeric()
  allCor<-numeric()
  for (cluster in clusters) {
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    for (i in keggList){
      if (i %nin% dfSel$Gene){
        curGene<-i
        curClust<-cluster
        curP<-1
        curCor<-0
        allGenes<-c(allGenes, curGene)
        allClust<-c(allClust, curClust)
        allP<-c(allP, curP)
        allCor<-c(allCor, curCor)
        # combine all vectors and make a dataframe from them
      }
    }
  }
  finalTable<-data.frame(Gene=allGenes, fdr_p=allP, Estimate=allCor, Cluster=allClust )
  return(finalTable)
}

missDat<-findMissing(dataTab =  corDf, clusters = clusters)

CorComplete<-rbind.fill(corDf, missDat)

# heatmap on all genes

CorComplete$Cluster = factor(CorComplete$Cluster, levels =  c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
# heatmap with Counts
curHeatMap= ggplot(CorComplete, aes(y=Gene, x=Cluster, fill=Estimate))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(discrete=F ) +
  theme(text = element_text(size = 24))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 1))

curHeatMap

targDir = "atacRna/PredictGeneActCor/Heatmaps/"
dir.create(targDir, recursive = T, showWarnings = F)

topN = "all"
png(filename = paste0(targDir,"Custom_", topN, "_PredGene_RNA_spearman_2023-10-11.png"), width = 20, height = 16, units = "in", res = 300)
print(curHeatMap)
dev.off()

## heatmap with top 50 genes overall genes
CorComplete = CorComplete[order(CorComplete$Estimate, decreasing = T),]
topGenes = unique(head(CorComplete$Gene, 57))
CorTop = CorComplete[CorComplete$Gene%in%topGenes,]

curHeatMap = ggplot(CorTop, aes(y=Gene, x=Cluster, fill=Estimate))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(discrete=F ) +
  theme(text = element_text(size = 24))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 18))

print(curHeatMap)

topN = "50"
png(filename = paste0(targDir,"Custom_", topN, "_PredGene_RNA_spearman_2023-10-11.png"), width = 20, height = 16, units = "in", res = 300)
print(curHeatMap)
dev.off()


# top 10 per cluster

getTop = function(corDf,  topN, cluster) {
  corDfSel = corDf[corDf$Cluster == cluster,]
  corDfSel = corDfSel[corDfSel$Estimate > 0 & corDfSel$fdr_p < 0.05,]
  corDfSel = corDfSel[order(corDfSel$Estimate, decreasing = T),]
  topGenes = head(corDfSel$Gene, topN)
  return(topGenes)
}

combGenes = character()
for (cluster in clusters ) {
  curTop = getTop(corDf=CorComplete,  topN=10, cluster=cluster)
  combGenes = c(combGenes, curTop)
}

combGenes = unique(combGenes)

CorTop = CorComplete[CorComplete$Gene%in%combGenes,]

curHeatMap = ggplot(CorTop, aes(y=Gene, x=Cluster, fill=Estimate))+
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(discrete=F ) +
  theme(text = element_text(size = 24))+
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 14))

print(curHeatMap)

topN = "10"
png(filename = paste0(targDir,"Custom_", topN, "_PerCluster_PredGene_RNA_spearman_2023-10-11.png"), width = 20, height = 16, units = "in", res = 300)
print(curHeatMap)
dev.off()