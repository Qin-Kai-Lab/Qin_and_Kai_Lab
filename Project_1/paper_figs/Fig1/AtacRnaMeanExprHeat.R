library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)
library(Seurat)
library(Signac)

corDf = read.csv("atacRna/PredictGeneActCor/Default/PredictGene_RNA_spearman_AllClusters_2023-10-23.csv")
clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'OPC', 'ODC', 'MG')
corDf = corDf[corDf$Cluster%in%clusters,]

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

DefaultAssay(atacFilt) = "PredictActivity"
atacFilt = ScaleData(atacFilt)


# DefaultAssay(atacFilt) = "RNA"
# atacFilt = ScaleData(atacFilt)

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

# get average RNA expression
 getMeanRNA = function(combGenes, cluster, atacFilt, slot) {
   objSub = subset(atacFilt, Annotations == cluster)
   curExpression = data.frame(AverageExpression(object=objSub,assays = "RNA", features = combGenes, group.by = "ident", slot = slot))
   colnames(curExpression)[1] = "Mean_RNA_Expression"
   curExpression$Gene = rownames(curExpression)
   rownames(curExpression) = NULL
   curExpression$Cluster = cluster
   return(curExpression)
 }
 
 combRNAExpr = data.frame()
 for (cluster in clusters ) {
   curExpr =  getMeanRNA(combGenes=combGenes, cluster=cluster, atacFilt=atacFilt, slot = "scale.data")
   combRNAExpr = rbind(combRNAExpr, curExpr)
 }
 
 # get average Predicted expression
 getMeanFrags = function(combGenes, cluster, atacFilt, slot) {
   objSub = subset(atacFilt, Annotations == cluster)
   curExpression = data.frame(AverageExpression(object=objSub,assays = "PredictActivity", features = combGenes, group.by = "ident", slot = slot))
   colnames(curExpression)[1] = "Mean_Fragments"
   curExpression$Gene = rownames(curExpression)
   rownames(curExpression) = NULL
   curExpression$Cluster = cluster
   return(curExpression)
 }
 
 combFrags = data.frame()
 for (cluster in clusters ) {
   curFras =   getMeanFrags(combGenes=combGenes, cluster=cluster, atacFilt=atacFilt, slot="scale.data")
   combFrags = rbind(combFrags, curFras)
 }
 
 # make heatmaps
 targDir = "atacRna/PredictGeneActCor/Heatmaps/MeanExpressions/"
 dir.create(targDir, recursive = T, showWarnings = F)
 
 topN = "10"
 
 
 curHeatMap = ggplot(combRNAExpr, aes(y=Gene, x=Cluster, fill=Mean_RNA_Expression))+
   geom_tile() +
   theme_classic() +
   scale_fill_viridis(discrete=F ) +
   theme(text = element_text(size = 24))+
   theme(legend.text = element_text(size = 18)) +
   theme(legend.title = element_text(size = 20)) +
   theme(axis.text.y = element_text(size = 14))
 
 print(curHeatMap)
 
 png(filename = paste0(targDir,"RNA_Expression_Scaled_Top_", topN, "_SpearCor_Default_2023-10-17.png"), width = 20, height = 16, units = "in", res = 300)
 print(curHeatMap)
 dev.off()
 
 curHeatMap = ggplot(combFrags, aes(y=Gene, x=Cluster, fill=Mean_Fragments))+
   geom_tile() +
   theme_classic() +
   scale_fill_viridis(discrete=F ) +
   theme(text = element_text(size = 24))+
   theme(legend.text = element_text(size = 18)) +
   theme(legend.title = element_text(size = 20)) +
   theme(axis.text.y = element_text(size = 14))
 
 print(curHeatMap)
 
 png(filename = paste0(targDir,"Fragments_Scaled_Top_", topN, "_SpearCor_Default_2023-10-17.png"), width = 20, height = 16, units = "in", res = 300)
 print(curHeatMap)
 dev.off()
 
 