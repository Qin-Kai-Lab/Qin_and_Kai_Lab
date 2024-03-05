library(AUCell)
library(limma)
library(Hmisc)
library(plyr)
library(ggplot2)
library(viridis)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

curDate<-Sys.Date()

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

# load enrichment results

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

targetDir<-'./Enrichment/AUcell/Kegg_ClustProfiler/pos_and_neg/'

dir.create(targetDir, recursive = T)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

# function to claculate differential gene expression
findMarkersGr<-function(dat, clust, pos){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$group
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox", ident.1 = "Stress", ident.2 = "Control")
    allMarkers$Genes<-rownames(allMarkers)
    write.csv(allMarkers, paste0(targetDir, i, '_StressVsContr_', curDate, '.csv'), row.names = F)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F)
write.csv(groupMarkers, paste0(targetDir, 'All_ContrVsStress_',curDate, '.csv'), row.names = F)

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
