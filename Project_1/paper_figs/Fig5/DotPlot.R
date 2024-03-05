library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)
library(Hmisc)
library(plyr)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig5/DotPlot/'

dir.create(targDir, recursive = T, showWarnings = F)

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DimPlot(RNA.combined.norm, group.by = "MonocClust")

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"

Idents(RNA.combined.norm) = RNA.combined.norm$newMonocClust

levels(x =RNA.combined.norm) <- c('OPC', 'Intermideate', 'ODC')

DimPlot(RNA.combined.norm, group.by = "newMonocClust")

groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

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

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

RNA.combined.norm$newMonocClust <- factor(RNA.combined.norm$newMonocClust, levels = c("OPC", "Intermideate", "ODC"))

RNA.combined.norm$Group_Cluster = paste(RNA.combined.norm$group, RNA.combined.norm$newMonocClust, sep = "_")

RNA.combined.norm$Group_Cluster <- factor(RNA.combined.norm$Group_Cluster, levels = c("Control_OPC", "Stress_OPC", "Control_Intermideate", "Stress_Intermideate", "Control_ODC", "Stress_ODC"))
  
dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = 'Group_Cluster')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size

dotPlot

png(paste0(targDir,'DotPlot_Top5Genes_GroupCluster_', curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
dotPlot
dev.off()

# custom dot plot
groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
topMark = getTopMarkers(df = groupMarkers, topNumb = 10)
groupMarkers$Cell_Type = as.character(groupMarkers$Cell_Type)
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '1'] = "ODC"
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '2'] = "OPC"
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '3'] = "Intermideate"


colnames(groupMarkers)[7] = "Description"
#topMarkDf = groupMarkers[(groupMarkers$Description%in%topMark),]

clusters = unique(groupMarkers$Cell_Type)

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

missDat<-findMissing(dataTab =  groupMarkers, clusters = clusters)
keggComplete<-rbind.fill(groupMarkers, missDat)
keggComplete = keggComplete[(keggComplete$Description%in%topMark),]
keggComplete$p_val_adj[is.na(keggComplete$p_val_adj)] = 1

keggComplete$PvalMod = keggComplete$p_val
keggComplete$PvalMod[keggComplete$PvalMod == 0] <-  1e-300
keggComplete$`-Log10Pval` = log10(keggComplete$PvalMod) * -1
keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =   c('OPC', 'Intermideate', 'ODC'))

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

png(paste0(targDir,'Custom_DotPlot_Top10Genes_RNA_', curDate, '.png'), height = 12, width = 18, units = 'in', res = 300)
print(custDotPlot)
dev.off()


# AU cell

# load enrichment results
opcOdcCells = rownames(RNA.combined.norm@meta.data)

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  AucSub = AUCmat[, c(opcOdcCells)]
  RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AucSub)
  DefaultAssay(RNA.combined.norm) <- 'AUC'
  RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))
  return(RNA.combined.norm)
}

# add AUcell data
RNA.combined.norm<-addEnrich(x='./cellsAUC_keggClustProf.RData')
DefaultAssay(RNA.combined.norm)


groupMarkers = read.csv("OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"
RNA.combined.norm$newMonocClust <- factor(RNA.combined.norm$newMonocClust, levels = c("OPC", "Intermideate", "ODC"))

RNA.combined.norm$Group_Cluster = paste(RNA.combined.norm$group, RNA.combined.norm$newMonocClust, sep = "_")

RNA.combined.norm$Group_Cluster <- factor(RNA.combined.norm$Group_Cluster, levels = c("Control_OPC", "Stress_OPC", "Control_Intermideate", "Stress_Intermideate", "Control_ODC", "Stress_ODC"))

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = "Group_Cluster")+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+ # all text elements size
  theme(axis.text = element_text(size = 16)) + 
  theme(axis.text.x = element_text(size = 14)) +
  scale_x_discrete(limits=rev)

dotPlot

png(paste0(targDir,'DotPlot_Top5Genes_AuCell_KeggClustProf_GroupCluster', curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
dotPlot
dev.off()

# custom dot Plot

groupMarkers = read.csv("OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
topMark = getTopMarkers(df = groupMarkers, topNumb = 10)
groupMarkers$Cell_Type = as.character(groupMarkers$Cell_Type)
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '1'] = "ODC"
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '2'] = "OPC"
groupMarkers$Cell_Type[groupMarkers$Cell_Type == '3'] = "Intermideate"
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
keggComplete$Cell_Type = factor(keggComplete$Cell_Type, levels =   c('OPC', 'Intermideate', 'ODC'))

# plot
custDotPlot = ggplot(keggComplete, aes(x = Description, y = Cell_Type, size = `-Log10Pval`, color = avg_log2FC)) +
  geom_point()+
  theme_minimal() + 
  scale_size_continuous(range = c(1, 16))+ 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=rev) +
  theme(text = element_text(size = 20, color = "black")) +
  theme(axis.text.x = element_text(size =16, color = "black")) +
  theme(axis.text.y = element_text(size =18, color = "black"))

custDotPlot

png(paste0(targDir,'Custom_DotPlot_Top10Genes_AuCell_KeggClustProf_', curDate, '.png'), height = 16, width = 18, units = 'in', res = 300)
print(custDotPlot)
dev.off()

