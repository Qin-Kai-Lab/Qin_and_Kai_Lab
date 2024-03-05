library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DimPlot(RNA.combined.norm, group.by = "MonocClust")

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3

DimPlot(RNA.combined.norm, group.by = "newMonocClust")

# markers

targetDir = "OPC_ODC/StressVsContr/RNA_DE/"
dir.create(targetDir, recursive = T)


clusters = unique(RNA.combined.norm$newMonocClust)

findMarkersGr<-function(dat, clust, pos, test, logtr){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = newMonocClust == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$group
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.1, logfc.threshold = logtr, test.use = test, ident.1 = "Stress", ident.2 = "Control")
    allMarkers$Genes<-rownames(allMarkers)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    write.csv(allMarkers, paste0(targetDir,'DiffExprMonoc3Clust_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F, test = "MAST", logtr = 0.25)
write.csv(groupMarkers, paste0(targetDir,'DiffExprMonoc3Clust_AllClust_ContrVsStress_', curDate, '.csv'), row.names = F)


# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$Cell_Type)
  markers = character()
  for (cluster in clusters) {
    dfSub = df[(df[['Cell_Type']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$AbsLog, decreasing = T),]
    topMark = head(dfOrd$Genes, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = 'newMonocClust')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size

dotPlot

ggsave(plot = dotPlot, file = paste0(targetDir, "dotPlot_combMonocClust_", curDate, ".jpeg"), height = 12, width = 24, units = 'in', dpi = 300)

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

targetDir = "OPC_ODC/StressVsContr/AUC_keggClustProf_DE/"
dir.create(targetDir, recursive = T)

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F, test = "wilcox", logtr = 0)
write.csv(groupMarkers, paste0(targetDir,'DiffExprMonoc3Clust_AllClust_ContrVsStress_', curDate, '.csv'), row.names = F)

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = 'newMonocClust')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+ # all text elements size
  theme(axis.text = element_text(size = 16)) + 
  theme(axis.text.x = element_text(size = 14)) +
  scale_x_discrete(limits=rev)

dotPlot

ggsave(plot = dotPlot, file = paste0(targetDir, "dotPlot_combMonocClust_", curDate, ".jpeg"), height = 12, width = 24, units = 'in', dpi = 300)
