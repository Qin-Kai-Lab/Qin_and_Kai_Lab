library(MAST)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

curDate<-Sys.Date()

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

# make plots
DimPlot(RNA.combined.norm, reduction = "umap", label = F)
ggsave( paste0('./bothCondUmapRenameClust0.2Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)

control<-subset(x = RNA.combined.norm, subset = group == "Control")
DimPlot(control, reduction = "umap", label = F)
ggsave( paste0('./controlUmapRenameClust0.2Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)

stress<-subset(x = RNA.combined.norm, subset = group == "Stress")
DimPlot(stress, reduction = "umap", label = F)
ggsave( paste0('./stressUmapRenameClust0.2Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)

# dot plot
genes<-c('Rbfox3', 'Snap25', 'Syt1', 'Fn1', 'Rxfp1', 'Meis2', 'Mpped1', 'Satb2', 
         'Cntn6', 'Kcnh5', 'Vwc2l', 'Nectin3', 'Trps1', 'Glis3', 'Prox1', 'Slc4a4', 'Gad1', 'Gad2', 
         'Ndnf', 'Reln', 'Calcrl', 'Cspg4', 'Vcan', 'Plp1', 'Ptgds', 'Trf', 'Cx3cr1', 'Hexb', 'P2ry12')

dimNames<-RNA.combined.norm@assays$RNA@data@Dimnames
all_genes<-dimNames[[1]]

markersPresent<-genes[(genes%in%all_genes)]
length(genes)
length(markersPresent)


dotPlot<-DotPlot(object = RNA.combined.norm, features = markersPresent, scale.max = 100, dot.scale = 16)+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size
#scale_y_discrete(limits=rev)

ggsave(paste0('./dotPlotRenameClust0.2Res_', curDate, '.jpeg'),plot = dotPlot,  
       height = 12, width = 24, units = 'in', dpi = 300)

# heatmap

allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

write.csv(allMarkers, paste0('allMarkersRenameClust0.2Res_',curDate, '.csv'), row.names = F)

allMarkersSign<-allMarkers[(allMarkers$p_val_adj<0.05),]

allMarkSignSum<-data.frame(table(allMarkersSign$cluster))
colnames(allMarkSignSum)<-c('Cluster', 'Significant_markers_number')
write.csv(allMarkSignSum, 'sign_all_markers_sum_RenameClust0.2Res.csv', row.names = F)

topMarkers<-character()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

extractTopMarkers<-function(x){
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$AbsLog2FC, decreasing = T),]
    top10<-head(dfSort, 10)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}
top10MarkersDf<-extractTopMarkers(allMarkersSign)

topMarkers<-unique(top10MarkersDf$gene)

allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers) + NoLegend()

ggsave(paste0('./heatMap_top10Genes_RenameClust0.2Res_', curDate, '.jpeg'), plot = allHeat, height = 16, width = 18, units = 'in', dpi = 300)


allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers)

ggsave(paste0('./heatMapLeg_top10Genes_RenameClust0.2Res_', curDate, '.jpeg'), plot = allHeat, height = 20, width = 20, units = 'in', dpi = 300)

rm(allHeat, allMarkers, allMarkersSign, top10MarkersDf, topMarkers, dotPlot, dimNames)

# only positive markers heatmap
allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
allMarkersSign<-allMarkers[(allMarkers$p_val_adj<0.05),]

topMarkers<-character()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

extractTopMarkers<-function(x){
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$AbsLog2FC, decreasing = T),]
    top10<-head(dfSort, 10)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}
top10MarkersDf<-extractTopMarkers(allMarkersSign)

topMarkers<-unique(top10MarkersDf$gene)

allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers) + NoLegend()

ggsave(paste0('./heatMap_top10PositiveGenes_RenameClust0.2Res_', curDate, '.jpeg'), plot = allHeat, height = 16, width = 18, units = 'in', dpi = 300)


allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers)

ggsave(paste0('./heatMapLeg_top10PositiveGenes_RenameClust0.2Res_', curDate, '.jpeg'), plot = allHeat, height = 20, width = 20, units = 'in', dpi = 300)

rm(allHeat, allMarkers, allMarkersSign, top10MarkersDf, topMarkers, dotPlot, dimNames)

# table cells per group
metadata<-data.frame(RNA.combined.norm@meta.data)
cellsPerClust<-data.frame(table(metadata$Annotations, metadata$group))
colnames(cellsPerClust)<-c('Cluster', 'Group', 'Cells_number')
write.csv(cellsPerClust, paste0('./Cells_per_cluster_group_res0.2_', curDate, '.csv'), row.names = F)