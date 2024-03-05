library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)

curDate<-Sys.Date()

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

RNA.combined.norm <- ScaleData(RNA.combined.norm, verbose = FALSE)

# make combinations of groups
group1<-unique(Idents(RNA.combined.norm))

co <- t(as.data.frame(combn(unique(as.character(group1)),2)))

co<- data.frame(co)

clusters<-group1

# all Markers

allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

targDir<-'./OPC_ODC/Concerved/'

dir.create(targDir, recursive = T)

write.csv(allMarkers, paste0(targDir, 'All_clusters.csv'), row.names = F)

for (i  in clusters) {
  selCl<-allMarkers[(allMarkers$cluster==i),]
  write.csv(selCl, paste0(targDir, 'Cluster_', i, '.csv'), row.names=F)
}

# pairwise comparisons
targDir<-'./OPC_ODC/Pairwise_clusters/'

dir.create(targDir, recursive = T)


allComparisons<-data.frame(matrix(nrow=0, ncol=0))
for ( i in 1:nrow(co)){
  groups<-co[i,]
  ident1<-groups[1,1]
  ident2<-groups[1,2]
  clustMark<- FindMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", ident.1 = ident1, ident.2 = ident2)
  Comparison<-paste0(ident1, '_vs_', ident2)
  clustMark$Comparison<-Comparison
  write.csv(clustMark, paste0(targDir, 'Clusters_', Comparison, '.csv'), row.names = F)
  
  allComparisons<-rbind(allComparisons, clustMark)
  
}

write.csv(allComparisons, paste0(targDir, 'All_clusters.csv'), row.names = F)

