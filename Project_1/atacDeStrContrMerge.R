library(Seurat)
library(Signac)
library(ggplot2)
library(gridExtra)
library(EnsDb.Mmusculus.v79)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

atacInt<-readRDS('atacMerged_macs2')

DefaultAssay(atacInt)

#q1<-atacInt@assays$Combined_peaks@counts

#q2<-atacInt@assays$Combined_peaks@data

#q3<-atacInt@assays$Combined_peaks@scale.data

#identical(q1,q2)

#identical(q3,q2)


clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

DimPlot(atacInt)

#
# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(atacInt) <- annotations

# all markers
#allMarkers <- FindAllMarkers(atacInt , only.pos = T, min.pct = 0.05, test.use = "LR", 
#latent.vars = 'atac_peak_region_fragments')
#closGenes<-ClosestFeature(atacInt, regions = rownames(allMarkers))

#allMarkEdit<-cbind(allMarkers, closGenes)

saveRDS(atacInt, file = 'atacMerged_macs2')


targDir<-'./atacMergeDE/'

dir.create(targDir)

#write.csv(allMarkEdit, file = paste0(targDir, 'clustSpecificMarkers.csv'), row.names = F)

#for ( i in clusters ){
#  dfSel<-allMarkEdit[(allMarkEdit$cluster==i),]
# write.csv(dfSel, file = paste0(targDir, i, '_clustSpecificMarkers.csv'), row.names = F)
#}

# function to claculate differential gene expression
findMarkersGr<-function(dat, clust, pos, latVar, test){
  combMarkers<-data.frame(matrix(ncol=0, nrow=0))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$dataset
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.05, test.use = test, 
                              latent.vars = latVar, 
                              ident.1 = "Stress", ident.2 = "Control")
    closGenes<-ClosestFeature(atacInt, regions = rownames(allMarkers))
    
    allMarkEdit<-cbind(allMarkers, closGenes)
    allMarkEdit<-tibble::add_column(allMarkEdit, Cell_Type=i, .before = 1)
    write.csv(allMarkEdit, paste0(targDir, 'atac_DE_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    combMarkers<-rbind(combMarkers, allMarkEdit)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=atacInt, clust=clusters, pos=F, latVar = 'atac_peak_region_fragments', test='LR') # atac_peak_region_fragments

write.csv(groupMarkers, paste0('allDe_atacMerge_StrVsContr_',curDate, '.csv'), row.names = F)