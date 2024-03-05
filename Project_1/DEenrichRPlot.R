library(enrichR)

library(MAST)
library(limma)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")


source('../../programs/renameClusters.R')

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

levels(RNA.combined.norm)

DefaultAssay(RNA.combined.norm)


df<-read.csv('./DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv')

databases<-read.table('../data/enrichR_databases.txt', F)
databases<-as.character(databases$V1)

# chack that all target databases are present
length(databases[(databases%in%dbs$libraryName)])

for (i in databases){
  dirPath<-paste0('Enrichment/', i)
  dir.create((dirPath))
}

rm(dirPath, i)

# DE plots
plotDE<-function(dat, clust, pos){
  for ( i in clust){
    for ( db in databases ) {
      # get reference information
      dfRef<-df[(df$Cell_Type==i),]
      dfRefFilt<-dfRef[(dfRef$p_val_adj < 0.05),]
      if (nrow(dfRefFilt > 0)) {
        genesNumber<-nrow(dfRefFilt)
        pCut<-max(dfRefFilt$p_val)
        # subset cell type from all data
        cellType<-subset(x = dat, subset = Annotations == i)
        # change identity from cell type to group
        Idents(cellType)<-cellType$group
        # calculate gene expression
        plot1<-DEenrichRPlot(object = cellType, ident.1 = 'Stress', ident.2 = 'Control', logfc.threshold = 0.25, test.use = "MAST", enrich.database = db,
                            return.gene.list = FALSE, min.pct = 0.1, p.val.cutoff = pCut, num.pathway = 50, max.genes = genesNumber)
        
        dirPath<-paste0('Enrichment/', db, '/')
      
        ggsave(paste0(dirPath,  i, "_DEenrichRPlot_", curDate, '.jpeg'), plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
      }
    }
  }
}

# De Tables

groupMarkers<-plotDE(dat=RNA.combined.norm, clust=clusters, pos=F)











# get reference databases
dbs <- listEnrichrDbs()

write.csv(dbs, 'EnrichR_databases.csv', row.names = F)

# C-R
plotCR<-function(){
  df<-read.csv('C-R_DE_StressVsControl_wilcox_fdr.csv')
  dfRef<-df[(df$Cell_Type== 'C-R'),]
  dfRefFilt<-dfRef[(dfRef$fdr < 0.05),]
  genesNumber<-nrow(dfRefFilt)
  pCut<-max(dfRefFilt$p_val)
  cellType<-subset(x = RNA.combined.norm, subset = Annotations == 'C-R')
  Idents(cellType)<-cellType$group
  plot1<-DEenrichRPlot(object = cellType, ident.1 = 'Stress', ident.2 = 'Control', logfc.threshold = 0.25, test.use = "wilcox", enrich.database = "KEGG_2019_Mouse",
                       return.gene.list = FALSE, min.pct = 0.1, p.val.cutoff = pCut, num.pathway = 50, max.genes = genesNumber)
  
  ggsave(paste0('C-R_wilcox_fdr', "_DEenrichRPlot_kegg_2019_mouse_", curDate, '.jpeg'), plot = plot1, height = 16, width = 20, units = 'in', dpi = 300)
}

plotCR()
