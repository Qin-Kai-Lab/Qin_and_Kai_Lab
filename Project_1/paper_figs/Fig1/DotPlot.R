targDir = './Paper_figs/Fig1/DotPlot/'

dir.create(targDir, recursive = T, showWarnings = F)

source('../programs/renameClusters.R')

curDate = Sys.Date()

makePlot<-function(x, genesSet, proj) {
  dotPlot<-DotPlot(object = proj, features = x, scale.max = 100, dot.scale = 16)+
    scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(text = element_text(size = 24))+ # all text elements size
    theme(axis.text = element_text(size = 24))  # axes text size
  #scale_y_discrete(limits=rev)
  
  ggsave(paste0(targDir,'DotPlot', genesSet, curDate, '.jpeg'), plot =  dotPlot, height = 12, width = 24, units = 'in', dpi = 300)
}


# top positive genes
allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

extractTopPosMarkers<-function(x, topNumb){
  topMarkers<-character()
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$avg_log2FC, decreasing = T),]
    top10<-head(dfSort, topNumb)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}

top10MarkersDf<-extractTopPosMarkers(x= allMarkersPct, topNumb = 3)

topPosGenes = unique(top10MarkersDf$gene)

makePlot(x = topPosGenes, genesSet = "_Top3PosMarkers_", proj = RNA.combined.norm)

## AUCell

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

allMarkers <- FindAllMarkers(RNA.combined.norm , assay = "AUC", only.pos = F, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox")
write.csv(allMarkers, paste0(targDir, 'AuCell_KeggClustProfAllMarkers.csv'), row.names = F)

makePlot<-function(x, genesSet, proj) {
  dotPlot<-DotPlot(object = proj, features = x , scale.max = 100, dot.scale = 14)+
    scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
    theme(text = element_text(size = 16))+ # all text elements size
    theme(axis.text = element_text(size = 16)) + 
    theme(axis.text.x = element_text(size = 14)) +
    scale_x_discrete(limits=rev)
  
  ggsave(paste0(targDir,'DotPlot', genesSet, curDate, '.jpeg'), plot =  dotPlot, height = 12, width = 24, units = 'in', dpi = 300)
}

allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]
topMarkersDf<-extractTopPosMarkers(x = allMarkersPct, topNumb = 3)
topPosGenes = unique(topMarkersDf$gene)
makePlot(x = topPosGenes, genesSet = "_AuCell_KeggClustProf_Top3PosMarkers_", proj = RNA.combined.norm)