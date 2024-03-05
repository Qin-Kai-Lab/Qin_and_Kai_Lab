
targDir = './Paper_figs/Fig1/'

dir.create(targDir, recursive = T, showWarnings = F)

source('../programs/renameClusters.R')

curDate = Sys.Date()

genes<-c('Rbfox3', 'Snap25', 'Syt1', 'Fn1', 'Rxfp1', 'Meis2', 'Mpped1', 'Satb2', 
         'Cntn6', 'Kcnh5', 'Vwc2l', 'Nectin3', 'Trps1', 'Glis3', 'Prox1', 'Slc4a4', 'Gad1', 'Gad2', 
         'Ndnf', 'Reln', 'Calcrl', 'Cspg4', 'Vcan', 'Plp1', 'Ptgds', 'Trf', 'Cx3cr1', 'Hexb', 'P2ry12')


makePlot<-function(x, genesSet, proj) {
  curHeat<-DoHeatmap(proj, features = x)+
    theme(text = element_text(size = 12)) +
    guides(color="none") +
    theme(legend.text = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16)) +
    theme(legend.title = element_text(size = 16))
  
  ggsave(paste0(targDir,'Heatmap', genesSet, curDate, '.jpeg'), plot =  curHeat, height = 16, width = 20, units = 'in', dpi = 300)
}

makePlot(x = genes, genesSet = "_CustomMarkers_", proj = RNA.combined.norm)

# top positive markers

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

top10MarkersDf<-extractTopPosMarkers(x = allMarkersPct, topNumb = 3)

topPosGenes = unique(top10MarkersDf$gene)

makePlot(x = topPosGenes, genesSet = "_Top3PosMarkers_", proj = RNA.combined.norm)

# AU cell

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

allMarkers = read.csv("/home/flyhunter/Wang/output/Paper_figs/Fig1/AuCell_KeggClustProfAllMarkers.csv" )
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]

top10MarkersDf<-extractTopPosMarkers(x = allMarkersPct, topNumb = 8)
topPosGenes = unique(top10MarkersDf$gene)
makePlot(x = topPosGenes, genesSet = "_Top8PosMarkers_AuCellKeggClustProf_", proj = RNA.combined.norm)