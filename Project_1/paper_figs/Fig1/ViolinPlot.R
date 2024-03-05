
targDir = './Paper_figs/Fig1/ViolinPlots/'

dir.create(targDir, recursive = T, showWarnings = F)

source('../programs/renameClusters.R')

curDate = Sys.Date()

# standard genes

# with standard genes too many plots to fit in one
genes<-c('Rbfox3', 'Snap25', 'Syt1', 'Fn1', 'Rxfp1', 'Meis2', 'Mpped1', 'Satb2', 
         'Cntn6', 'Kcnh5', 'Vwc2l', 'Nectin3', 'Trps1', 'Glis3', 'Prox1', 'Slc4a4', 'Gad1', 'Gad2', 
         'Ndnf', 'Reln', 'Calcrl', 'Cspg4', 'Vcan', 'Plp1', 'Ptgds', 'Trf', 'Cx3cr1', 'Hexb', 'P2ry12')
# split genes in groups of 8
markersList<-split(genes, ceiling(seq_along(genes)/6))

makePlots<-function(x, genesSet) {
  for (i in 1:length(x)) {
    genes<-x[[i]]
    vnPlot<-VlnPlot(
      object = RNA.combined.norm,
      features = genes,
      pt.size = 0.2
    )
    ggsave(paste0(targDir,'Violin_Plots_', i, genesSet, curDate, '.jpeg'), plot = vnPlot, height = 12, width = 18, units = 'in', dpi = 300)
  }
}

makePlots(x = markersList, genesSet = "_CustomMarkers_")

# plots with top markers

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv" )
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]


clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

extractTopPosMarkers<-function(x){
  topMarkers<-character()
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$avg_log2FC, decreasing = T),]
    top10<-head(dfSort, 3)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}

top10MarkersDf<-extractTopPosMarkers(allMarkersPct)

topUnique = unique(top10MarkersDf$gene)

topPosGenes = split(topUnique, ceiling(seq_along(topUnique)/6))

makePlots(x = topPosGenes, genesSet = "_Top3PosMarkers_")