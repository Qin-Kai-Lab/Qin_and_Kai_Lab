

source('../programs/renameClusters.R')

targDir = './Paper_figs/Fig1/FeatuePlot/Top/'

dir.create(targDir, recursive = T, showWarnings = F)

curDate = Sys.Date()

# plot function
exprPlot<-function(x, genesSet, curProj) {
  for (i in 1:length(x)) {
    genes<-x[[i]]
    fPlot<-FeaturePlot(curProj, features = genes, min.cutoff = "q9")
    
    ggsave(paste0(targDir,'FeaturePlot_', i, genesSet, curDate, '.jpeg'), plot = fPlot, height = 10, width = 12, units = 'in', dpi = 300)
  }
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

top10MarkersDf<-extractTopPosMarkers(x = allMarkersPct, topNumb = 5)

topUnique = unique(top10MarkersDf$gene)

topPosGenes = split(topUnique, ceiling(seq_along(topUnique)/6))

exprPlot(x = topPosGenes, genesSet = "_Top5PosGenes_", curProj = RNA.combined.norm)

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
allMarkersPct = allMarkers[(allMarkers$pct.1 >0.25) | (allMarkers$pct.2 >0.25), ]
topMarkersDf<-extractTopPosMarkers(x = allMarkersPct, topNumb = 5)
topUnique = unique(topMarkersDf$gene)
topPosGenes = split(topUnique, ceiling(seq_along(topUnique)/6))

exprPlot(x = topPosGenes, genesSet = "_AuCell_KeggClustProf_Top5PosGenes_", curProj = RNA.combined.norm)