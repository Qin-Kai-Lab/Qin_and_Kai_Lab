library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)
library(Hmisc)
library(plyr)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig5/'

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

keggdb = read.csv("keggMouse_clustprof_db.csv")

unique(groupMarkers$Cell_Type)

keggMark = read.csv("OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

unique(keggMark$Cell_Type)

keggGenes = unique(keggdb$Gene)
length(keggGenes)

myGenes = unique(groupMarkers$Genes)
length(myGenes)

length(myGenes[myGenes%in%keggGenes])
length(keggGenes[keggGenes%in%myGenes])

allGenes = unique(row.names(RNA.combined.norm))

length(keggGenes[keggGenes%in%allGenes])

keggAbs = data.frame(keggGenes[!(keggGenes%in%allGenes)])

# manually adjust p values
adjPval = function(df, cluctCol, pCol) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(df[, cluctCol])
  for (cluster in clusters) {
    dfSub = df[(df[cluctCol] == cluster),]
    dfSub$fdr_p = p.adjust(dfSub[, pCol], method = 'fdr')
    dfSub$bonferroni_p = p.adjust(dfSub[, pCol], method = 'bonferroni')
    combDat = rbind(combDat, dfSub)
  }
  return(combDat)
}

markAdj = adjPval(df = groupMarkers, cluctCol = "Cell_Type", pCol = "p_val")
write.csv(markAdj, paste0(targDir, 'DiffExprMonoc3ClustAdjP_', curDate, ".csv"), row.names = F)


# adjust p values based on kegg pathways
grep("Gfus", allGenes, ignore.case = T)

allCurRes = unique(groupMarkers$Genes)

getKeggSignGenes = function(df, db, cluctCol, pCol, genesDf) {
  curFile <- file(paste0(targDir, 'checkKeggGenesCurDE_', curDate, ".txt"), "w")
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(df[, cluctCol])
  for (cluster in clusters) {
    dfSub = df[(df[cluctCol] == cluster) & (df[pCol] < 0.05) ,]
    keggList = unique(dfSub$Genes)
    dbSub = db[(db$Kegg_path %in% keggList),]
    genesList = unique(dbSub$Gene)
    genesDfSub = genesDf[(genesDf$Genes %in% genesList) & (genesDf$Cell_Type == cluster),]
    
    line1 = paste("Cluster", cluster, 'need ', length(genesList), sep = "  ")
    line2 = paste("Cluster", cluster, 'found in current results', length(unique(genesDfSub$Genes)), sep = "  ")
    line3 = paste("Cluster", cluster, 'found in current all Genes', length(unique(genesList[genesList%in% allGenes])), sep = "  ")
    line4 = paste("Cluster", cluster, 'found in current all Rsults', length(unique(genesList[genesList%in% allCurRes])), sep = "  ")
    # write in the file
    print(line1)
    print(line2)
    print(line3)
    print(line4)
    lines <- c(line1, line2, line3, line4)
    writeLines(lines, curFile)
    curData = data.frame(genesList)
    curData$Cluster = cluster
    combDat = rbind(combDat, genesDfSub)
  }
  close(curFile)
  return(combDat)
}

genesOfInt = getKeggSignGenes(df = keggMark, db = keggdb, cluctCol = "Cell_Type", pCol = "p_val", genesDf = groupMarkers)

# adjust pvalues
markAdj = adjPval(df = genesOfInt, cluctCol = "Cell_Type", pCol = "p_val")
write.csv(markAdj, paste0(targDir, 'DiffExprMonoc3ClustAdjP_SignKeggGenes_', curDate, ".csv"), row.names = F)


# add kegg pathways to markers
list.files(targDir, pattern = '.csv')
curDf = read.csv(paste0(targDir, "DiffExprMonoc3ClustAdjP_SignKeggGenes_2023-06-12.csv"))
colnames(keggdb)[1] = "Genes"
joinDat = unique(plyr::join(curDf, keggdb, by = "Genes", type = "left", match = "all"))
unique(is.na(combDat$Kegg_path))

filterPath = function(df ,refDf) {
  combDat = data.frame(matrix(nrow = 0 , ncol = 0))
  clusters = unique(df$Cell_Type)
  for (cluster in clusters) {
    refSub = refDf[(refDf$Cell_Type == cluster) & (refDf$p_val < 0.05),]
    refPath = unique(refSub$Genes)
    dfSub = df[(df$Cell_Type == cluster),]
    dfSub = dfSub[(dfSub$Kegg_path%in%refPath),]
    combDat = rbind(combDat, dfSub)
  }
  combUniq = unique(combDat)
  return(combUniq)
}

combUnique = filterPath(df = joinDat, refDf = keggMark)

write.csv(combUnique, paste0(targDir, "DiffExprMonoc3ClustAdjP_SignKeggGenes_2023-06-20.csv"), row.names = F)
