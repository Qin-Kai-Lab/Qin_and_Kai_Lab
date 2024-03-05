source('../programs/renameClusters.R')

library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)

targDir = './Paper_figs/Fig2/'
dir.create(targDir, recursive = T, showWarnings = F)

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

curDate<-Sys.Date()

levels(RNA.combined.norm)

DefaultAssay(RNA.combined.norm)

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")
keggMark = read.csv('AUCell_KeggClustProF_All_ContrVsStress_2022-11-28.csv')
keggdb = read.csv("keggMouse_clustprof_db.csv")

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
write.csv(markAdj, paste0(targDir, 'allDiffExprLogfc0.25_ContrVsStress_AdjP_', curDate, ".csv"), row.names = F)

allGenes = unique(row.names(RNA.combined.norm))
allCurRes = unique(groupMarkers$Genes)

getKeggSignGenes = function(df, db, cluctCol, pCol, genesDf) {
  #curFile <- file(paste0(targDir, 'checkKeggGenesCurDE_', curDate, ".txt"), "w")
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
    #line3 = paste("Cluster", cluster, 'found in current all Genes', length(unique(genesList[genesList%in% allGenes])), sep = "  ")
    #line4 = paste("Cluster", cluster, 'found in current all Rsults', length(unique(genesList[genesList%in% allCurRes])), sep = "  ")
    # write in the file
    print(line1)
    print(line2)
    #print(line3)
    #print(line4)
    #lines <- c(line1, line2, line3, line4)
    #writeLines(lines, curFile)
    curData = data.frame(genesList)
    curData$Cluster = cluster
    combDat = rbind(combDat, genesDfSub)
  }
  #close(curFile)
  return(combDat)
}

genesOfInt = getKeggSignGenes(df = keggMark, db = keggdb, cluctCol = "Cell_Type", pCol = "p_val", genesDf = groupMarkers)

# adjust pvalues
markAdj = adjPval(df = genesOfInt, cluctCol = "Cell_Type", pCol = "p_val")
#write.csv(markAdj, paste0(targDir, 'DiffExpr_SignKeggGenes_', curDate, ".csv"), row.names = F)

# atac genesScore matrix
groupMarkers = read.csv("/home/flyhunter/Wang/output/archContrStressRna/atacRna/StressVsControl/StressVsControl_GeneScoremat_allClusters2023-03-20.csv")
keggMark = read.csv('AUCell_KeggClustProF_All_ContrVsStress_2022-11-28.csv')
keggdb = read.csv("keggMouse_clustprof_db.csv")

colnames(groupMarkers)[7] = "Cell_Type"
genesOfInt = getKeggSignGenes(df = keggMark, db = keggdb, cluctCol = "Cell_Type", pCol = "p_val", genesDf = groupMarkers)

# manually adjust p values
adjPval = function(df, cluctCol) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(df[, cluctCol])
  for (cluster in clusters) {
    dfSub = df[(df[cluctCol] == cluster),]
    dfSub$FDR_P = p.adjust(dfSub$Pval, method = 'fdr')
    dfSub$BONFERRONI_P = p.adjust(dfSub$Pval, method = 'bonferroni')
    combDat = rbind(combDat, dfSub)
  }
  return(combDat)
}

markAdj = adjPval(df = genesOfInt, cluctCol = "Cell_Type")

nrow(markAdj[(markAdj$BONFERRONI_P < 0.05),])
nrow(groupMarkers[(groupMarkers$bonferonni_p < 0.05),])

nrow(markAdj[(markAdj$FDR_P < 0.05),])
nrow(groupMarkers[(groupMarkers$fdr_p < 0.05),])

write.csv(markAdj, paste0(targDir, 'archRGenecoreMat_SignKeggGenes_', curDate, ".csv"), row.names = F)

# atac macs2 matrix
groupMarkers = read.csv("/home/flyhunter/Wang/output/archContrStressRna/atacRna/StressVsControl/Macs2PerCluster/Annotated/StressVsControl_Macs2_allClusters2023-03-28.csv")
colnames(groupMarkers)[11] = "Cell_Type"
colnames(groupMarkers)[13] = "Genes"

genesOfInt = getKeggSignGenes(df = keggMark, db = keggdb, cluctCol = "Cell_Type", pCol = "p_val", genesDf = groupMarkers)
markAdj = adjPval(df = genesOfInt, cluctCol = "Cell_Type")

nrow(markAdj[(markAdj$BONFERRONI_P < 0.05),])
nrow(groupMarkers[(groupMarkers$bonferonni_p < 0.05),])

nrow(markAdj[(markAdj$FDR_P < 0.05),])
nrow(groupMarkers[(groupMarkers$fdr_p < 0.05),])

write.csv(markAdj, paste0(targDir, 'archRMacs2Mat_SignKeggGenes_', curDate, ".csv"), row.names = F)