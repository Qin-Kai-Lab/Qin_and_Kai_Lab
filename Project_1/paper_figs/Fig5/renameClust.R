setwd("~/Wang/output")

targDir = './Paper_figs/Fig5/'


list.files("/home/flyhunter/Wang/output/OPC_ODC/StressVsContr/AUC_keggClustProf_DE", pattern = ".csv")
auDf = read.csv("/home/flyhunter/Wang/output/OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
rnaDf = read.csv("/home/flyhunter/Wang/output/OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")
aucDf = read.csv(paste0(targDir, "DiffExprMonoc3ClustAdjP_SignKeggGenes_2023-05-22.csv"))

# rename cluster
renameClusters = function(df, clustCol) {
  newClust = character()
  for ( i in 1:nrow(df)) {
    if (df[i, clustCol] == 1) {
      curClust = "ODC"
    } else if (df[i, clustCol] == 2) {
      curClust = "OPC"
    } else if ((df[i, clustCol] == 3) | (df[i, clustCol] == 4)) {
      curClust = "Intermideate"
    }
    newClust = c(newClust, curClust)
  }
  df$Cluster = newClust
  return(df)
}

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

rnaDf = renameClusters(df = rnaDf, clustCol = 'Cell_Type')
rnaDf = adjPval(df = rnaDf, cluctCol = 'Cell_Type', pCol = "p_val" )

write.csv(rnaDf, paste0(targDir, "DiffExprRNAMonoc3ClustAdjP_2023-06-08.csv"), row.names = F)

auDf = renameClusters(df = auDf, clustCol = 'Cell_Type')
auDf = adjPval(df = auDf, cluctCol = 'Cell_Type', pCol = "p_val" )

write.csv(auDf, paste0(targDir, "DiffExpr_AUC_keggClustProf_Monoc3ClustAdjP_2023-06-09.csv"), row.names = F)

aucDf = renameClusters(df = aucDf, clustCol = 'Cell_Type')
write.csv(aucDf, paste0(targDir, "DiffExprMonoc3ClustAdjP_SignKeggGenes_2023-06-12.csv"), row.names = F)