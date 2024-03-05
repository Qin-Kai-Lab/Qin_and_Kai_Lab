dfDef = read.csv("Paper_figs/TasksList/P26/PredictGene_RNA_spearman_AllClusters_SepGroups_2024-01-03.csv")
df20 = read.csv("Paper_figs/TasksList/P26/PredictGene_20KB_RNA_spearman_AllClusters_SepGroups_2024-01-03.csv")

dfDef$ClustMeth = paste0(dfDef$Cluster, "_2KB")
dfDef$Method = "2KB"
df20$Method = "20KB"

df20$ClustMeth = paste0(df20$Cluster, "_20KB")

combDat = rbind(dfDef, df20)

minVal = min(combDat$Pval[combDat$Pval != 0]) * 0.5
combDat$Pval[combDat$Pval == 0] = minVal

combDat$`-Log10P` = log10(combDat$Pval) * -1

sumDat = data.frame(table(combDat$Group, combDat$Method))

sigDat = combDat[combDat$Pval < 0.05,]

sumSig = data.frame(table(sigDat$Group, sigDat$Method))

sumSig$Percent = sumSig$Freq / sumDat$Freq * 100

# mm = read.csv("~/Wang/output/Paper_figs/TasksList/P19/peakGeneCor_AllClusters_500K_Top10_MouseMotifs_2023-12-19.csv")
# mm = mm[mm$Cluster == "OPC",]
# 
# vert = read.csv("~/Wang/output/Paper_figs/TasksList/P19/peakGeneCor_AllClusters_500K_Top10_VertMotifs_2023-12-19.csv")
# vert = vert[vert$Cluster == "OPC",]