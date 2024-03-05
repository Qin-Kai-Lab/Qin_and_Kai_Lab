
targDir = './Paper_figs/Fig2/HeatMap/'

dir.create(targDir, recursive = T, showWarnings = F)

source('../programs/renameClusters.R')

allMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

# get RNA seq UMAP
source("~/Wang/programs/Utilities/CombGroups.R")
combGroup = combGroups(gr1 = groups, gr2 = clusters)

RNA.combined.norm$Cluster_Group = paste(RNA.combined.norm$Annotations, RNA.combined.norm$group, sep = "_")
Idents(RNA.combined.norm) = RNA.combined.norm$Cluster_Group

RNA.combined.norm$Cluster_Group = factor(RNA.combined.norm$Cluster_Group, levels = c(combGroup))

source("~/Wang/programs/Utilities/getTop.R")

topMark = getTopMarkers(df = allMarkers, topNumb = 8, clustCol = "Cell_Type")

curHeat<-DoHeatmap(RNA.combined.norm, features = topMark,  label = F)+
  theme(text = element_text(size = 12)) +
  #guides(color="none") +
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16))



curHeat 

png(filename = paste0(targDir, "HeatMap_Top8Genes_StrContr.png"), width = 20, height = 16, units = "in", res = 300)
print(curHeat)
dev.off()

# alternative 

curHeat<-DoHeatmap(RNA.combined.norm, features = topMark,  label = T, size = 2.8)+
  theme(text = element_text(size = 12)) +
  guides(color="none") +
  theme(legend.text = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16))


#curHeat 

png(filename = paste0(targDir, "HeatMap_Top8Genes_labels_StrContr.png"), width = 28, height = 20, units = "in", res = 300)
print(curHeat)
dev.off()
