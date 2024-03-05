library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

targDir = './Paper_figs/Fig5/'

markers = read.csv(paste0(targDir, "DiffExprRNAMonoc3ClustAdjP_2023-06-08.csv"))

clusters= unique(markers$Cluster)

for ( cluster in clusters ) {
  dfSub = markers[markers$Cluster == cluster,]
  dfSub = dfSub[dfSub$fdr_p < 0.05,]
  dfSub = dfSub[order(dfSub$fdr_p, decreasing = F),]
  write.csv(dfSub, paste0(targDir, "DiffExprRNA_", cluster, "_2023-10-16.csv"), row.names = F)
}

##
