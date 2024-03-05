
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

kegg<-download_KEGG(species= "mmu", keggType = "KEGG", keyType = "kegg")



entr<-as.data.frame(kegg$KEGGPATHID2EXTID)

keggPath<-as.data.frame(kegg$KEGGPATHID2NAME)

symbols <- mapIds(org.Mm.eg.db, keys = entr$to, keytype = "ENTREZID", column="SYMBOL")

entrezNames<-as.data.frame(names(symbols))

keggGenes<-data.frame(symbols)

keggGen<-cbind(keggGenes, entrezNames)

# rename and combine

colnames(entr)<-c('kegg_id', 'entrez')

colnames(keggGen)<-c('Gene', 'entrez')

allGenes<-plyr::join(keggGen, entr, by='entrez', type='left', match='all')

# add pathway names

colnames(keggPath)<-c('kegg_id', 'Kegg_path')

allKegg<-plyr::join(allGenes, keggPath, by='kegg_id', type='left', match='all')

new_DF <- allKegg[rowSums(is.na(allKegg)) > 0,]

write.csv(allKegg, 'keggMouse_clustprof_db.csv', row.names = F)