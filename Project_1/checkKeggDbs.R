


allGenes<-rownames(countsMat)

rm(AUCmat, cells_AUC, countsMatq1, RNA.combined.norm)

kegg2019<-read.delim('../data/kegg_mouse_2019.tsv', T, sep='\t')

kegg_path<-kegg2019

kegg_list<-character()

for ( i in 1:nrow(kegg_path)){
  genes<-strsplit(kegg_path$Genes[i], "\t")
  
  genesFromat<-str_to_title(genes[[1]])
  
  kegg_list<-c(kegg_list, genesFromat)
}

presKegg2019<-allGenes[allGenes%in%kegg_list]

keggClust<-read.csv('../data/keggMouse_clustprof_db.csv')

clustGenes<-keggClust$Gene

presKeggClust<-allGenes[allGenes%in%clustGenes]


load('./cellsAUC_keggClustProf.RData')

clustMat <- AUCell::getAUC(cells_AUC)

load('./cellsAUC_kegg2019Mouse.RData')

kegg2019Mat <- AUCell::getAUC(cells_AUC)

ncol(clustMat)

ncol(kegg2019Mat)

nrow(clustMat)

nrow(kegg2019Mat)

length(presKegg2019)
length(presKeggClust)