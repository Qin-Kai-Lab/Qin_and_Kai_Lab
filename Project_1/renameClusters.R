library(Seurat)
library(ggplot2)
library(stringr)


setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integRNADoublFilt')
DefaultAssay(RNA.combined.norm) <- "RNA"
DefaultAssay(RNA.combined.norm)
Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

currentClustNames<-levels(RNA.combined.norm)
newClustNames<-currentClustNames

for ( i in 1:length(newClustNames)){
  if (newClustNames[i] == 1 | newClustNames[i] == 11){
    newClustNames[i]<-'CA1'
  } else if ( newClustNames[i] == 2){
    newClustNames[i]<-'CA3'
  } else if ( newClustNames[i] == 4 ) {
    newClustNames[i]<-'GABA'
  } else if ( newClustNames[i] == 6 ) {
    newClustNames[i]<-'OPC'
  } else if ( newClustNames[i] == 12 ) {
    newClustNames[i]<-'C-R'
  } else if ( newClustNames[i] == 8 ) {
    newClustNames[i]<-'SUB'
  } else if ( newClustNames[i] == 7 ) {
    newClustNames[i]<-'CA2'
  } else if ( newClustNames[i] == 0 | newClustNames[i] == 3 | newClustNames[i] == 10 ) {
    newClustNames[i]<-'DG'
  } else if ( newClustNames[i] == 5 ) {
    newClustNames[i]<-'ODC'
  } else if ( newClustNames[i] == 9 ) {
    newClustNames[i]<-'MG'
  }
}


names(newClustNames) <- levels(RNA.combined.norm)
RNA.combined.norm <- RenameIdents(RNA.combined.norm, newClustNames)
RNA.combined.norm[["Annotations"]] <- Idents(object = RNA.combined.norm)


levels(x =RNA.combined.norm) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

#q1<-data.frame(RNA.combined.norm@meta.data$integrated_snn_res.0.2)
#q2<-RNA.combined.norm@meta.data$Annotations
#q3<-data.frame(cbind(q1,q2))
