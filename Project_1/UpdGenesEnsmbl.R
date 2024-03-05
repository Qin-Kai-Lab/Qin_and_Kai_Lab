control<-read.table('../data/control/extracted/features.tsv', F, sep= '\t')

stress<-read.table('../data/stress/extracted/features.tsv', F, sep= '\t')

contr<-control[, 1:2]

strs<-stress[, 1:2]

dictEns<-unique(rbind(contr, strs))

df<-read.csv('./DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv')

colnames(dictEns)<-c('Ens', 'Genes')

dfSelUpd<-plyr::join(df, dictEns, by='Genes', type='left', match='first')

dfNa<-dfSelUpd[is.na(dfSelUpd$Ens),]

'Pakap1'%in%dictEns$Genes

'Pakap'%in%dictEns$Genes


editDf<-dictEns[(dictEns$Genes=='Pakap'),]

for (i in 1:nrow(dfSelUpd)){
  if (dfSelUpd$Genes[i] == 'Pakap.1'){
    dfSelUpd$Ens[i]<-'ENSMUSG00000090053'
  }
}

dfNa<-dfSelUpd[is.na(dfSelUpd$Ens),]


write.csv(dfSelUpd, './DE_StressVsControl/mast_logfc0.25/allDiffExprLogfc0.25_ContrVsStress_2022-10-25_Ens.csv', row.names = F)