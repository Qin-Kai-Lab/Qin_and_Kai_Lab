setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

df<-data.frame(RNA.combined.norm@reductions$umap@cell.embeddings)
df$Cell_id<-row.names(df)

metadata<-data.frame(RNA.combined.norm@meta.data)
metadata$Cell_id<-row.names(metadata)

metaSel<-metadata[, c("Cell_id", "integrated_snn_res.0.2", "group")]

combDat<-plyr::join(df, metaSel, by="Cell_id", type="left", match='all')

row.names(combDat)<-combDat$Cell_id

write.csv(combDat, 'rnaSeqNorm_loopBrowser.csv')

##

df<-read.csv('rnaSeqNorm_loopBrowser.csv')



df$Cell_id<-gsub("\\_.*", "", df$Cell_id)

df$X<-df$Cell_id

# control
control<-df[(df$group=='Control'),]

controlDat<-control[, 1:3]

colnames(controlDat)[1]<-'Barcode'

controlMeta<-control[, 4:5]

colnames(controlMeta)[1]<-'Barcode'

write.csv(controlDat, 'Opc_Odc_controlUmap.csv', row.names = F)

write.csv(controlMeta, 'Opc_Odc_controlMeta.csv', row.names = F)

#write.csv(control, 'Opc_Odc_controlLoop.csv', row.names = F)

# stress
stress<-df[(df$group=='Stress'),]

stressDat<-stress[, 1:3]

colnames(stressDat)[1]<-'Barcode'

stressMeta<-stress[, 4:5]

colnames(stressMeta)[1]<-'Barcode'

write.csv(stressDat, 'Opc_Odc_stressUmap.csv', row.names = F)

write.csv(stressMeta, 'Opc_Odc_stressMeta.csv', row.names = F)
