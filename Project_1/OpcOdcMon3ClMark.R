library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

curDate<-Sys.Date()

set.seed(23)

targDir <- 'OPC_ODC/Monocle3/86PC/'
dir.create(targDir)

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

#Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)

cds <- as.cell_data_set(RNA.combined.norm)


cds <- preprocess_cds(cds, num_dim = 86)

cds <- reduce_dimension(cds)

#cds <- preprocess_cds(cds, num_dim = 100)

#plot_pc_variance_explained(cds)


#cds <- reduce_dimension(cds)

#plot_cells(cds)

#plot_cells(cds, color_cells_by="integrated_snn_res.0.2")

#colnames(cds@colData)

metadat<-data.frame(cds@colData)

cds <- align_cds(cds, num_dim = 86, alignment_group = "group")
cds <- reduce_dimension(cds)
conStr<-plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE, cell_size=1)

#ggsave(paste0(targDir,'OPC_ODC_MonocClust_86PC_Root2_StressVsContr', curDate, '.jpeg'), plot = conStr, height = 12, width = 12, units = 'in', dpi = 300)

#cds <- cluster_cells(cds, resolution=2e-3)
cds <- cluster_cells(cds)
#plot_cells(cds)

plot_cells(cds, color_cells_by = "cluster")

q1<-cds@clusters$UMAP

clustID<-q1$clusters

cellNames<-names(clustID)
cellClusters<-as.character(clustID)
cellClustInf<-data.frame(cbind(cellNames, cellClusters))

### pseudotime
# learn trajectory
cds <- learn_graph(cds, use_partition = T)

regPlot<-monocle3::plot_cells(cds, group_label_size = 5, graph_label_size = 4)

regPlot

cds<-order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) =='2']))

psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                     color_cells_by = 'pseudotime',
                     label_branch_points = F,
                     label_roots = F,
                     label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

psedPlot

plotClust<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                      color_cells_by = 'cluster',
                      label_branch_points = F,
                      label_roots = F,
                      label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

plotClust

ggsave(paste0(targDir, 'OPC_ODC_MonocClust_86PC_Root2_ColClust_3', curDate, '.jpeg'), plot = plotClust, height = 12, width = 12, units = 'in', dpi = 300)

cds$monocle3_pseudotime<- pseudotime(cds)

data.pseudo<- data.frame(colData(cds))

identical(row.names(data.pseudo), cellClustInf$cellNames)

data.pseudo$monocClust<-cellClustInf$cellClusters

psedBox<-ggplot(data.pseudo, aes(monocle3_pseudotime, cellClusters, fill= cellClusters))+
  geom_boxplot()+
  theme(text = element_text(size = 20)) 

psedBox

medTimes<-data.pseudo%>%
  group_by(monocClust)%>% 
  summarise(Mean=mean(monocle3_pseudotime), Median=median(monocle3_pseudotime))

ggsave(paste0(targDir, 'OPC_ODC_MonocClust_86PC_Root2_Clust3', curDate, '.jpeg'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_86PC_Root2_Box_Clust3', curDate, '.jpeg'), plot = psedBox, height = 12, width = 12, units = 'in', dpi = 300)
write.csv(medTimes, paste0(targDir, 'OPC_ODC_MonocClust_86PC_Root2_sum_', curDate, '.csv'), row.names = F)

monocUmap<-plot_cells(cds, label_branch_points = F,label_roots = F, cell_size=1, group_label_size = 4)+
  theme(text = element_text(size = 20)) 
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_UMAP_', curDate, '.jpeg'), plot = monocUmap, height = 12, width = 12, units = 'in', dpi = 300)

oldLabel<-plot_cells(cds, color_cells_by="integrated_snn_res.0.2", label_roots = F, label_leaves = F, group_label_size = 4, cell_size=1)+
  theme(text = element_text(size = 20)) 
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_UMAP_OldLabels', curDate, '.jpeg'), plot = oldLabel, height = 12, width = 12, units = 'in', dpi = 300)

###
traj.plot <- plot_cells(cds)
point.data <- ggplot_build(traj.plot)[["plot"]][["data"]]

# find markers
marker_test_res <- top_markers(
  cds,
  group_cells_by = "cluster",
  genes_to_test_per_group = 100,
  reduction_method = "UMAP",
  marker_sig_test = TRUE,
  speedglm.maxiter = 25,
  cores = 4
)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

marker_test_group <- top_markers(
  cds,
  group_cells_by = 'group',
  genes_to_test_per_group = 100,
  reduction_method = "UMAP",
  marker_sig_test = TRUE,
  speedglm.maxiter = 25,
  cores = 4
)


write.csv(marker_test_res, paste0(targDir, 'OPC_ODC_MonocClust_86PC_clusterMarkersTop100_', curDate, '.csv'), row.names = F)
write.csv(marker_test_group, paste0(targDir,'OPC_ODC_MonocClust_86PC_groupMarkersTop100_', curDate, '.csv'), row.names = F)

# find correlations between gene expression and pseudo time
cds$monocle3_pseudotime<- pseudotime(cds)
metadat<-data.frame(cds@colData)
clustGenes<-unique(marker_test_res$gene_id)

gene_fits <- fit_models(cds, model_formula_str = "~monocle3_pseudotime", expression_family="negbinomial")
fit_coefs <- coefficient_table(gene_fits)
coefFiltr<-data.frame(fit_coefs[(fit_coefs$status == 'OK') & (fit_coefs$term == 'monocle3_pseudotime'),])
coefFiltr<-coefFiltr[, c(1:2, 5:ncol(coefFiltr))]

write.csv(coefFiltr, paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_', curDate, '.csv'), row.names = F)
# take group into account
genes_adjGroup <- fit_models(cds, model_formula_str = "~monocle3_pseudotime + group", expression_family="negbinomial")
fit_adjGroup <- coefficient_table(genes_adjGroup)
coefAdjGrFiltr<-fit_adjGroup[(fit_adjGroup$status == 'OK') & (fit_adjGroup$term == 'monocle3_pseudotime'),]
coefAdjGrFiltr<-coefAdjGrFiltr[, c(1:2, 5:ncol(coefAdjGrFiltr))]
write.csv(coefAdjGrFiltr, paste0(targDir, 'OPC_ODC_MonocClust_86PC_top100G_adjGr_TimeCor_', curDate, '.csv'), row.names = F)
# compare models
modelComp<-compare_models(genes_adjGroup, gene_fits)

write.csv(modelComp, paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_ModelComp_', curDate, '.csv'), row.names = F)

saveRDS(cds, 'OpcOdcInt_MonocClust_86PC')