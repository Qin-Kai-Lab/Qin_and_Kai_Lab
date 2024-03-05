library(monocle3)
library(Seurat)
library(SeuratWrappers)

curDate<-Sys.Date()

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

source('../../programs/renameClusters.R')

cds <- as.cell_data_set(RNA.combined.norm)

colData(cds)

fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partition<- c(rep(1, length(cds@colData@rownames)))

names(recreate.partition) <- cds@colData@rownames

recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

list_cluster<- RNA.combined.norm@active.ident

cds@clusters$UMAP$clusters<-list_cluster

cds@int_colData@listData$reducedDims$UMAP<-RNA.combined.norm$umap@cell.embeddings

monocle3::plot_cells(cds)

# learn trajectory
cds <- learn_graph(cds)

regPlot<-monocle3::plot_cells(cds, group_label_size = 5, graph_label_size = 4)

ggsave(paste0('combinedMonoc_NotOrd_', curDate, '.jpeg'), plot = regPlot, height = 12, width = 12, units = 'in', dpi = 300)

cds<-order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) =='CA1']))

psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
           color_cells_by = 'pseudotime',
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)

ggsave(paste0('combinedMonoc_CA1_Root_', curDate, '.jpeg'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)

cds$monocle3_pseudotime<- pseudotime(cds)

data.pseudo<- data.frame(colData(cds))

psedBox<-ggplot(data.pseudo, aes(monocle3_pseudotime, Annotations, fill= Annotations))+
  geom_boxplot()



ggsave(paste0('combinedMonocBox_CA1_Root_', curDate, '.jpeg'), plot = psedBox, height = 12, width = 12, units = 'in', dpi = 300)