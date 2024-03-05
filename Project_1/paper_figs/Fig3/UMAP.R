library(monocle3)
library(ggplot2)



cds<-readRDS('OpcOdcInt_MonocClust_86PC')

rnaDat = readRDS('integ_OPC_ODC')


psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                     color_cells_by = 'pseudotime',
                     label_branch_points = F,
                     label_roots = F,
                     label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

psedPlot

plotClust<-plot_cells(cds, graph_label_size = 4, group_label_size =0,
                      color_cells_by = 'cluster',
                      label_branch_points = F,
                      label_roots = F,
                      label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

plotClust

conStr<-plot_cells(cds, label_cell_groups=FALSE, graph_label_size = 4, 
                   color_cells_by = 'group',
                   label_branch_points = F,
                   label_roots = F,
                   label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20)) 

targDir = './Paper_figs/Fig3/'

ggsave(paste0(targDir, 'Umap_MonocClust_86PC_Clusters_2023-12-13.png'), plot = plotClust, height = 12, width = 12, units = 'in', dpi = 300)
ggsave(paste0(targDir, 'Umap_MonocClust_86PC_Group_2023-12-13.png'), plot = conStr, height = 12, width = 12, units = 'in', dpi = 300)
ggsave(paste0(targDir, 'Umap_MonocClust_86PC_Time_2023-12-13.png'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)


# get umap with seurat clusters


uMap = data.frame(cds@reduce_dim_aux@listData[["UMAP"]]@listData[["model"]]@listData[["umap_model"]][["embedding"]])

identical(colnames(rnaDat), colnames(cds))


##
source('../programs/renameClusters.R')

RNA.combined.norm$CurCells = colnames(RNA.combined.norm)

rnaFilt = subset(RNA.combined.norm, subset=CurCells %in% colnames(rnaDat))

identical(colnames(rnaFilt), colnames(rnaDat))
identical(colnames(rnaFilt), colnames(cds))
cds$Seurat_Clusters = rnaFilt$Annotations


surPlot = plot_cells(cds, graph_label_size = 4, label_cell_groups=FALSE,
           color_cells_by = 'Seurat_Clusters',
           label_branch_points = F,
           label_roots = F,
           label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

surPlot
targDir = './Paper_figs/Fig3/'

ggsave(paste0(targDir, 'Umap_MonocClust_86PC_Seurat_Clusters_2023-12-17.png'), plot = surPlot, height = 12, width = 12, units = 'in', dpi = 300)