# single-cell analysis package
library(Seurat)
library(harmony)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(irlba)
library(corrplot)
library(igraph)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 6)

curDate<-Sys.Date()

source('../programs/renameClusters.R')

RNA.combined.norm <- RNA.combined.norm %>% 
  RunHarmony("group", plot_convergence = TRUE)

seurat_obj = subset(RNA.combined.norm, subset = Annotations == "CA1")
#rm(RNA.combined.norm)
gc()


seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("group"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 20, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'group', # set the Idents of the metacell seurat object
  #min_cells = 89
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- FindVariableFeatures(object = seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj) )

seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj), fastpath=FALSE)
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='group')
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)

#p1 <- DimPlotMetacells(seurat_obj, group.by='Annotations') + umap_theme() + ggtitle("Annotations")
p2 <- DimPlotMetacells(seurat_obj, group.by='group') + umap_theme() + ggtitle("group")




seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Stress"), # the name of the group of interest in the group.by column
  group.by='group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)


seurat_obj <- ConstructNetwork(
  seurat_obj,
  setDatExpr=FALSE,
  tom_name = 'CA1_Stress', overwrite_tom = TRUE)

PlotDendrogram(seurat_obj, main='CA1_Stress')

TOM <- GetTOM(seurat_obj)

###
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="group"
)

hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

###
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'group', group_name = 'Stress'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "CA1_Stress"
)

p <- PlotKMEs(seurat_obj, ncol=5)

modules <- GetModules(seurat_obj)

hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)


plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

wrap_plots(plot_list, ncol=6)

ModuleCorrelogram(seurat_obj)

#
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)


p <- DotPlot(seurat_obj, features=mods, group.by = 'group')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

ModuleNetworkPlot(seurat_obj)
