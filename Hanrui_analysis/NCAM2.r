.libPaths(c(
  LIBPATH,
  .libPaths()
))
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(harmony)
library(SeuratDisk)
library(WGCNA); allowWGCNAThreads()
library(readxl)
library(readr)
library(writexl)
library(tibble)
library(scCustomize)
library(tibble)
library(jsonlite)
library(dplyr)
library(tibble)
library(writexl)
library(stringr)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(cowplot)
library(harmony)

set.seed(123)

setwd(WD)
# test Rstudio plots
options(bitmapType = "cairo")
capabilities()[c("png", "cairo", "X11")]
plot(rnorm(50), rnorm(50)) 

# h5Seurat object from Green et al
astro <- LoadH5Seurat(
  "astrocytes.h5Seurat",
  assays = list(RNA = "counts"),
  reductions = FALSE,
  graphs     = FALSE,
  images     = FALSE,
  meta.data  = TRUE
)

astro
GetAssayData(astro,layer = "counts") 
DefaultAssay(astro)
xlsx_path <- "ROSMAP_ID_Status.xlsx"

samp <- read_excel("ROSMAP_ID_Status.xlsx",sheet = 1, na = c("NA","N/A",""))

samp$individualID <- trimws(as.character(samp$individualID))
astro@meta.data$individualID <- trimws(as.character(astro@meta.data$individualID))

samp$individualID_excel = samp$individualID
samp$individualID = NULL

names(samp)[!names(samp) %in% "individualID_excel"] <-
  paste0("clin_", names(samp)[!names(samp) %in% "individualID_excel"])

cell_barcodes <- rownames(astro@meta.data)

# Sample selection
astro@meta.data <- astro@meta.data %>%
  mutate(individualID_cell = individualID) %>%
  left_join(samp, by = c("individualID_cell" = "individualID_excel"))

rownames(astro@meta.data) <- cell_barcodes
head(astro@meta.data)
astro@meta.data <- astro@meta.data %>%
  mutate(
    AD_status = case_when(
      clin_braaksc %in% c(0, 1, 2) &
        clin_ceradsc == 4 &
        clin_cogdx == 1 &
        clin_dcfdx_lv == 1 ~ "Control",
      
      clin_braaksc %in% c(5, 6) &
        clin_ceradsc == 1 &
        clin_cogdx == 4 &
        clin_dcfdx_lv == 4 ~ "AD",
      
      TRUE ~ NA_character_  
    )
  )
astro_selected_samples <- subset(astro, subset = AD_status %in% c("AD","Control"))
head(astro_selected_samples@meta.data)
astro = NULL

oli <- LoadH5Seurat(
  "oligodendroglia.h5Seurat",
  assays     = list(RNA = "counts"),
  reductions = FALSE,
  graphs     = FALSE,
  images     = FALSE,
  meta.data  = TRUE
)

GetAssayData(oli,layer = "counts") 
DefaultAssay(oli)
xlsx_path <- "ROSMAP_ID_Status.xlsx"
samp <- read_excel("ROSMAP_ID_Status.xlsx",sheet = 1, na = c("NA","N/A",""))

samp$individualID <- trimws(as.character(samp$individualID))
oli@meta.data$individualID <- trimws(as.character(oli@meta.data$individualID))

samp$individualID_excel = samp$individualID
samp$individualID = NULL

names(samp)[!names(samp) %in% "individualID_excel"] <-
  paste0("clin_", names(samp)[!names(samp) %in% "individualID_excel"])

cell_barcodes <- rownames(oli@meta.data)

oli@meta.data <- oli@meta.data %>%
  mutate(individualID_cell = individualID) %>%
  left_join(samp, by = c("individualID_cell" = "individualID_excel"))

rownames(oli@meta.data) <- cell_barcodes

oli@meta.data <- oli@meta.data %>%
  mutate(
    AD_status = case_when(
      clin_braaksc %in% c(0, 1, 2) &
        clin_ceradsc == 4 &
        clin_cogdx == 1 &
        clin_dcfdx_lv == 1 ~ "Control",
      
      clin_braaksc %in% c(5, 6) &
        clin_ceradsc == 1 &
        clin_cogdx == 4 &
        clin_dcfdx_lv == 4 ~ "AD",
      TRUE ~ NA_character_  
    )
  )

oli_selected_samples <- subset(oli, subset = AD_status %in% c("AD","Control"))
oli = NULL
head(oli_selected_samples@meta.data)

# Then start processing pipeline
oli_ODC <- subset(oli_selected_samples, subset = subset == "Oligodendrocytes")
oli_selected_samples = NULL
oli_cells_sel = sample(Cells(oli_ODC), size = 1500)
objs_Oli_1500 <- subset(oli_ODC, cells = oli_cells_sel)
objs_merge <- merge(astro_selected_samples, y = objs_Oli_1500)
objs_merge <- Convert_Assay(objs_merge, assay = "RNA", convert_to = "Assay5")
class(objs_merge[["RNA"]])
objs_merge <- JoinLayers(objs_merge, assay = "RNA")
DefaultAssay(objs_merge) <- "RNA"


all_genes <- rownames(objs_merge)

genes_all <- all_genes
anno <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys=all_genes,
                        keytype="SYMBOL", columns=c("SYMBOL","CHR"))
)

anno2 <- anno %>%
  filter(!is.na(CHR)) %>%
  mutate(CHR = gsub("^chr","", CHR, ignore.case = TRUE)) %>%
  filter(CHR %in% as.character(1:22)) %>%
  distinct(SYMBOL, .keep_all = TRUE)
auto_genes <- intersect(all_genes, anno2$SYMBOL)
sex_or_nonauto <- setdiff(all_genes, auto_genes)
length(auto_genes); length(sex_or_nonauto); head(sex_or_nonauto)
objs_merge <- subset(objs_merge, features = auto_genes)

objs_merge <- NormalizeData(
  objs_merge,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = TRUE
  )

objs_merge <- FindVariableFeatures(
  objs_merge,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)


top100 <- head(VariableFeatures(objs_merge), 100)
top200 <- head(VariableFeatures(objs_merge), 200)
top500 <- head(VariableFeatures(objs_merge), 500)

objs_merge<- ScaleData(object = objs_merge,features = all_genes)
max_pcs <- 50
objs_merge <- RunPCA(
  object   = objs_merge,
  features = VariableFeatures(objs_merge),
  npcs     = max_pcs,
  verbose  = TRUE
)

objs_merge <- RunHarmony(
  object        = objs_merge,
  group.by.vars = "individualID",  
  dims.use          = 1:max_pcs,
  assay         = DefaultAssay(objs_merge),
  verbose       = TRUE
)

pcs_use <- 40

k_param <- 40
objs_merge <- FindNeighbors(
  object = objs_merge,
  reduction="harmony",
  dims = 1:pcs_use,
  k.param = k_param,
  verbose = TRUE
)


# res scan
res_grid <- seq(0.4, 2.8, by = 0.2)

objs_merge <- FindClusters(
  object = objs_merge,
  resolution = res_grid,
  verbose = TRUE
)
objs_merge@meta.data$DoubletFinder.score

res_cols <- grep(pattern = "_snn_res\\.", x = colnames(objs_merge@meta.data), value = TRUE)

cluster_summary <- lapply(res_cols, function(colnm) {
  Idents(objs_merge) <- objs_merge@meta.data[[colnm]]
  data.frame(
    resolution = sub(".*_snn_res\\.", "", colnm),
    n_clusters = length(levels(Idents(objs_merge))),
    min_size   = min(table(Idents(objs_merge))),
    median_size= stats::median(as.numeric(table(Idents(objs_merge)))),
    max_size   = max(table(Idents(objs_merge)))
  )
}) %>% bind_rows()
print(cluster_summary)

# decided to use 2.6
best_res <- 2.6

objs_merge <- FindClusters(
  object = objs_merge,
  resolution = best_res,
  verbose = TRUE
)

Idents(objs_merge) <- "seurat_clusters"
table(objs_merge@meta.data$individualID, objs_merge@meta.data$seurat_clusters)[,1:11]
table(objs_merge@meta.data$individualID, objs_merge@meta.data$seurat_clusters)[,12:25]
table(objs_merge@meta.data$subset, objs_merge@meta.data$seurat_clusters)

objs_merge <- Seurat::RunUMAP(
  object = objs_merge,
  reduction = "harmony",
  dims = 1:pcs_use
)

p_umap_subset <- Seurat::DimPlot(
  objs_merge,
  reduction = "umap",
  group.by = "subset",
  pt.size = 0.1
) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)

p_umap_subset
outdir = OUTDIR
ggsave(
      filename = file.path(outdir, paste0("1stCluster_subset_umap_ALL",".png")),
      plot = p_umap_subset,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )

p_umap_cluster <- Seurat::DimPlot(
  objs_merge,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)

p_umap_cluster

ggsave(
      filename = file.path(outdir, paste0("1stCluster_umap_ALL",".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )

for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_merge, subset = AD_status == grp)
    p_umap_subset <- Seurat::DimPlot(
  seu_grp,
  reduction = "umap",
  group.by = "subset",
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
  ggsave(
      filename = file.path(outdir, paste0("1stCluster_subset_umap_",grp,".png")),
      plot = p_umap_subset,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
  p_umap_cluster <- Seurat::DimPlot(
  seu_grp,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
  ggsave(
      filename = file.path(outdir, paste0("1stCluster_umap_",grp,".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
}


genes <- c("GFAP", "OSMR","SORL1","NCAM2")
outdir <- OUTDIR

for (gene in genes) {
  expr_all <- GetAssayData(
    object = objs_merge,
    assay = "RNA",
    slot = "data"
  )[gene, ]
  expr_all <- as.numeric(expr_all)
  gene_min <- min(expr_all, na.rm = TRUE)
  gene_max <- max(expr_all, na.rm = TRUE)
  message(gene, ": min = ", gene_min, ", max = ", gene_max)
  ## All cells
  p_all <- FeaturePlot(
    object = objs_merge,
    features = gene,
    reduction = "umap",
    order = TRUE,
    pt.size = 0.45,
    slot = "data",
    min.cutoff = gene_min,
    max.cutoff = gene_max
  ) +
    ggtitle(paste(gene, "All cells", sep = " - ")) +
    ppt_feature_theme+  
    ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
  print(p_all)
  ggsave(
    filename = file.path(outdir, paste0("052726_", gene, "_AllCells_FeaturePlot.png")),
    plot = p_all,
    width = 7,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  ## subset
  for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_merge, subset = AD_status == grp)
    p <- FeaturePlot(
      object = seu_grp,
      features = gene,
      reduction = "umap",
      order = TRUE,
      pt.size = 0.45,
      slot = "data",
      min.cutoff = gene_min,
      max.cutoff = gene_max
    ) +
      ggtitle(paste(gene, grp, sep = " - ")) +
      ppt_feature_theme + 
      ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
    print(p)
    ggsave(
      filename = file.path(outdir, paste0("052726_", gene, "_", grp, "_FeaturePlot.png")),
      plot = p,
      width = 7,
      height = 6,
      dpi = 600,
      bg = "white"
    )
  }
}


# Exclude clusters 16 and 20


Idents(objs_merge) <- "seurat_clusters"
objs_merge_Ast <- subset(objs_merge, idents = c("16","20"), invert = TRUE)

outdir = OUTDIR

p_umap_subset <- Seurat::DimPlot(
  objs_merge_Ast,
  reduction = "umap",
  group.by = "subset",
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
p_umap_subset
ggsave(
      filename = file.path(outdir, paste0("1stCluster_subset_umap_OliRemoved_ALL",".png")),
      plot = p_umap_subset,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
p_umap_cluster <- Seurat::DimPlot(
  objs_merge_Ast,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
p_umap_cluster
ggsave(
      filename = file.path(outdir, paste0("1stCluster_umap_OliRemoved_ALL",".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )



for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_merge_Ast, subset = AD_status == grp)
    p_umap_subset <- Seurat::DimPlot(
  seu_grp,
  reduction = "umap",
  group.by = "subset",
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
  ggsave(
      filename = file.path(outdir, paste0("1stCluster_subset_umap_OliRemoved_",grp,".png")),
      plot = p_umap_subset,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
  p_umap_cluster <- Seurat::DimPlot(
  seu_grp,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
  ggsave(
      filename = file.path(outdir, paste0("1stCluster_umap_OliRemoved_",grp,".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
}



genes <- GENELIST
outdir <- OUTDIR

for (gene in genes) {
  expr_all <- GetAssayData(
    object = objs_merge_Ast,
    assay = "RNA",
    slot = "data"
  )[gene, ]
  expr_all <- as.numeric(expr_all)
  gene_min <- min(expr_all, na.rm = TRUE)
  gene_max <- max(expr_all, na.rm = TRUE)
  message(gene, ": min = ", gene_min, ", max = ", gene_max)
  ## All cells
  p_all <- FeaturePlot(
    object = objs_merge,
    features = gene,
    reduction = "umap",
    order = TRUE,
    pt.size = 0.45,
    slot = "data",
    min.cutoff = gene_min,
    max.cutoff = gene_max
  ) +
    ggtitle(paste(gene, "All cells", sep = " - ")) +
    ppt_feature_theme+  
    ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
  print(p_all)
  ggsave(
    filename = file.path(outdir, paste0(gene, "_AllCells_FeaturePlot_OliRemoved.png")),
    plot = p_all,
    width = 7,
    height = 6,
    dpi = 600,
    bg = "white"
  )

  for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_merge, subset = AD_status == grp)
    p <- FeaturePlot(
      object = seu_grp,
      features = gene,
      reduction = "umap",
      order = TRUE,
      pt.size = 0.45,
      layer = "data",
      min.cutoff = gene_min,
      max.cutoff = gene_max
    ) +
      ggtitle(paste(gene, grp, sep = " - ")) +
      ppt_feature_theme + 
      ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
    print(p)
    ggsave(
      filename = file.path(outdir, paste0("052726_", gene, "_", grp, "_FeaturePlot_OliRemoved.png")),
      plot = p,
      width = 7,
      height = 6,
      dpi = 600,
      bg = "white"
    )
  }
}


Idents(objs_merge_Ast) <- "seurat_clusters"
DefaultAssay(objs_merge_Ast) <- "RNA"
avg_list <- AggregateExpression(
  objs_merge_Ast,
  assays = DefaultAssay(objs_merge_Ast),
  slot   = "data",
  group.by = "seurat_clusters",
  return.seurat = FALSE
)
# GDA > 0.8 genes
vvip_genes <- c(
  "APP","PSEN1","APOE","MAPT","PSEN2","GRN","ACE","TOMM40","NECTIN2","CLU",
  "PPARG","CD33","ABCA7","ACHE","BIN1","TREM2","GSK3B","PICALM","APOC1","BCHE",
  "TNF","BACE1","BDNF","ADAM10","IGF1","IL1B","SORL1","AGER","SNCA","TTR",
  "HFE","NOS3","A2M","MAOB","TFAM","VEGFA","CYP46A1","BACE2","ADAM17","KLC1",
  "NTF3","HLA-DRB1","MME","INS","PILRA","LEP","PLCG2","APH1B","CHRNA7","CRH",
  "ESR1","IDE","ABI3","ICAM1","TLR4","HTR6","IL10","MAOA","NFE2L2","CDK5",
  "CHAT","IL4","LRP1","IGFBP3","PLAU","EIF2AK2","S100B","UCHL1","IRS1","LRP8",
  "CSF1R","MPO","MS4A4A","CYP2D6","ENO1","NCSTN","SOD2","EIF2S1","PRNP","APOA1",
  "LDLR","ESR2","ECE1","ABCA1","IL2","CD2AP","CASS4","EPHA1","CST3","IGF2",
  "INPP5D","WWOX","INSR","RELN","BCL2","CASP3","HMOX1","IGF1R","NPY","BAX",
  "PTK2B","APBB1","NOS1","PON1","GAPDH","GCG","STAT3","PARP1","ADRA1A","AGTR1",
  "HMGCR","PPARA","CTSD","CYP19A1","RCAN1","GRIN2B","C3","HSPA1A","FYN","GLUL",
  "NTRK2","SERPINE2","KLK6","CLOCK","DLST","APOA4","IL1A","IL6R","LPL","MT3",
  "NGFR","NTRK1","VLDLR"
)
# GDA between 0.6 and 0.8 genes
vip_genes <- c(
  "VCP","MIR146A","TF","TSPAN14","MTHFR","DHCR24","DPYSL2","ARC","MIR124-3","NCK2",
  "F2","VSNL1","GSTT1","EPO","FAS","ALDH2","GSTP1","GSTO1","LRPAP1","CALM1",
  "ND2","UQCRC1","PON3","ND1","ECE2","ATP5F1A","CHRNB2","GAPDHS","ADAMTS1","IGF2R",
  "PCDH11X","TPP1","IQCK","MCM2","PCK1","PON2","TPH1","WT1","ABAT","CDK5R1",
  "GSTM3","LIPC","MBL2","NEFM","PIK3R1","PPP3R1","CASP7","IREB2","PPP2R2B","PYY",
  "TPI1","PGRMC1","SLC30A6","SLC2A4","MIR375","SLC30A4","APBB2","PNMT","GSTO2","GRIN2C",
  "HSPA1B"
)



# Keep only genes expressing by > 10% of the cells

assay_use <- "RNA"
layer_use  <- "data"
expr_cut  <- 0
pct_cut   <- 0.1

DefaultAssay(objs_merge_Ast) <- assay_use
vip_genes_in <- intersect(vip_genes, rownames(objs_merge_Ast))
mat <- GetAssayData(objs_merge_Ast, assay = assay_use, layer = layer_use)
pct_expr <- Matrix::rowSums(mat[vip_genes_in, , drop = FALSE] > expr_cut) / ncol(objs_merge_Ast)
vip_keep <- names(pct_expr)[pct_expr > pct_cut]   # 只保留 >10%
vip_pct_table <- data.frame(gene = names(pct_expr), pct_expressed = as.numeric(pct_expr)) |>
  dplyr::arrange(dplyr::desc(pct_expressed))
vip_keep
vip_pct_table

vvip_genes_in <- intersect(vvip_genes, rownames(objs_merge_Ast))
mat <- GetAssayData(objs_merge_Ast, assay = assay_use, layer = layer_use)
pct_expr <- Matrix::rowSums(mat[vvip_genes_in, , drop = FALSE] > expr_cut) / ncol(objs_merge_Ast)
vvip_keep <- names(pct_expr)[pct_expr > pct_cut]   # 只保留 >10%
vvip_pct_table <- data.frame(gene = names(pct_expr), pct_expressed = as.numeric(pct_expr)) |>
  dplyr::arrange(dplyr::desc(pct_expressed))
vvip_keep


avg_mat <- as.matrix(avg_list$RNA)     # genes x subclusters
hvgs <- VariableFeatures(objs_merge_Ast,n=3000)

# 2) Intersect with genes actually present in the AverageExpression matrix
#avg_mat <- as.matrix(avg_list$RNA)           # genes x subclusters
genes_in_avg <- hvgs
prio_all     <- unique(c(vvip_genes, vip_genes))

# union
genes_use <- unique(c(genes_in_avg, prio_all))
genes_use <- intersect(rownames(avg_mat), genes_use)


X  <- avg_mat[genes_use, , drop=FALSE]
Xz <- t(scale(t(X)))
Dat <- t(Xz)

sft <- pickSoftThreshold(
  data         = Dat,
  dataIsExpr   = TRUE,
  networkType  = "unsigned",
  corFnc       = "bicor",
  corOptions   = list(use="pairwise.complete.obs", maxPOutliers=0.05),
  powerVector  = c(2:20),
  verbose      = 1)

beta <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 8

S_cor <- suppressWarnings(
  WGCNA::bicor(t(Xz), use="pairwise.complete.obs", maxPOutliers=0.05)
)
S <- switch("unsigned",
            "signed"        = (S_cor+1)/2,
            "signed hybrid" = pmax(S_cor,0),
            "unsigned"      = abs(S_cor))

S[!is.finite(S)] <- 0
diag(S) <- 0

A <- S^beta
diag(A) <- 0

TOM <- WGCNA::TOMsimilarity(A, TOMType = "unsigned")
dimnames(TOM) <- dimnames(A)
dissTOM <- 1 - TOM


W <- TOM; diag(W) <- 0

library(igraph)
g0 <- graph_from_adjacency_matrix(W, mode="undirected", weighted=TRUE, diag=FALSE)

seed <- setNames(rep(0, nrow(W)), rownames(W))
seed[intersect(names(seed), vip_genes)]  <- 6
seed[intersect(names(seed), vvip_genes)] <- 8

pr <- page_rank(g0, personalized=seed, weights=E(g0)$weight, damping=0.85)$vector
p  <- (pr - min(pr)) / (quantile(pr, 0.95) - min(pr) + 1e-8)
p  <- pmin(p, 1)

alpha <- 0.4
W2 <- W * (1 + alpha * outer(p, p))

s <- 1/sqrt(pmax(rowSums(W2), 1e-8))
names(s) <- rownames(W2)

Wn <- sweep(W2, 1, s, "*")
Wn <- sweep(Wn, 2, s, "*")


Wn <- (Wn + t(Wn))/2
diag(Wn) <- 0

k <- 50
dn <- rownames(Wn); n <- nrow(Wn)


A_knn <- matrix(0, n, n, dimnames=list(dn, dn))
nbrs  <- lapply(1:n, function(i) head(order(Wn[i,], decreasing=TRUE), k+1)[-1])
for (i in 1:n) A_knn[i, nbrs[[i]]] <- 1
A_mut <- (A_knn>0) & (t(A_knn)>0)
Wk    <- Wn * A_mut


diss_k <- 1 - Wk
hc_k   <- if (requireNamespace("fastcluster", quietly=TRUE))
            fastcluster::hclust(as.dist(diss_k), "average") else hclust(as.dist(diss_k), "average")
cl_k   <- cutreeDynamic(dendro=hc_k, distM=as.matrix(diss_k),
                        deepSplit=4, minClusterSize=100, pamRespectsDendro=TRUE)
modules <- split(dn, labels2colors(cl_k))
sizes   <- sort(sapply(modules, length), decreasing=TRUE)
print(head(sizes, 10)); summary(as.numeric(sizes))

hub_brown = modules$brown
brown_genes = modules$brown



# Here we insert the network expension code:
library(matrixStats)
get_module_eigengene <- function(mod_genes, expr_mat, min_genes = 5) {
  g <- intersect(mod_genes, rownames(expr_mat))
  X <- t(scale(t(expr_mat[g, , drop = FALSE]))
  pca <- prcomp(t(X), center = FALSE, scale. = FALSE)
  eig <- pca$x[, 1]
  names(eig) <- colnames(expr_mat)
  return(eig)
}


target_size <- 500

modules_expanded <- modules
mods_to_expand  <- setdiff(names(modules), "grey")

for (mod_name in mods_to_expand) {
  core <- modules[[mod_name]]
  core <- intersect(core, rownames(avg_mat))
  cat("Module", mod_name, "core size =", length(core), "\n")
  # stop if already >500 genes
  if (length(core) >= target_size) {
    modules_expanded[[mod_name]] <- core
    next
  }
  #eigengene
  g <- intersect(core, rownames(avg_mat))

  X <- t(scale(t(avg_mat[g, , drop = FALSE])))
  pca <- prcomp(t(X), center = FALSE, scale. = FALSE)
  ME <- pca$x[, 1]
  names(ME) <- colnames(avg_mat)

  genes_all <- rownames(avg_mat)
  cor_with_ME <- sapply(genes_all, function(gene) {
    x <- avg_mat[gene, ]
    if (all(is.na(x))) return(NA_real_)
    suppressWarnings(
      cor(x, ME, use="pairwise.complete.obs", method="spearman")
    )
  })

  cand <- setdiff(genes_all, core)
  cor_df <- data.frame(
    gene = cand,
    rho  = cor_with_ME[cand],
    stringsAsFactors = FALSE
  )
  cor_df <- subset(cor_df, !is.na(rho))

  cor_df <- cor_df[order(-cor_df$rho), ]

  n_need <- target_size - length(core)
  n_need <- max(0, n_need)
  extra <- head(cor_df$gene, n_need)
  cat("  Add", length(extra), "genes, final size =", length(core) + length(extra), "\n")
  modules_expanded[[mod_name]] <- c(core, extra)
}
sizes_expanded <- sapply(modules_expanded, length)
sizes_expanded
summary(as.numeric(sizes_expanded))



# Module scoring part, I will use "core" genes instead of core+expended genes!
mods_test <- names(modules)
gene_sets <- lapply(mods_test, function(m) intersect(modules[[m]], rownames(objs_merge_Ast)))
names(gene_sets) <- mods_test
gene_sets <- Filter(length, gene_sets)

objs_merge_Ast$AD <- factor(objs_merge_Ast$AD_status, levels=c("Control","AD"))
objs_merge_Ast <- AddModuleScore(objs_merge_Ast, features = gene_sets, name = "MS_")
score_cols <- paste0("MS_", seq_along(gene_sets))
names(score_cols) <- names(gene_sets)

# Wilcoxon（YES vs NO）
ad_fac <- objs_merge_Ast$AD
res <- lapply(names(score_cols), function(m){
  v <- objs_merge_Ast@meta.data[[score_cols[m]]]
  x <- v[ad_fac=="AD"]; y <- v[ad_fac=="Control"]
  data.frame(
    module   = m,
    size     = length(gene_sets[[m]]),
    mean_AD = mean(x, na.rm=TRUE),
    mean_control  = mean(y, na.rm=TRUE),
    diff     = mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE),
    pval     = suppressWarnings(wilcox.test(x, y)$p.value),
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, res)
tab$FDR <- p.adjust(tab$pval, method = "BH")
tab <- tab[order(tab$FDR, -abs(tab$diff)), ]
print(head(tab, 20))


library(writexl)


names(modules_expanded)
sapply(modules_expanded, length)


module_order <- c("blue", "brown", "green", "turquoise", "yellow")
module_order <- intersect(module_order, names(modules_expanded))

modules_to_export <- modules_expanded[module_order]

module_sheets <- lapply(names(modules_to_export), function(m) {
  data.frame(
    module = m,
    gene = as.character(modules_to_export[[m]]),
    stringsAsFactors = FALSE
  )
})
names(module_sheets) <- names(modules_to_export)


summary_sheet <- data.frame(
  module = names(modules_to_export),
  n_genes = as.integer(sapply(modules_to_export, length)),
  stringsAsFactors = FALSE
)

module_sheets <- c(
  list(Summary = summary_sheet),
  module_sheets
)

writexl::write_xlsx(
  module_sheets,
  path = "coexpression_module_genes_expend.xlsx"
)

#Calculate hub genes from brown module



library(WGCNA)
allowWGCNAThreads()

brown_genes <- modules[["brown"]]
brown_genes <- intersect(brown_genes, rownames(avg_mat))  # 保证在 avg_mat 里

Xb <- avg_mat[brown_genes, , drop=FALSE]

sdv <- matrixStats::rowSds(as.matrix(Xb))
Xb <- Xb[sdv > 0, , drop=FALSE]

Xbz <- t(scale(t(Xb)))
Xbz[!is.finite(Xbz)] <- 0
Dat_b <- t(Xbz)

sft_b <- pickSoftThreshold(
  data = Dat_b, dataIsExpr=TRUE,
  networkType="unsigned",
  corFnc="bicor",
  corOptions=list(use="pairwise.complete.obs", maxPOutliers=0.05),
  powerVector=2:20, verbose=0
)
beta_b <- if (!is.na(sft_b$powerEstimate)) sft_b$powerEstimate else 8

S_cor_b <- suppressWarnings(
  WGCNA::bicor(t(Xbz), use="pairwise.complete.obs", maxPOutliers=0.05)
)
S_b <- abs(S_cor_b); S_b[!is.finite(S_b)] <- 0; diag(S_b) <- 0
A_b <- S_b^beta_b; diag(A_b) <- 0

TOM_b <- WGCNA::TOMsimilarity(A_b, TOMType="unsigned")
diag(TOM_b) <- 0
rownames(TOM_b) <- colnames(TOM_b) <- rownames(A_b)

kWithin_brown <- rowSums(TOM_b, na.rm=TRUE)
hub_brown <- sort(kWithin_brown, decreasing=TRUE)
head(hub_brown, 80)

head(names(hub_brown), 80)



# Then the second part
avg_mat <- as.matrix(avg_list$RNA)

genes_all <- unique(res$gene)

module_genes <- unique(modules$brown)
module_genes <- intersect(module_genes, rownames(objs_merge))
library(Seurat)
library(dplyr)
DefaultAssay(objs_merge_Ast) <- "RNA"

objs_merge_Ast <- NormalizeData(objs_merge_Ast, normalization.method = "LogNormalize", scale.factor = 1e4)
objs_merge_Ast <- FindVariableFeatures(objs_merge_Ast, selection.method = "vst", nfeatures = 5000)

max_pcs <- 50
objs_merge_Ast <- ScaleData(objs_merge_Ast, features = rownames(objs_merge_Ast))
objs_merge_Ast <- RunPCA(objs_merge_Ast, features = VariableFeatures(objs_merge_Ast), npcs = max_pcs)

mods_to_pick <- c("brown")
mods_pick <- c("brown")

genes_mod <- unique(unlist(modules_expanded[mods_pick]))
length(genes_mod)

hvg1 <- VariableFeatures(objs_merge_Ast)
features_use1 <- unique(c(genes_mod, hvg1))[1:2500]
setwd("/scratch/hanwu/110525_large_dataset/051926_NCAM2RE_coexp_brown")

length(features_use1)

objs_r1 <- objs_merge_Ast
objs_r1 <- ScaleData(objs_r1, features = features_use1, verbose = FALSE)
objs_r1 <- RunPCA(objs_r1, features = features_use1, npcs = 50)

objs_r1 <- RunHarmony(
  objs_r1,
  group.by.vars = "individualID",
  dims.use = 1:max_pcs,
  assay = DefaultAssay(objs_r1),
  verbose = FALSE
)
pcs_use <- 40
objs_r1 <- FindNeighbors(objs_r1, reduction = "harmony", dims = 1:pcs_use, k.param = 40)
objs_r1 <- FindClusters(objs_r1, resolution = 0.8,verbose = FALSE)  # 0.8 we decided to use 0.8 as the resolution
objs_r1 <- RunUMAP(objs_r1, reduction = "harmony", dims = 1:pcs_use)


objs_merge_Ast$cl_brown_2500 <- Idents(objs_r1)[Cells(objs_merge_Ast)]

Idents(objs_r1) <- "seurat_clusters"
print(unique(objs_r1@meta.data$seurat_clusters))
# Keep in mind we are using CORE module genes for scoring, not extended module.
gene_sets <- modules[mods_pick]
gene_sets <- lapply(gene_sets, function(gs) intersect(gs, rownames(objs_r1)))
objs_r1 <- AddModuleScore(objs_r1, features = gene_sets, name = "MS_")
score_cols_r1 <- paste0("MS_", seq_along(gene_sets))
meta1 <- objs_r1@meta.data %>%
  mutate(cluster = Idents(objs_r1),
         MS_mean = rowMeans(across(all_of(score_cols_r1)))) %>%
  group_by(cluster) %>%
  summarise(n = n(),
            MS_mean = mean(MS_mean),
            across(all_of(score_cols_r1), mean, .names = "avg_{.col}")) %>%
  arrange(desc(MS_mean))
print(meta1)
meta1_filt <- meta1 %>%
  filter(n > 500)
meta1_filt

print(meta1_filt,n=15)



p_umap_cluster <- Seurat::DimPlot(
  objs_r1,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)

p_umap_cluster

ggsave(
      filename = file.path(outdir, paste0("2ndCluster_umap_ALL",".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )

for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_r1, subset = AD_status == grp)
  p_umap_cluster <- Seurat::DimPlot(
  seu_grp,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  ppt_feature_theme +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(title = NULL)
  ggsave(
      filename = file.path(outdir, paste0("2ndCluster_umap_",grp,".png")),
      plot = p_umap_cluster,
      width = 8,
      height = 6,
      dpi = 600,
      bg = "white"
    )
}



genes <- GENELIST

outdir <- OUTDIR

for (gene in genes) {
  expr_all <- GetAssayData(
    object = objs_r1,
    assay = "RNA",
    layer = "data"
  )[gene, ]
  expr_all <- as.numeric(expr_all)
  gene_min <- min(expr_all, na.rm = TRUE)
  gene_max <- max(expr_all, na.rm = TRUE)
  message(gene, ": min = ", gene_min, ", max = ", gene_max)
  p_all <- FeaturePlot(
    object = objs_r1,
    features = gene,
    reduction = "umap",
    order = TRUE,
    pt.size = 0.45,
    layer = "data",
    min.cutoff = gene_min,
    max.cutoff = gene_max
  ) +
    ggtitle(paste(gene, "All cells", sep = " - ")) +
    ppt_feature_theme+  
    ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
  print(p_all)
  ggsave(
    filename = file.path(outdir, paste0("2ndCluster_", gene, "_AllCells_FeaturePlot.png")),
    plot = p_all,
    width = 7,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  ## subset
  for (grp in c("AD", "Control")) {
    seu_grp <- subset(objs_r1, subset = AD_status == grp)
    p <- FeaturePlot(
      object = seu_grp,
      features = gene,
      reduction = "umap",
      order = TRUE,
      pt.size = 0.45,
      layer = "data",
      min.cutoff = gene_min,
      max.cutoff = gene_max
    ) +
      ggtitle(paste(gene, grp, sep = " - ")) +
      ppt_feature_theme + 
      ggplot2::scale_color_gradientn(
        colors = c("grey90", "blue"),
        limits = c(gene_min, gene_max),
        oob = scales::squish,
        name = gene
        )
    print(p)
    ggsave(
      filename = file.path(outdir, paste0(gene, "_", grp, "_2ndCluster", "_FeaturePlot.png")),
      plot = p,
      width = 7,
      height = 6,
      dpi = 600,
      bg = "white"
    )
  }
}

# Violin Plots

outdir <- OUTDIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(objs_r1) <- "RNA"
fill_cols <- c(
  AD = "#D55E00",
  Control = "#0072B2",
  ALL = "#009E73"
)

for (gene in genes) {
  
  expr_all <- as.numeric(expr_mat[gene, colnames(objs_r1)])
  
  gene_min <- min(expr_all, na.rm = TRUE)
  gene_max <- max(expr_all, na.rm = TRUE)
  
  if (gene_min == gene_max) {
    gene_max <- gene_max + 0.1
  }
  
  message(gene, ": min = ", gene_min, ", max = ", gene_max)
  
  for (grp in c("AD", "Control", "ALL")) {
    
    seu_use <- plot_list[[grp]]
    seu_use$plot_group <- factor(seu_use$plot_group, levels = grp)
    
    p <- VlnPlot(
      object = seu_use,
      features = gene,
      assay = "RNA",
      layer = expr_layer,
      group.by = "plot_group",
      pt.size = 0.08,          # 显示细胞小点
      alpha = 0.25,            # 点透明一点，避免太黑
      cols = fill_cols[grp],
      combine = FALSE
    )[[1]] +
      geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        fill = "white",
        color = "black",
        linewidth = 0.3,
        alpha = 0.85
      ) +
      ggtitle(paste(gene, grp, sep = " - ")) +
      xlab(NULL) +
      ylab("Expression") +
      coord_cartesian(ylim = c(gene_min, gene_max)) +
      theme_use +
      NoLegend() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 18,
          face = "bold"
        ),
        axis.title.y = element_text(
          size = 16,
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 13,
          color = "black"
        ),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
    
    print(p)
    
    ggsave(
      filename = file.path(outdir, paste0(gene, "_", grp, "_2ndCluster_ViolinPlot.png")),
      plot = p,
      width = 4,
      height = 5.5,
      dpi = 600,
      bg = "white"
    )
  }
}


# Then split clusters and plot

outdir_cluster <- OUTDIR
dir.create(outdir_cluster, recursive = TRUE, showWarnings = FALSE)


DefaultAssay(objs_r1) <- "RNA"
expr_layer <- "data"

objs_r1$AD_status <- factor(
  objs_r1$AD_status,
  levels = c("AD", "Control")
)

cluster_levels <- sort(unique(as.character(objs_r1$seurat_clusters)))
if (all(grepl("^[0-9]+$", cluster_levels))) {
  cluster_levels <- as.character(sort(as.numeric(cluster_levels)))
}

objs_r1$seurat_clusters <- factor(
  as.character(objs_r1$seurat_clusters),
  levels = cluster_levels
)


expr_mat <- GetAssayData(
  object = objs_r1,
  assay = "RNA",
  layer = expr_layer
)

fill_cols <- c(
  "AD" = "#D55E00",
  "Control" = "#0072B2"
)

theme_use <- if (exists("ppt_feature_theme")) {
  ppt_feature_theme
} else {
  theme_classic(base_size = 16)
}


for (gene in genes) {
  
  expr_all <- as.numeric(expr_mat[gene, colnames(objs_r1)])
  
  gene_min <- min(expr_all, na.rm = TRUE)
  gene_max <- max(expr_all, na.rm = TRUE)
  
  if (gene_min == gene_max) {
    gene_max <- gene_max + 0.1
  }
  
  message(gene, ": min = ", gene_min, ", max = ", gene_max)
  
  # boxplot
  df_box <- data.frame(
    expr = expr_all,
    seurat_clusters = factor(
      as.character(objs_r1$seurat_clusters),
      levels = cluster_levels
    ),
    AD_status = factor(
      as.character(objs_r1$AD_status),
      levels = c("AD", "Control")
    )
  )
  
  p <- VlnPlot(
    object = objs_r1,
    features = gene,
    assay = "RNA",
    layer = expr_layer,
    group.by = "seurat_clusters",
    split.by = "AD_status",
    split.plot = FALSE,
    pt.size = 0.08,
    alpha = 0.25,
    cols = fill_cols,
    combine = FALSE
  )[[1]] +

    geom_boxplot(
      data = df_box,
      aes(
        x = seurat_clusters,
        y = expr,
        group = interaction(seurat_clusters, AD_status)
      ),
      position = position_dodge(width = 0.9),
      width = 0.12,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.3,
      alpha = 0.85,
      inherit.aes = FALSE
    ) +
    ggtitle(paste(gene)) +
    xlab("Seurat cluster") +
    ylab("Expression") +
    coord_cartesian(ylim = c(gene_min, gene_max)) +
    theme_use +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 12,
        color = "black",
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    filename = file.path(
      outdir_cluster,
      paste0(gene, "_2ndCluster_splitAD_ViolinPlot.png")
    ),
    plot = p,
    width = 9,
    height = 5.8,
    dpi = 600,
    bg = "white"
  )
}


##############################

objs <- objs_r1

objs$cluster2 <- as.character(Idents(objs))

objs$AD_status2 <- tolower(as.character(objs$AD_status))

assay_use <- DefaultAssay(objs)
layer_use  <- "data"

print(table(objs$AD_status2))
print(length(unique(objs$cluster2)))


hub_genes_in <- features_use1
mat <- GetAssayData(objs_merge_Ast, assay = assay_use, layer = layer_use)
pct_expr <- Matrix::rowSums(mat[hub_genes_in, , drop = FALSE] > 0) / ncol(objs_merge_Ast)
#hub_keep <- names(pct_expr)[pct_expr > 0]   # 只保留 >10%
avg_expr <- Matrix::rowMeans(mat[hub_genes_in, , drop = FALSE])
hub_pct_table <- data.frame(
  gene = hub_genes_in,
  pct_expressed = as.numeric(pct_expr[hub_genes_in]),
  avg_expression_level = as.numeric(avg_expr[hub_genes_in]),
  row.names = NULL
) |>
  dplyr::arrange(dplyr::desc(pct_expressed))

#hub_keep
head(hub_pct_table)

hub_pct_table2_v1 <- hub_pct_table |>
  dplyr::mutate(
    connectivity = as.numeric(hub_brown[gene])
  ) |>
  dplyr::arrange(dplyr::desc(connectivity))
head(hub_pct_table2_v1)

hub_pct_filt <- hub_pct_table |>
  dplyr::filter(pct_expressed > 0.05)

genes_keep <- hub_pct_filt$gene
length(genes_keep)
head(genes_keep)

hub_genes_in <- genes_keep


get_avg_mat_cells <- function(obj, cells = NULL){
  if (!is.null(cells)) {
    obj <- subset(obj, cells = cells)
  }
  out <- AggregateExpression(
    obj,
    assays = assay_use,
    layer   = layer_use,
    group.by = "cluster2",
    return.seurat = FALSE
  )[[assay_use]]
  as.matrix(out)
}
cells_AD  <- WhichCells(objs, expression = AD_status2 == "ad")
cells_CTL <- WhichCells(objs, expression = AD_status2 == "control")

avg_all <- get_avg_mat_cells(objs, cells = NULL)
avg_AD  <- get_avg_mat_cells(objs, cells = cells_AD)
avg_CTL <- get_avg_mat_cells(objs, cells = cells_CTL)

max(abs(avg_AD - avg_CTL))

common_cl <- Reduce(intersect, list(colnames(avg_all), colnames(avg_AD), colnames(avg_CTL)))

avg_all <- avg_all[, common_cl, drop=FALSE]
avg_AD  <- avg_AD[,  common_cl, drop=FALSE]
avg_CTL <- avg_CTL[, common_cl, drop=FALSE]

cat("common_cl =", length(common_cl), "\n")
cat("genes_keep =", length(genes_keep), "\n")


Xb_all <- avg_all[genes_keep, , drop=FALSE]
Xbz_all <- t(scale(t(Xb_all)))
Xbz_all[!is.finite(Xbz_all)] <- 0
Dat_all <- t(Xbz_all)  # cluster x gene

sft <- pickSoftThreshold(
  data = Dat_all, dataIsExpr=TRUE,
  networkType="unsigned",
  corFnc="bicor",
  corOptions=list(use="pairwise.complete.obs", maxPOutliers=0.05),
  powerVector=2:20, verbose=0
)
beta <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 8
cat("beta =", beta, "\n")

calc_kWithin <- function(avg_mat, genes_keep, beta){
  genes_used <- intersect(genes_keep, rownames(avg_mat))
  Xb <- avg_mat[genes_used, , drop=FALSE]
  if (ncol(Xb) < 3) stop("nClusters < 3")
  if (nrow(Xb) < 2) stop("nGenes < 2")

  Xbz <- t(scale(t(as.matrix(Xb))))
  Xbz[!is.finite(Xbz)] <- 0

  S_cor <- suppressWarnings(
    WGCNA::bicor(t(Xbz), use="pairwise.complete.obs", maxPOutliers=0.05)
  )
  S <- abs(as.matrix(S_cor))
  S[!is.finite(S)] <- 0
  diag(S) <- 0

  A <- S^beta
  diag(A) <- 0
  rownames(A) <- colnames(A) <- rownames(S)

  TOM <- WGCNA::TOMsimilarity(A, TOMType="unsigned")

  TOM <- tryCatch(as.matrix(TOM), error=function(e) TOM)
  if (is.null(dim(TOM))) TOM <- matrix(TOM, nrow=nrow(A), ncol=ncol(A))
  storage.mode(TOM) <- "double"
  TOM[!is.finite(TOM)] <- 0
  diag(TOM) <- 0
  rownames(TOM) <- colnames(TOM) <- rownames(A)

  k <- rowSums(TOM, na.rm=TRUE)
  names(k) <- rownames(TOM)
  k <- sort(k, decreasing=TRUE)
  return(k)
}

calc_pearson_sum_abs <- function(avg_mat, genes_keep){
  genes_used <- intersect(genes_keep, rownames(avg_mat))
  Xb <- avg_mat[genes_used, , drop=FALSE]  # gene x cluster
  if (ncol(Xb) < 3) stop("nClusters < 3")
  if (nrow(Xb) < 2) stop("nGenes < 2")

  Xbz <- t(scale(t(as.matrix(Xb))))
  Xbz[!is.finite(Xbz)] <- 0

  Dat <- t(Xbz)

  Cor_p <- suppressWarnings(cor(Dat, method = "pearson", use="pairwise.complete.obs"))
  Cor_p <- as.matrix(Cor_p)
  Cor_p[!is.finite(Cor_p)] <- 0
  diag(Cor_p) <- 0

  k_p <- rowSums(abs(Cor_p), na.rm = TRUE)
  names(k_p) <- colnames(Cor_p)
  k_p <- sort(k_p, decreasing = TRUE)
  return(k_p)
}


hub_all <- calc_kWithin(avg_all, genes_keep, beta)
hub_AD_r1  <- calc_kWithin(avg_AD,  genes_keep, beta)
hub_CTL_r1 <- calc_kWithin(avg_CTL, genes_keep, beta)

pear_all_abs <- calc_pearson_sum_abs(avg_all, genes_keep)
pear_AD_abs  <- calc_pearson_sum_abs(avg_AD,  genes_keep)
pear_CTL_abs <- calc_pearson_sum_abs(avg_CTL, genes_keep)

mat <- GetAssayData(objs, assay = assay_use, layer = layer_use)

pct_expr <- Matrix::rowSums(mat[hub_genes_in, , drop=FALSE] > 0) / ncol(objs)
avg_expr <- Matrix::rowMeans(mat[hub_genes_in, , drop=FALSE])
#genes_target <- hub_pct_table2_v1$gene
genes_target = hub_genes_in
add_r1 <- tibble::tibble(
  gene = genes_target,
  connectivity_all = as.numeric(hub_all[genes_target]),
  connectivity_AD      = as.numeric(hub_AD_r1[genes_target]),
  connectivity_control = as.numeric(hub_CTL_r1[genes_target]),

  pearson_sum_abs_all = as.numeric(pear_all_abs[genes_target]),
  pearson_sum_abs_AD  = as.numeric(pear_AD_abs[genes_target]),
  pearson_sum_abs_control = as.numeric(pear_CTL_abs[genes_target])
) |>
  dplyr::mutate(
    delta_connectivity = connectivity_AD - connectivity_control,
    log2ratio_connectivity = log2((connectivity_AD + 1e-6) / (connectivity_control + 1e-6)),

    delta_pearson_sum_abs = pearson_sum_abs_AD - pearson_sum_abs_control,
    log2ratio_pearson_sum_abs = log2((pearson_sum_abs_AD + 1e-6) / (pearson_sum_abs_control + 1e-6))
  )

add_r1
hub_pct_plus <- hub_pct_filt |>
  dplyr::left_join(add_r1, by = "gene")

hub_pct_plus_sorted <- hub_pct_plus |>
  dplyr::arrange(dplyr::desc(pearson_sum_abs_all))

genes_keep <- hub_pct_filt$gene
length(genes_keep)
head(genes_keep)

objs <- objs_r1


objs$cluster2 <- as.character(Idents(objs))

objs$AD_status2 <- tolower(as.character(objs$AD_status))

assay_use <- DefaultAssay(objs) 
layer_use  <- "data"

print(table(objs$AD_status2))
print(length(unique(objs$cluster2)))


get_avg_mat <- function(obj){
  out <- AggregateExpression(
    obj,
    assays = assay_use,
    layer   = layer_use,
    group.by = "cluster2",
    return.seurat = FALSE
  )[[assay_use]]
  as.matrix(out)
}

avg_all <- get_avg_mat(objs)
avg_AD  <- get_avg_mat(subset(objs, subset = AD_status2 == "ad"))
avg_CTL <- get_avg_mat(subset(objs, subset = AD_status2 == "control"))

common_cl <- Reduce(intersect, list(colnames(avg_all), colnames(avg_AD), colnames(avg_CTL)))

avg_all <- avg_all[, common_cl, drop=FALSE]
avg_AD  <- avg_AD[,  common_cl, drop=FALSE]
avg_CTL <- avg_CTL[, common_cl, drop=FALSE]

cat("common_cl =", length(common_cl), "\n")

brown_genes = names(hub_brown)
#brown_genes = hub_brown
brown_genes <- intersect(brown_genes, rownames(avg_all))

sd_all <- matrixStats::rowSds(as.matrix(avg_all[brown_genes, , drop=FALSE]))
genes_keep <- brown_genes[sd_all > 0]

cat("genes_keep =", length(genes_keep), "\n")


Xb_all <- avg_all[genes_keep, , drop=FALSE]
Xbz_all <- t(scale(t(Xb_all)))
Xbz_all[!is.finite(Xbz_all)] <- 0
Dat_all <- t(Xbz_all)  # cluster x gene

sft <- pickSoftThreshold(
  data = Dat_all, dataIsExpr=TRUE,
  networkType="unsigned",
  corFnc="bicor",
  corOptions=list(use="pairwise.complete.obs", maxPOutliers=0.05),
  powerVector=2:20, verbose=0
)
beta <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 8
cat("beta =", beta, "\n")

# calculate connectivity
calc_kWithin <- function(avg_mat, genes_keep, beta){
  genes_used <- intersect(genes_keep, rownames(avg_mat))
  Xb <- avg_mat[genes_used, , drop=FALSE]
  if (ncol(Xb) < 3) stop("nClusters < 3")
  if (nrow(Xb) < 2) stop("nGenes < 2")

  Xbz <- t(scale(t(as.matrix(Xb))))
  Xbz[!is.finite(Xbz)] <- 0

  S_cor <- suppressWarnings(
    WGCNA::bicor(t(Xbz), use="pairwise.complete.obs", maxPOutliers=0.05)
  )
  S <- abs(as.matrix(S_cor))
  S[!is.finite(S)] <- 0
  diag(S) <- 0

  A <- S^beta
  diag(A) <- 0
  rownames(A) <- colnames(A) <- rownames(S)

  TOM <- WGCNA::TOMsimilarity(A, TOMType="unsigned")

  TOM <- tryCatch(as.matrix(TOM), error=function(e) TOM)
  if (is.null(dim(TOM))) TOM <- matrix(TOM, nrow=nrow(A), ncol=ncol(A))
  storage.mode(TOM) <- "double"
  TOM[!is.finite(TOM)] <- 0
  diag(TOM) <- 0
  rownames(TOM) <- colnames(TOM) <- rownames(A)

  k <- rowSums(TOM, na.rm=TRUE)
  names(k) <- rownames(TOM)
  k <- sort(k, decreasing=TRUE)
  return(k)
}

hub_all <- calc_kWithin(avg_all, genes_keep, beta)
hub_AD_r1  <- calc_kWithin(avg_AD,  genes_keep, beta)
hub_CTL_r1 <- calc_kWithin(avg_CTL, genes_keep, beta)

# tables
hub_genes_in <- head(names(hub_all), 196)
hub_genes_in <- intersect(hub_genes_in, rownames(GetAssayData(objs, assay=assay_use, layer=layer_use)))

mat <- GetAssayData(objs, assay = assay_use, layer = layer_use)

pct_expr <- Matrix::rowSums(mat[hub_genes_in, , drop=FALSE] > 0) / ncol(objs)
avg_expr <- Matrix::rowMeans(mat[hub_genes_in, , drop=FALSE])
genes_target <- hub_pct_table2_v1$gene

add_r1 <- tibble::tibble(
  gene = genes_target,
  connectivity_r1_AD      = as.numeric(hub_AD_r1[genes_target]),
  connectivity_r1_control = as.numeric(hub_CTL_r1[genes_target])
) |>
  dplyr::mutate(
    delta_connectivity_r1 = connectivity_r1_AD - connectivity_r1_control,
    log2ratio_connectivity_r1 = log2((connectivity_r1_AD + 1e-6) / (connectivity_r1_control + 1e-6))
  )

hub_pct_table2_v1_plus <- hub_pct_table2_v1 |>
  dplyr::left_join(add_r1, by = "gene")

head(hub_pct_table2_v1_plus)


hub_pct_table2_v1_plus_clean <- hub_pct_table2_v1_plus %>%
  dplyr::select(
    -connectivity
  ) %>%
  dplyr::rename(
    connectivity_AD = connectivity_r1_AD,
    connectivity_control = connectivity_r1_control,
    delta_connectivity = delta_connectivity_r1,
    log2ratio_connectivity = log2ratio_connectivity_r1
  )


hub_connectivity_export <- hub_pct_table2_v1_plus %>%
  dplyr::select(
    -dplyr::any_of("connectivity")
  ) %>%
  dplyr::rename(
    connectivity_AD = connectivity_r1_AD,
    connectivity_control = connectivity_r1_control,
    delta_connectivity = delta_connectivity_r1,
    log2ratio_connectivity = log2ratio_connectivity_r1
  ) %>%
  dplyr::mutate(
    gene_class = dplyr::case_when(
      gene %in% vvip_genes ~ "GDA ≥ 0.8",
      gene %in% vip_genes ~ "0.6 < GDA < 0.8",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::relocate(
    gene_class,
    .after = gene
  )

# write data in xlsx format

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(
  wb,
  sheetName = "AD_control_connectivity"
)

openxlsx::writeData(
  wb,
  sheet = "AD_control_connectivity",
  x = hub_connectivity_export,
  withFilter = TRUE
)


header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  fgFill = "#FFFFFF",
  border = "Bottom",
  halign = "center",
  valign = "center"
)

vvip_style <- openxlsx::createStyle(
  fgFill = "#F4A3B4",   # soft pink/red
  fontColour = "#000000"
)

vip_style <- openxlsx::createStyle(
  fgFill = "#F6D365",   # soft yellow
  fontColour = "#000000"
)

num_style <- openxlsx::createStyle(
  numFmt = "0.000"
)

sheet_name <- "AD_control_connectivity"

openxlsx::addStyle(
  wb,
  sheet = sheet_name,
  style = header_style,
  rows = 1,
  cols = 1:ncol(hub_connectivity_export),
  gridExpand = TRUE,
  stack = TRUE
)


gene_col <- which(colnames(hub_connectivity_export) == "gene")

vvip_rows <- which(hub_connectivity_export$gene %in% vvip_genes) + 1
vip_rows <- which(
  hub_connectivity_export$gene %in% vip_genes &
    !(hub_connectivity_export$gene %in% vvip_genes)
) + 1

if (length(vvip_rows) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = vvip_style,
    rows = vvip_rows,
    cols = gene_col,
    gridExpand = TRUE,
    stack = TRUE
  )
}

if (length(vip_rows) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = vip_style,
    rows = vip_rows,
    cols = gene_col,
    gridExpand = TRUE,
    stack = TRUE
  )
}

numeric_cols <- which(sapply(hub_connectivity_export, is.numeric))

if (length(numeric_cols) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = num_style,
    rows = 2:(nrow(hub_connectivity_export) + 1),
    cols = numeric_cols,
    gridExpand = TRUE,
    stack = TRUE
  )
}

openxlsx::freezePane(
  wb,
  sheet = sheet_name,
  firstRow = TRUE
)

openxlsx::setColWidths(
  wb,
  sheet = sheet_name,
  cols = 1:ncol(hub_connectivity_export),
  widths = "auto"
)

openxlsx::saveWorkbook(
  wb,
  file = "hub_pct_table2_v1_plus_AD_control_connectivity_colored.xlsx",
  overwrite = TRUE
)


# count cell number from 2nd round clusters
## cluster x condition count table
cluster_count_mat <- table(
  objs@meta.data$cluster2,
  objs@meta.data$AD_status2
)

cluster_count_mat
cluster_enrichment_simple <- lapply(rownames(cluster_count_mat), function(cl) {
  
  in_AD <- cluster_count_mat[cl, "ad"]
  in_control <- cluster_count_mat[cl, "control"]
  
  out_AD <- sum(cluster_count_mat[, "ad"]) - in_AD
  out_control <- sum(cluster_count_mat[, "control"]) - in_control
  
  fisher_mat <- matrix(
    c(
      in_AD, out_AD,
      in_control, out_control
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(fisher_mat) <- c("AD", "Control")
  colnames(fisher_mat) <- c("In_cluster", "Not_in_cluster")
  
  fisher_res <- fisher.test(fisher_mat)
  
  data.frame(
    cluster = cl,
    AD_cells = as.integer(in_AD),
    control_cells = as.integer(in_control),
    odds_ratio = as.numeric(fisher_res$estimate),
    p_value = fisher_res$p.value,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    enrichment = dplyr::case_when(
      odds_ratio > 1 ~ "AD enriched",
      odds_ratio < 1 ~ "Control enriched",
      TRUE ~ "No bias"
    )
  ) |>
  dplyr::arrange(dplyr::desc(odds_ratio))

cluster_enrichment_simple

names(hub_brown)


# 2nd cluster
objs <- objs_r1
objs$cluster2 <- as.character(Idents(objs))
objs$AD_status2 <- tolower(as.character(objs$AD_status))

assay_use <- DefaultAssay(objs)
layer_use  <- "data"

table(objs$AD_status2)
length(unique(objs$cluster2))

get_avg_mat <- function(obj){
  as.matrix(AggregateExpression(
    obj,
    assays = assay_use,
    layer   = layer_use,
    group.by = "cluster2",
    return.seurat = FALSE
  )[[assay_use]])
}

avg_AD  <- get_avg_mat(subset(objs, subset = AD_status2 == "ad"))
avg_CTL <- get_avg_mat(subset(objs, subset = AD_status2 == "control"))


avg_ALL <- get_avg_mat(objs)

avg_AD  <- avg_AD[,  common_cl, drop=FALSE]
avg_CTL <- avg_CTL[, common_cl, drop=FALSE]
avg_ALL <- avg_ALL[, common_cl, drop=FALSE]   # 新增这一行（保持同一批 cluster）

cat("common clusters =", length(common_cl), "\n")


genes <- names(hub_brown)

cat("genes used =", length(genes), "\n")
genes_keep = genes

cat("genes_keep =", length(genes_keep), "\n")

beta <- 8

calc_TOM <- function(avg_mat, genes_keep, beta){
  X <- avg_mat[genes_keep, , drop=FALSE]
  Xz <- t(scale(t(as.matrix(X))))
  Xz[!is.finite(Xz)] <- 0

  R <- suppressWarnings(WGCNA::bicor(
    t(Xz), use="pairwise.complete.obs", maxPOutliers=0.05
  ))
  R <- as.matrix(R)
  R[!is.finite(R)] <- 0
  diag(R) <- 0
  rownames(R) <- colnames(R) <- rownames(X)

  A <- abs(R)^beta
  diag(A) <- 0
  rownames(A) <- colnames(A) <- rownames(R)

  TOM <- WGCNA::TOMsimilarity(A, TOMType="unsigned")
  TOM <- as.matrix(TOM)
  TOM[!is.finite(TOM)] <- 0
  diag(TOM) <- 0
  rownames(TOM) <- colnames(TOM) <- rownames(A)

  list(R=R, A=A, TOM=TOM)
}

net_AD  <- calc_TOM(avg_AD,  genes_keep, beta)
net_CTL <- calc_TOM(avg_CTL, genes_keep, beta)
net_ALL <- calc_TOM(avg_ALL, genes_keep, beta) 

TOM_AD  <- net_AD$TOM
TOM_CTL <- net_CTL$TOM
TOM_ALL <- net_ALL$TOM
dTOM    <- TOM_AD - TOM_CTL

# not removing any genes
pct_min <- 0

expr_tbl <- hub_pct_table2_v1_plus %>%
  dplyr::select(gene, pct_expressed, avg_expression_level)

genes_plot <- expr_tbl %>%
  filter(pct_expressed >= pct_min) %>%
  pull(gene)


genes_plot <- intersect(genes_plot, rownames(TOM_ALL))
genes_plot <- intersect(genes_plot, genes_keep)


TOM_AD_full  <- TOM_AD
TOM_CTL_full <- TOM_CTL
TOM_ALL_full <- TOM_ALL
dTOM_full    <- dTOM


TOM_AD  <- TOM_AD[genes_plot, genes_plot, drop=FALSE]
TOM_CTL <- TOM_CTL[genes_plot, genes_plot, drop=FALSE]
TOM_ALL <- TOM_ALL[genes_plot, genes_plot, drop=FALSE]
dTOM    <- TOM_AD - TOM_CTL
diag(dTOM) <- 0

node_tbl_all <- tibble(
  gene = rownames(TOM_ALL),
  connectivity_ALL = rowSums(TOM_ALL, na.rm=TRUE)
) %>%
  left_join(expr_tbl, by="gene") %>%
  arrange(desc(connectivity_ALL))

ut_all <- which(upper.tri(TOM_ALL), arr.ind = TRUE)
g_all  <- rownames(TOM_ALL)

edge_tbl_all <- tibble(
  gene1 = g_all[ut_all[,1]],
  gene2 = g_all[ut_all[,2]],
  tom_ALL = as.numeric(TOM_ALL[ut_all])
) %>%
  arrange(desc(tom_ALL))

head(node_tbl_all, 10)
head(edge_tbl_all, 10)

node_tbl <- tibble(
  gene = rownames(dTOM),
  connectivity_AD      = rowSums(TOM_AD,  na.rm=TRUE),
  connectivity_control = rowSums(TOM_CTL, na.rm=TRUE)
) |>
  mutate(
    delta_connectivity = connectivity_AD - connectivity_control,
    log2ratio_connectivity = log2((connectivity_AD + 1e-6) / (connectivity_control + 1e-6))
  ) |>
  arrange(desc(abs(delta_connectivity)))

head(node_tbl, 20)


ut <- which(upper.tri(dTOM), arr.ind = TRUE)
g  <- rownames(dTOM)

edge_tbl <- tibble(
  gene1 = g[ut[,1]],
  gene2 = g[ut[,2]],
  tom_AD  = TOM_AD[ut],
  tom_control = TOM_CTL[ut],
  tom_ALL = TOM_ALL[ut],
  delta_tom = dTOM[ut],
  abs_delta_tom = abs(dTOM[ut])
) |>
  arrange(desc(abs_delta_tom))

head(edge_tbl, 20)


rewire_strength <- rowSums(abs(dTOM), na.rm=TRUE)
thr <- 0.10
rewire_count <- rowSums(abs(dTOM) >= thr, na.rm=TRUE)

rewire_tbl <- tibble(
  gene = names(rewire_strength),
  rewire_strength = as.numeric(rewire_strength),
  rewire_count = as.numeric(rewire_count)
) |>
  arrange(desc(rewire_strength))

head(rewire_tbl, 20)

node_tbl2 <- node_tbl |>
  left_join(rewire_tbl, by="gene") |>
  arrange(desc(rewire_strength), desc(abs(log2ratio_connectivity)))

head(node_tbl2, 30)


dTOM <- TOM_AD - TOM_CTL
diag(dTOM) <- 0

ut <- which(upper.tri(dTOM), arr.ind = TRUE)
g  <- rownames(dTOM)

edge_rank <- tibble(
  gene1 = g[ut[,1]],
  gene2 = g[ut[,2]],
  tom_AD  = as.numeric(TOM_AD[ut]),
  tom_control = as.numeric(TOM_CTL[ut]),
  tom_ALL = as.numeric(TOM_ALL[ut]),
  delta_tom = as.numeric(dTOM[ut]),
  direction = ifelse(delta_tom > 0, "AD_stronger", "control_stronger")
) %>%
  arrange(desc(delta_tom))

head(edge_rank, 30)

expr1 <- hub_pct_table2_v1_plus
colnames(expr1)[colnames(expr1) == "gene"] <- "gene1"
colnames(expr1)[colnames(expr1) == "pct_expressed"] <- "pct1"
colnames(expr1)[colnames(expr1) == "avg_expression_level"] <- "avg1"

expr2 <- hub_pct_table2_v1_plus
colnames(expr2)[colnames(expr2) == "gene"] <- "gene2"
colnames(expr2)[colnames(expr2) == "pct_expressed"] <- "pct2"
colnames(expr2)[colnames(expr2) == "avg_expression_level"] <- "avg2"


edge_rank_w <- merge(edge_rank, expr1, by = "gene1", all.x = TRUE)
edge_rank_w <- merge(edge_rank_w, expr2, by = "gene2", all.x = TRUE)

edge_rank_w$edge_score <- edge_rank_w$delta_tom * sqrt(edge_rank_w$pct1 * edge_rank_w$pct2)
edge_rank_w <- edge_rank_w[order(edge_rank_w$edge_score, decreasing = TRUE), ]

head(edge_rank_w, 20)
nrow(edge_rank_w)

write_xlsx(edge_rank_w,
           "edge_rank_w_052026.xlsx")


edge_rank_w2 <- edge_rank_w %>%
  dplyr::mutate(
    abs_edge_score = abs(edge_score)
  ) %>%
  dplyr::arrange(dplyr::desc(abs_edge_score))
head(edge_rank_w2, 40)



# below is the code for creating the network plot for top 100 edges

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(igraph)
library(scales)


top_n <- 100
out_prefix <- paste0("top", top_n, "_abs_edge_score_network")

## Edge colors
ad_edge_col  <- "#D55E00"
ctl_edge_col <- "#0072B2"

## Node colors: priority class
## node colors
vvip_col  <-"#F94C66"
vip_col   <- "#53BF9D"  
other_col <- "#E5E7EB"
## Layout parameters
layout_niter <- 3000
graphopt_charge <- 0.04
graphopt_mass <- 70
graphopt_spring_length <- 3.2
graphopt_spring_constant <- 0.6

set.seed(123)


vip_genes <- unique(as.character(vip_genes))
vvip_genes <- unique(as.character(vvip_genes))


top_edges <- edge_rank_w2 %>%
  dplyr::mutate(
    direction = trimws(as.character(direction)),
    direction = factor(direction, levels = c("AD_stronger", "control_stronger")),
    abs_edge_score = abs(edge_score),
    edge_name = paste(gene1, gene2, sep = "--")
  ) %>%
  dplyr::arrange(dplyr::desc(abs_edge_score)) %>%
  dplyr::slice_head(n = top_n)

cat("Number of selected edges:", nrow(top_edges), "\n")
cat("Direction counts:\n")
print(table(top_edges$direction, useNA = "ifany"))

node_names <- unique(c(top_edges$gene1, top_edges$gene2))
cat("Unique genes:", length(node_names), "\n")


nodes <- tibble::tibble(
  name = node_names
) %>%
  dplyr::left_join(
    node_tbl2 %>%
      dplyr::select(
        gene,
        connectivity_AD,
        connectivity_control,
        delta_connectivity
      ),
    by = c("name" = "gene")
  ) %>%
  dplyr::mutate(
    connectivity_ALL = (connectivity_AD + connectivity_control) / 2,

    node_class = dplyr::case_when(
      name %in% vvip_genes ~ "GDA ≧ 0.8 genes",
      name %in% vip_genes  ~ "0.6 < GDA < 0.8 genes",
      TRUE ~ "Other"
    ),
    node_class = factor(node_class, levels = c("GDA ≧ 0.8 genes", "0.6 < GDA < 0.8 genes", "Other"))
  )

degree_tbl <- dplyr::bind_rows(
  top_edges %>%
    dplyr::count(gene1, name = "degree") %>%
    dplyr::rename(name = gene1),
  top_edges %>%
    dplyr::count(gene2, name = "degree") %>%
    dplyr::rename(name = gene2)
) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(degree_top = sum(degree), .groups = "drop")

nodes <- nodes %>%
  dplyr::left_join(degree_tbl, by = "name") %>%
  dplyr::mutate(
    degree_top = ifelse(is.na(degree_top), 0, degree_top),
    connectivity_ALL = ifelse(
      is.na(connectivity_ALL),
      median(connectivity_ALL, na.rm = TRUE),
      connectivity_ALL
    )
  )

edges_for_graph <- top_edges %>%
  dplyr::transmute(
    from = gene1,
    to = gene2,
    gene1 = gene1,
    gene2 = gene2,
    edge_name = edge_name,
    tom_AD = tom_AD,
    tom_control = tom_control,
    tom_ALL = tom_ALL,
    delta_tom = delta_tom,
    direction = direction,
    edge_score = edge_score,
    abs_edge_score = abs_edge_score
  )



g <- igraph::graph_from_data_frame(
  d = edges_for_graph %>% dplyr::select(from, to),
  vertices = nodes,
  directed = FALSE
)

layout_mat <- igraph::layout_with_graphopt(
  graph = g,
  niter = layout_niter,
  charge = graphopt_charge,
  mass = graphopt_mass,
  spring.length = graphopt_spring_length,
  spring.constant = graphopt_spring_constant
)

node_coord <- tibble::tibble(
  name = igraph::V(g)$name,
  x = layout_mat[, 1],
  y = layout_mat[, 2]
) %>%
  dplyr::left_join(nodes, by = "name")


edge_plot <- edges_for_graph %>%
  dplyr::left_join(
    node_coord %>% dplyr::select(name, x, y),
    by = c("gene1" = "name")
  ) %>%
  dplyr::rename(x1 = x, y1 = y) %>%
  dplyr::left_join(
    node_coord %>% dplyr::select(name, x, y),
    by = c("gene2" = "name")
  ) %>%
  dplyr::rename(x2 = x, y2 = y)


plot_network_absEdgeScore <- function(edge_df,
                                      node_df,
                                      title_text = NULL) {
  
  edge_df2 <- edge_df %>%
    dplyr::mutate(
      direction = factor(
        as.character(direction),
        levels = c("AD_stronger", "control_stronger")
      )
    )
  
  node_df2 <- node_df %>%
    dplyr::mutate(
      node_class_plot = dplyr::case_when(
        node_class %in% c("VVIP", "GDA ≧ 0.8 genes", "GDA >= 0.8 genes", "GDA > 0.8") ~ "GDA ≧ 0.8 genes",
        node_class %in% c("VIP", "0.6 < GDA < 0.8 genes", "0.6 < GDA <= 0.8 genes", "0.6 < GDA ≤ 0.8 genes") ~ "0.6 < GDA < 0.8 genes",
        TRUE ~ "Other"
      ),
      node_class_plot = factor(
        node_class_plot,
        levels = c("GDA ≧ 0.8 genes", "0.6 < GDA < 0.8 genes", "Other")
      )
    )
  
  n_vvip_plot <- sum(node_df2$node_class_plot == "GDA ≧ 0.8 genes", na.rm = TRUE)
  n_vip_plot  <- sum(node_df2$node_class_plot == "0.6 < GDA < 0.8 genes", na.rm = TRUE)
  
  node_class_labels <- c(
    "GDA ≧ 0.8 genes" = paste0("GDA ≧ 0.8 genes (n = ", n_vvip_plot, ")"),
    "0.6 < GDA < 0.8 genes" = paste0("0.6 < GDA < 0.8 genes (n = ", n_vip_plot, ")"),
    "Other" = "Other"
  )
  
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df2,
      ggplot2::aes(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2,
        linewidth = abs_edge_score,
        color = direction
      ),
      alpha = 0.72,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = node_df2,
      ggplot2::aes(
        x = x,
        y = y,
        size = degree_top,
        fill = node_class_plot
      ),
      shape = 21,
      color = "black",
      stroke = 0.7
    ) +
    ggrepel::geom_text_repel(
      data = node_df2,
      ggplot2::aes(
        x = x,
        y = y,
        label = name
      ),
      size = 10.8,
      max.overlaps = Inf,
      box.padding = 1.2,
      point.padding = 0.6,
      force = 5,
      force_pull = 0.25,
      segment.size = 0.25,
      min.segment.length = 0,
      seed = 123
    ) +
    ggplot2::scale_color_manual(
      name = "Direction",
      values = c(
        "AD_stronger" = ad_edge_col,
        "control_stronger" = ctl_edge_col
      ),
      breaks = c("AD_stronger", "control_stronger"),
      labels = c("AD stronger", "Control stronger"),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = "Node class",
      values = c(
        "GDA ≧ 0.8 genes" = vvip_col,
        "0.6 < GDA < 0.8 genes" = vip_col,
        "Other" = other_col
      ),
      breaks = c(
        "GDA ≧ 0.8 genes",
        "0.6 < GDA < 0.8 genes",
        "Other"
      ),
      labels = node_class_labels,
      drop = FALSE
    ) +
    ggplot2::scale_linewidth(
      name = "|edge score|",
      range = c(0.25, 4.2)
    ) +
    ggplot2::scale_size(
      name = "Top-edge degree",
      range = c(4.2, 12)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          linewidth = 8,
          alpha = 1
        ),
        keywidth = grid::unit(2.4, "cm"),
        keyheight = grid::unit(1.2, "cm"),
        order = 1
      ),
      fill = ggplot2::guide_legend(
        override.aes = list(
          size = 18,
          shape = 21,
          color = "black",
          stroke = 1
        ),
        keywidth = grid::unit(1.6, "cm"),
        keyheight = grid::unit(1.3, "cm"),
        order = 2
      ),
      linewidth = ggplot2::guide_legend(
        keywidth = grid::unit(2.4, "cm"),
        keyheight = grid::unit(1.2, "cm"),
        order = 3
      ),
      size = ggplot2::guide_legend(
        override.aes = list(
          size = c(10, 14, 18)
        ),
        keywidth = grid::unit(1.6, "cm"),
        keyheight = grid::unit(1.3, "cm"),
        order = 4
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = NULL) +
    ggplot2::theme_void(base_size = 24) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = 30,
        color = "black",
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 27,
        color = "black"
      ),
      legend.key.size = grid::unit(1.5, "cm"),
      legend.spacing.y = grid::unit(1.0, "cm"),
      legend.box.spacing = grid::unit(0.8, "cm"),
      legend.margin = ggplot2::margin(12, 12, 12, 12),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

p_absEdgeScore <- plot_network_absEdgeScore(
  edge_df = edge_plot,
  node_df = node_coord,
  title_text = paste0("Top ", top_n, " rewired edges")
)

p_absEdgeScore

ggsave(
  filename = paste0(out_prefix, "_absEdgeScore_VIPnode_network_100.png"),
  plot = p_absEdgeScore,
  width = 20,
  height = 15,
  units = "in",
  dpi = 600,
  bg = "white"
)
