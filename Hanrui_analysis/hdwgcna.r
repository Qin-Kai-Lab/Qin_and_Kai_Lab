.libPaths(LIBPATH)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(harmony)
library(SeuratDisk)
library(readxl)
library(readr)
library(writexl)
library(tibble)
library(scCustomize)
library(jsonlite)
library(stringr)
library(purrr)
library(hdWGCNA)
library(gridExtra)
library(grid)
library(cowplot)
library(openxlsx)
library(UCell)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(AnnotationHub)
library(ensembldb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)

setwd(WORKDIR)

astro <- LoadH5Seurat(
  "astrocytes.h5Seurat",
  assays     = list(RNA = "counts"),
  reductions = FALSE,
  graphs     = FALSE,
  images     = FALSE,
  meta.data  = TRUE
)

astro
GetAssayData(astro,layer = "counts") 
DefaultAssay(astro)
head(astro@meta.data)
dim(astro)

xlsx_path <- "ROSMAP_ID_Status.xlsx" 

samp <- read_excel("ROSMAP_ID_Status.xlsx",sheet = 1, na = c("NA","N/A",""))

samp$individualID <- trimws(as.character(samp$individualID))
astro@meta.data$individualID <- trimws(as.character(astro@meta.data$individualID))

samp$individualID_excel = samp$individualID
samp$individualID = NULL
names(samp)[!names(samp) %in% "individualID_excel"] <-
  paste0("clin_", names(samp)[!names(samp) %in% "individualID_excel"])
cell_barcodes <- rownames(astro@meta.data)

astro@meta.data <- astro@meta.data %>%
  mutate(individualID_cell = individualID) %>%
  left_join(samp, by = c("individualID_cell" = "individualID_excel"))


rownames(astro@meta.data) <- cell_barcodes


# Sample selection:
astro@meta.data <- astro@meta.data %>%
  mutate(
    AD_status = case_when(
      clin_braaksc %in% c(0, 1, 2) &
        clin_ceradsc == 4 &
        clin_cogdx == 1 &
        clin_dcfdx_lv == 1 ~ "control",
      
      clin_braaksc %in% c(5, 6) &
        clin_ceradsc == 1 &
        clin_cogdx == 4 &
        clin_dcfdx_lv == 4 ~ "AD",
    
      TRUE ~ NA_character_  
    )
  )

astro_selected_samples <- subset(astro, subset = AD_status %in% c("AD","control"))

astro = astro_selected_samples

all_genes <- rownames(astro)

genes_all <- all_genes
anno <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys=all_genes,
                        keytype="SYMBOL", columns=c("SYMBOL","CHR"))
)
# Remove unnecessary genes
anno2 <- anno %>%
  filter(!is.na(CHR)) %>%
  mutate(CHR = gsub("^chr","", CHR, ignore.case = TRUE)) %>%
  filter(CHR %in% as.character(1:22)) %>%
  distinct(SYMBOL, .keep_all = TRUE)
auto_genes <- intersect(all_genes, anno2$SYMBOL)
sex_or_nonauto <- setdiff(all_genes, auto_genes)  # only keeps 1-22 chr genes
length(auto_genes); length(sex_or_nonauto); head(sex_or_nonauto)
astro <- subset(astro, features = auto_genes)

# cluster Fisher enrichment
cluster_fisher_table <- function(seu,
                                 cluster_col = "seurat_clusters",
                                 status_col  = "AD_status",
                                 case = "AD",
                                 ctrl = "control") {
  md <- seu@meta.data
  stopifnot(cluster_col %in% colnames(md), status_col %in% colnames(md))

  md <- md %>% dplyr::filter(.data[[status_col]] %in% c(case, ctrl))
  md[[cluster_col]] <- as.character(md[[cluster_col]])
  clus <- sort(unique(md[[cluster_col]]))

  out <- lapply(clus, function(cl) {
    in_cl <- md[[cluster_col]] == cl
    a <- sum(in_cl  & md[[status_col]] == case)
    b <- sum(in_cl  & md[[status_col]] == ctrl)
    c <- sum(!in_cl & md[[status_col]] == case)
    d <- sum(!in_cl & md[[status_col]] == ctrl)

    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    data.frame(
      cluster = cl,
      AD = a,
      control = b,
      AD_pct = ifelse((a + b) > 0, round(a/(a + b), 3), NA_real_),
      OR = unname(ft$estimate),
      p = ft$p.value,
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()

  out$FDR <- p.adjust(out$p, method = "fdr")
  out %>% dplyr::arrange(FDR, dplyr::desc(OR))
}

set.seed(123)

objs_astro <- astro
DefaultAssay(objs_astro) <- "RNA"

objs_astro <- NormalizeData(
  objs_astro,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

objs_astro <- FindVariableFeatures(
  objs_astro,
  selection.method = "vst",
  nfeatures = 5000,
  verbose = FALSE
)

hvg5000 <- VariableFeatures(objs_astro)

objs_astro <- ScaleData(
  objs_astro,
  features = hvg5000,
  verbose = FALSE
)

objs_astro <- RunPCA(
  objs_astro,
  features = hvg5000,
  npcs = 80,
  seed.use = 123,
  verbose = FALSE
)

## SCAN PCS
pcs_vec <- seq(30, 60, by = 5)
k_param <- 40
res_use <- 0.8

pdf_file <- sprintf("PC_scan_res%.1f_k%d_noSex_031826_seed123_final.pdf", res_use, k_param)
pdf(pdf_file, width = 14, height = 6, onefile = TRUE)

for (pcs_use in pcs_vec) {
  message("Running pcs_use = ", pcs_use)
  set.seed(123)
  seu <- objs_astro
  ## harmony
  seu <- harmony::RunHarmony(
    object = seu,
    group.by.vars = "individualID",
    reduction.use = "pca",
    dims.use = 1:pcs_use,
    reduction.save = "harmony",
    project.dim = TRUE,
    verbose = FALSE
  )

  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:pcs_use, k.param = k_param, verbose = FALSE)
  seu <- FindClusters(seu, resolution = res_use, random.seed = 123, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:pcs_use,seed.use = 123, verbose = FALSE)
  print(table(seu@meta.data$seurat_clusters,seu@meta.data$AD_status))

  p_umap <- DimPlot(
    seu,
    reduction = "umap",
    split.by = "AD_status",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0(
      "pcs_use=", pcs_use,
      " | res=", res_use,
      " | k=", k_param,
      " | clusters=", length(unique(seu$seurat_clusters)),
      " | min_cluster_n=", min(table(seu$seurat_clusters))
    )) +
    theme(plot.title = element_text(size = 12))
  # also add summary
  tab <- cluster_fisher_table(
    seu,
    cluster_col = "seurat_clusters",
    status_col  = "AD_status",
    case = "AD",
    ctrl = "control"
  ) %>%
    dplyr::mutate(dplyr::across(c(OR, p, FDR), ~ signif(.x, 3)))
  tg <- gridExtra::tableGrob(tab, rows = NULL)
  title_grob <- grid::textGrob(
    "Cell-level Fisher (enrichment)",
    gp = grid::gpar(fontsize = 11, fontface = "bold")
  )
  right_panel <- gridExtra::arrangeGrob(
    title_grob, tg,
    ncol = 1,
    heights = grid::unit(c(0.8, 9.2), "null")
  )
  combined <- gridExtra::arrangeGrob(
    ggplotGrob(p_umap),
    right_panel,
    ncol = 2,
    widths = grid::unit(c(3.2, 1.8), "null")
  )
  grid::grid.newpage()
  grid::grid.draw(combined)
}
dev.off()
cat("Wrote PDF:", pdf_file, "\n")

# we decided to use PC = 50
# PC = 50
set.seed(123)
seu <- objs_astro
k_param <- 40
res_use <- 0.8
pcs_use = 50

seu <- harmony::RunHarmony(
    object = seu,
    group.by.vars = "individualID",
    reduction.use = "pca",
    dims.use = 1:pcs_use,
    reduction.save = "harmony",
    project.dim = TRUE,
    verbose = FALSE
  )

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:pcs_use, k.param = k_param, verbose = FALSE)
seu <- FindClusters(seu, resolution = res_use, random.seed = 123,verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:pcs_use, seed.use = 123,verbose = FALSE)

seu

p_a <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "AD_status",
  label = TRUE,
  repel = TRUE
)

p_b <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

pdf("DimPlot_umap_seed123_final.pdf", width = 12, height = 6)
print(p_a)
print(p_b)
dev.off()

# Then DEG

clusters <- sort(unique(seu@meta.data$seurat_clusters))
wb <- createWorkbook()

for (cl in clusters) {
  obj_cl <- subset(seu, subset = seurat_clusters == cl)
  Idents(obj_cl) <- "AD_status"
  deg <- FindMarkers(
    obj_cl,
    ident.1 = "AD",
    ident.2 = "control",
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  # first column gene 
  deg <- deg %>%
    rownames_to_column("gene") %>%
    mutate(FDR_BH = p.adjust(p_val, method = "BH")) %>%
    relocate(gene, .before = 1)
  # deg <- deg %>% rename(logFC = avg_log2FC)

  sheet_name <- paste0("cluster_", cl)
  sheet_name <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", sheet_name), 1, 31)

  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, deg, tableStyle = "TableStyleLight9")

  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(deg), widths = "auto")
}
saveWorkbook(wb, file = "DEG_by_cluster_AD_vs_control_NoSexMT_pc50_res08_k40_031826Seed123.xlsx", overwrite = TRUE)


# update seurat object to WGCNA format

seu$ALL <- "ALL"
seurat_obj <- SeuratObject::UpdateSeuratObject(seu)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Astrocyte031826"
)


seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters", "AD_status","ALL"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "seurat_clusters"
)

seurat_obj <- NormalizeMetacells(seurat_obj)
mc_ALL <- seurat_obj

sparse_rowVars <- function(M) {
  mu <- Matrix::rowMeans(M)
  ex2 <- Matrix::rowMeans(M^2)
  as.numeric(ex2 - mu^2)
}

mat_ALL <- GetAssayData(mc_ALL, assay = "RNA", layer = "data")
det_ALL <- Matrix::rowMeans(mat_ALL > 0)
var_ALL <- sparse_rowVars(mat_ALL)

keep <- (det_ALL >= 0.05) & (var_ALL > 0)

genes_keep <- rownames(mat_ALL)[keep]

mc_ALL <- subset(mc_ALL, features = genes_keep)

mc_ALL <- SetDatExpr(
  mc_ALL,
  group_name = "ALL",
  group.by   = "ALL",
  assay      = "RNA",
  layer      = "data",
  wgcna_name = "Astrocyte031826"
)


mc_ALL <- TestSoftPowers(mc_ALL, networkType = "signed")
mc_ALL <- ConstructNetwork(
  mc_ALL,
  overwrite_tom = TRUE,
  tom_name = "ALL"
)

png("hdWGCNA_dendrogram_ALL_pc50_res08_k40_032526.png", width = 10, height = 7, onefile = TRUE)
PlotDendrogram(mc_ALL, main='ALL hdWGCNA Dendrogram')

dev.off()



mc_ALL <- ModuleEigengenes(mc_ALL)
ME <- GetMEs(mc_ALL)
# compute eigengene-based connectivity (kME):

mc_ALL <- ModuleConnectivity(
  mc_ALL,
  group.by = 'ALL', group_name = 'ALL'
)


modules <- GetModules(mc_ALL) %>% subset(module != 'grey')
# keeps only non-grey modules
modules <- modules[modules$module != "grey", ]

# gene sets：list(module -> genes)
module_genes <- split(modules$gene, modules$module)

# get hub genes
hub_df <- GetHubGenes(mc_ALL, n_hubs = 10)


mc_ALL <- ModuleExprScore(
  mc_ALL,
  n_genes = 25,
  method='UCell'
)

vip  <- vip_genes
vvip <- vvip_genes

mods <- GetModules(mc_ALL)
mods$module[is.na(mods$module)] <- "grey"

group1 <- mc_ALL@meta.data %>% subset(AD_status == 'AD') %>% rownames
group2 <- mc_ALL@meta.data %>% subset(AD_status == 'control') %>% rownames

head(group1)

# This is gene list level expression comparasion, not cell level
DMEs <- FindDMEs(
  mc_ALL,
  barcodes1 = group1,
  barcodes2 = group2,
  features = "MEs",
  test.use='wilcox',
  wgcna_name='Astrocyte031826'
)
DMEs

# TF analysis with hdWGCNA

ah <- AnnotationHub()
qry <- query(ah, c("EnsDb", "Homo sapiens", "105"))
mcols(qry)[, c("title", "species")]
edb105 <- qry[[1]]
edb105

JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

mc_ALL <- ModuleExprScore(
  mc_ALL,
  n_genes = 25,
  method='UCell'
)

mc_ALL <- MotifScan(
  mc_ALL,
  species_genome = "hg38",
  pfm = pfm_core,
  EnsDb = edb105
)

# get the motif df:
motif_df <- GetMotifs(mc_ALL)

# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(mc_ALL)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
mc_ALL <- SetWGCNAGenes(mc_ALL, genes_use)
mc_ALL <- SetDatExpr(
  mc_ALL,
  group_name = "ALL",       
  group.by   = "ALL",
  assay      = "RNA"
)
# define model params
model_params <- list(
    objective = 'reg:squarederror',
    max_depth = 1,
    eta = 0.1,
    nthread=32,
    alpha=0.5
)

# construct the TF network success！
mc_ALL <- ConstructTFNetwork(mc_ALL, model_params)

results <- GetTFNetwork(mc_ALL)
head(results)
write.csv(results,
          file = "031826_TF_results.csv",
          row.names = FALSE, quote = FALSE)

obj_ad  <- subset(mc_ALL, cells = colnames(mc_ALL)[mc_ALL$AD_status == "AD"])
obj_ctl <- subset(mc_ALL, cells = colnames(mc_ALL)[mc_ALL$AD_status == "control"])
p_AD <- DimPlot(
  obj_ad,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("AD")

p_control <- DimPlot(
  obj_ctl,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("control")

p_all <- DimPlot(
  mc_ALL,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("ALL")

ggsave("p_AD_dim.png", p_AD, width=8, height=8, dpi=600)
ggsave("p_control_dim.png", p_control, width=8, height=8, dpi=600)
ggsave("p_ALL_dim.png", p_all, width=8, height=8, dpi=600)

