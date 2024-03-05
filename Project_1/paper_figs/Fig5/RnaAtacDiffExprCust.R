library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(MAST)
library(ggplot2)
library(data.table)

setwd("/home/flyhunter/Wang/output")

opcObj = readRDS("OPC_CombPeaks_macs2_2.rds")

opcObj <- RegionStats(opcObj, genome = BSgenome.Mmusculus.UCSC.mm10)

custom_colors <- c("red", "blue", "green", "orange", "violet", "black", "brown")

CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

opcObj = LinkPeaks(
  opcObj,
  peak.assay = "combined_peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  distance = 5e+05,
  #distance = 5000,
  min.cells = 3,
  genes.use = "Gpr17",
  n_sample = 200,
  pvalue_cutoff = 1,
  score_cutoff = 0,
  gene.id = FALSE,
  verbose = TRUE)


p2g = data.frame(Links(opcObj))
endPost = as.numeric(as.character(gsub("^.*\\-","", p2g$peak )))
endPos1 = 31949636- endPost
p2g$Distance_TSS = endPos1 

#write.csv(p2g, "~/Wang/output/Paper_figs/Fig5/OPC_Peaks_Gpr17_Cor_500KB_ALL_2024-02-20.csv", row.names = F)

p3 = CoveragePlot(
  object = opcObj ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p3)

p4 = CoveragePlot(
  object = opcObj ,
  region = "Gpr17",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000, group.by = "Group"
) & theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

print(p4)

p5 = CoveragePlot(
  object = opcObj ,
  region = "chr18-31918192-31918651",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  50,
  extend.downstream = 50, group.by = "Group", window = 100, annotation = F, assay.scale = "common", show.bulk = T) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

p5

curPeak = opcObj@assays$combined_peaks@counts["chr18-32162139-32163924",]
q1 = curPeak[curPeak > 0]

contrOpen = opcContr@assays$combined_peaks@counts["chr18-31918192-31918651",]
stressOpen = opcStress@assays$combined_peaks@counts["chr18-31918192-31918651",]
mean(contrOpen)
mean(stressOpen)
sum(contrOpen)
sum(stressOpen)
hist(contrOpen)
hist(stressOpen)

# split groups
# Control
opcContr = subset(opcObj, Group == "Control")
opcContr = LinkPeaks(
  opcContr,
  peak.assay = "combined_peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  distance = 5e+05,
  #distance = 5000,
  min.cells = 3,
  genes.use = "Gpr17",
  n_sample = 200,
  pvalue_cutoff = 1,
  score_cutoff = 0,
  gene.id = FALSE,
  verbose = TRUE)


p2g_contr = data.frame(Links(opcContr))
endPost = as.numeric(as.character(gsub("^.*\\-","", p2g_contr$peak )))
endPos1 = 31949636- endPost
p2g_contr$Distance_TSS = endPos1 

# Stress
opcStress = subset(opcObj, Group == "Stress")
opcStress = LinkPeaks(
  opcStress,
  peak.assay = "combined_peaks",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  distance = 5e+05,
  #distance = 5000,
  min.cells = 3,
  genes.use = "Gpr17",
  n_sample = 200,
  pvalue_cutoff = 1,
  score_cutoff = 0,
  gene.id = FALSE,
  verbose = TRUE)


p2g_Stress = data.frame(Links(opcStress))
endPost = as.numeric(as.character(gsub("^.*\\-","", p2g_Stress$peak )))
endPos1 = 31949636- endPost
p2g_Stress$Distance_TSS = endPos1 

p2g_contr = p2g_contr[(p2g_contr$pvalue < 0.1 & p2g_contr$score > 0),]
p2g_Stress = p2g_Stress[(p2g_Stress$pvalue < 0.1 & p2g_Stress$score > 0),]

stress_peak_uniq = p2g_Stress[!p2g_Stress$peak%in%p2g_contr$peak,]
stress_peak_uniq$peak
write.csv(stress_peak_uniq, "~/Wang/output/Paper_figs/Fig5/OPC_Peaks_LinkExpr_Gpr17_500KB_StressUniq_2024-02-22.csv", row.names = F)

p5 = CoveragePlot(
  object = opcObj ,
  region = "chr18-32391152-32391674",
  features = "Gpr17",
  expression.assay = "RNA",
  extend.upstream =  50,
  extend.downstream = 50, group.by = "Group", window = 100, annotation = F, assay.scale = "common", show.bulk = T) & 
  theme(text = element_text(size = 22)) & scale_fill_manual(values = custom_colors) 

p5

contrOpen = opcContr@assays$combined_peaks@counts["chr18-32391152-32391674",]
stressOpen = opcStress@assays$combined_peaks@counts["chr18-32391152-32391674",]
mean(contrOpen)
mean(stressOpen)
sum(contrOpen)
sum(stressOpen)

# diff expression
compOpen = function(atacFiltSub, curName, curCluster, extUp, extDown) {
  if (curCluster == "all") {
    subObj = atacFiltSub
  } else {
    subObj = subset(atacFiltSub, newMonocClust == curCluster)
  }
  annotation <- Annotation(object = subObj[["combined_peaks"]])
  curAnnot = annotation[annotation$gene_name == curName,]
  transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
  curStrand = as.character(transcripts@strand)
  if ( curStrand == "-") {
    gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = end(transcripts), end = end(transcripts)), strand = strand("-"), gene_name = curName)
  } else {
    gr <- GRanges(seqnames = seqnames(transcripts), ranges = IRanges(start = start(transcripts), end = start(transcripts)), strand = strand("+"), gene_name = curName)
  }
  
  transcripts <- Extend(x = gr, upstream = extUp, 
                        downstream = extDown)
  
  frags <- Fragments(object = subObj[["combined_peaks"]])
  cells <- colnames(x = subObj[["combined_peaks"]])
  counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                          process_n = 2000, cells = cells, verbose = T)
  gene.key <- transcripts$gene_name
  names(x = gene.key) <- GRangesToString(grange = transcripts)
  rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
  counts <- counts[rownames(x = counts) != "", ]
  
  contrGr = subObj@meta.data$Merge_CellID[subObj@meta.data$Group == "Control"]
  strGr = subObj@meta.data$Merge_CellID[subObj@meta.data$Group == "Stress"]
  
  contrDat = counts[names(counts)%in%contrGr]
  strDat = counts[names(counts)%in%strGr]
  
  curTest = wilcox.test(strDat, contrDat, alternative = "greater")
  #curTest = t.test(strDat, contrDat, alternative = "greater")
  #curTest = t.test(strDat, contrDat)
  print(curTest)
  #print(mean(strDat))
  #print(mean(contrDat))
  curList = list(strDat, contrDat)
  return(curTest)
}

atacExpr = compOpen(atacFiltSub=opcObj, curName="Gpr17", curCluster="all", extUp=2500, extDown=500)

rnaExpr = FindMarkers(opcObj, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "MAST", assay = "RNA", group.by = "Group")
#curExpr1 = FindMarkers(curClust, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
rnaExpr$Gene = rownames(rnaExpr)
Gpr17 = rnaExpr[rnaExpr$Gene == "Gpr17",]

rnaExpr2 = FindMarkers(opcObj, ident.1 = "Stress", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox", assay = "RNA", group.by = "Group")
rnaExpr2$Gene = rownames(rnaExpr2)
Gpr17_2 = rnaExpr2[rnaExpr2$Gene == "Gpr17",]

# diff expression Str Contr
da_peaks <- FindMarkers(
  object = opcObj,
  ident.1 = "Stress",
  ident.2 = "Control",
  test.use = 'LR',
  latent.vars = 'nCount_combined_peaks',
  logfc.threshold = 0, only.pos = F, group.by = "Group")

da_peaks$Peak = rownames(da_peaks)

genes = ClosestFeature(opcObj, regions = da_peaks$Peak)

combDat = cbind(da_peaks, genes)
combDat$fdr_p = p.adjust(combDat$p_val, method = "fdr")
combDat$bonferroni_p = p.adjust(combDat$p_val, method = "bonferroni")

targDir = './Paper_figs/Fig5/'

write.csv(combDat, 'Paper_figs/Fig5/DiffExprStrContr_OpcPeaks_2024-02-20.csv', row.names = F)

gpr = combDat[combDat$gene_name == "Gpr17",]

annotation <- Annotation(object = opcObj[["combined_peaks"]])
curAnnot = annotation[annotation$gene_name == "Gpr17",]
transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
transcripts <- Extend(x = transcripts, upstream = 5000, 
                      downstream = 5000)

gr <- GRanges(
  seqnames = sapply(strsplit(combDat$Peak, "-"), `[`, 1),
  ranges = IRanges(
    start = as.numeric(sapply(strsplit(combDat$Peak, "-"), `[`, 2)),
    end = as.numeric(sapply(strsplit(combDat$Peak, "-"), `[`, 3))
  )
)


#
# Extract overlapping regions from the original GRanges object
overlapping_regions <- subsetByOverlaps(gr,transcripts)

gprOverlapPeaks = c("chr18-31943368-31943733", "chr18-31950032-31950341", "chr18-31938975-31939403")

gprOverlap = combDat[combDat$Peak%in%gprOverlapPeaks,]