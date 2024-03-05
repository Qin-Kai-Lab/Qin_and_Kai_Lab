library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)

editClusterNames = function(meta) {
  meta$newMonocClust = as.character(meta$newMonocClust)
  meta$Cluster = meta$newMonocClust
  for ( i in 1:nrow(meta) ) {
    if ( meta$Cell_ID[i] %in%opc_inter  ) {
      meta$Cluster[i] = "OPC_Inter"
    } else if (meta$Cell_ID[i]%in%odc) {
      meta$Cluster[i] = "ODC"
    } else {
      meta$Cluster[i] = meta$newMonocClust[i]
    }
  }
  return(meta)
}

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

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

opcOdc =  readRDS('integ_OPC_ODC')
opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

opcOdc$Cell_ID = rownames(opcOdc@meta.data)

atacFiltSub = subset(atacFilt, Merged_CellName%in%rownames(opcOdc@meta.data))
opcOdcSub = subset(opcOdc, Cell_ID%in%rownames(atacFiltSub@meta.data))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

opc_inter = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "OPC" | opcOdcSub@meta.data$newMonocClust == "Intermideate"]
odc = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "ODC" ]

opcOdcSub@meta.data = editClusterNames(opcOdcSub@meta.data)

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

atacFiltSub$newMonocClust = opcOdcSub$newMonocClust
atacFiltSub$Cluster = opcOdcSub$Cluster

#
compOpen = function(atacFiltSub, curName, curCluster, extUp, extDown) {
  if (curCluster == "all") {
    subObj = atacFiltSub
  } else {
    subObj = subset(atacFiltSub, newMonocClust == curCluster)
  }
  annotation <- Annotation(object = subObj[["Combined_peaks"]])
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
  
  frags <- Fragments(object = subObj[["Combined_peaks"]])
  cells <- colnames(x = subObj[["Combined_peaks"]])
  counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                          process_n = 2000, cells = cells, verbose = T)
  gene.key <- transcripts$gene_name
  names(x = gene.key) <- GRangesToString(grange = transcripts)
  rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
  counts <- counts[rownames(x = counts) != "", ]
  
  contrGr = subObj@meta.data$Merged_CellName[subObj@meta.data$dataset == "Control"]
  strGr = subObj@meta.data$Merged_CellName[subObj@meta.data$dataset == "Stress"]
  
  contrDat = counts[names(counts)%in%contrGr]
  strDat = counts[names(counts)%in%strGr]
  
  curTest = wilcox.test(strDat, contrDat)
  #curTest = t.test(strDat, contrDat)
  #print(curTest)
  #print(mean(strDat))
  #print(mean(contrDat))
  #curList = list(strDat, contrDat)
  pval = curTest$p.value
  return(pval)
}

# kb1 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 1000, extDown = 1000)
# kb_2_1 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 2000, extDown = 1000)
# kb_2_2 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 2000, extDown = 2000)
# kb_3_0.5 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 3000, extDown = 500)
# b_200 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 200, extDown = 200)
# kb_2_0.5 = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = 2000, extDown = 500)
# 
# compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "Intermideate", extUp = 2000, extDown = 1000)
# compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "ODC", extUp = 2000, extDown = 1000)
# compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "all", extUp = 2000, extDown = 1000)

my_vector <- seq(from = 0, to = 5000, by = 1000)

combP = numeric()
combUp = numeric()
combDown = numeric()
for ( curUp in my_vector ) {
  for (curDown in my_vector ) {
    curPval = compOpen(atacFiltSub = atacFiltSub, curName = "Gpr17", curCluster = "OPC", extUp = curUp, extDown = curDown)
    #print(paste(curPval, curUp, curDown, sep = "__"))
    #if ( curPval < 0.1 ) {
    combP = c(combP, curPval)
    combUp = c(combUp, curUp)
    combDown = c(combDown, curDown)
    #}
  }
}

combPDf_3 = data.frame(Pval = combP, Ext_Up = combUp, Ext_Down = combDown)



