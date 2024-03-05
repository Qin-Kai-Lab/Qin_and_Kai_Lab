library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)

setwd("~/Wang/output")

my_vector <- seq(from = 0, to = 5000, by = 500)
outFile = "Paper_figs/Fig5/Pvals_Gpr17_5KB_by500_Sum_2024-02-20.csv"

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

subObj = subset(atacFiltSub, newMonocClust == "OPC")

rm(atacFilt, atacFiltSub, opcOdc)
gc()

# get the TSS genomic range


curName = "Gpr17"

compFrags = function(curName, extUp, extDown, subObj) {
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
  
  ## 
  controlFile = "/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz"
  stressFile = "/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz"
  
  curChr = as.character(seqnames(transcripts))
  curStart = as.numeric(start(transcripts))
  curEnd = as.numeric(end(transcripts))
  
  getGeneRows = function(startPos, endPos, chr, curFile) {
    curRange = paste0(chr, ":", startPos, "-", endPos)
    cmd_str = paste0("tabix ", curFile, " ", curRange, " > curFragment.tsv")
    curList = system(cmd_str, intern = F)
    combDf = read.table("curFragment.tsv", F, sep = "\t")
    file.remove("curFragment.tsv")
    return(combDf)
  }
  
  
  curChrDfContr = getGeneRows(startPos=curStart, endPos=curEnd, chr=curChr, curFile = controlFile)
  curChrDfStress = getGeneRows(startPos=curStart, endPos=curEnd, chr=curChr, curFile = stressFile)
  
  # check that all regions overlap 
  contrRange = GRanges(seqnames = curChrDfContr$V1, ranges = IRanges(start = curChrDfContr$V2, end = curChrDfContr$V3), strand = strand("-"), gene_name = curName)
  overlapping_regions <- subsetByOverlaps(contrRange,transcripts)
  length(overlapping_regions) == nrow(curChrDfContr)
  
  #filter by cell ID
  
  curChrDfContr$Merged_CellName = paste0(curChrDfContr$V4, "_", 1)
  curChrDfStress$Merged_CellName = paste0(curChrDfStress$V4, "_", 2)
  
  filtContr = curChrDfContr[curChrDfContr$Merged_CellName%in%subObj@meta.data$Merged_CellName,]
  filtStress = curChrDfStress[curChrDfStress$Merged_CellName%in%subObj@meta.data$Merged_CellName,]
  
  # check using other cell ID
  contrCells = subObj@meta.data$gex_barcode[subObj@meta.data$dataset == "Control"]
  stressCells = subObj@meta.data$gex_barcode[subObj@meta.data$dataset == "Stress"]
  filtContr1 = curChrDfContr[curChrDfContr$V4%in%contrCells,]
  filtStress1 = curChrDfStress[curChrDfStress$V4%in%stressCells,]
  
  
  # sum fragments per cell
  contrSum = aggregate(V5 ~ V4, data = filtContr1, FUN = sum)
  stressSum = aggregate(V5 ~ V4, data = filtStress1, FUN = sum)
  
  # add 0s
  addMissingCells = function(df, cellList) {
    missCells = character()
    for ( curCell in cellList ) {
      if ( !curCell%in%df$V4 ) {
        missCells = c(missCells, curCell)
      }
    }
    missDf = data.frame(V4 = missCells, V5 = 0)
    combDf = rbind(df, missDf)
    return(combDf)
  }
  
  combContr = addMissingCells(df=contrSum, cellList=contrCells)
  nrow(combContr) == length(contrCells)
  combStress = addMissingCells(df=stressSum, cellList=stressCells)
  nrow(combStress) == length(stressCells)
  
  curTest = wilcox.test(combContr$V5, combStress$V5, alternative = "less")
  #curTest = t.test(combContr$V5, combStress$V5)
  pval = curTest$p.value
  return(pval)
}

combP = numeric()
combUp = numeric()
combDown = numeric()
for ( curUp in my_vector ) {
  for (curDown in my_vector ) {
    curPval = compFrags(curName=curName, extUp=curUp, extDown=curDown, subObj = subObj)
    #print(paste(curPval, curUp, curDown, sep = "__"))
    #if ( curPval < 0.1 ) {
      combP = c(combP, curPval)
      combUp = c(combUp, curUp)
      combDown = c(combDown, curDown)
    #}
  }
}

combPDf = data.frame(Pval = combP, Ext_Up = combUp, Ext_Down = combDown)
combPDf = combPDf[order(combPDf$Pval, decreasing = F),]

#write.csv(combPDf, outFile, row.names = F)