library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(MAST)
library(ggplot2)
library(data.table)

setwd("/home/flyhunter/Wang/output")

opcObj = readRDS("OPC_CombPeaks_macs2_2.rds")

curName = 'chr18-31788658-31789364'

# diff expression
compOpen = function(atacFiltSub, curName, curCluster) {
  if (curCluster == "all") {
    subObj = atacFiltSub
  } else {
    subObj = subset(atacFiltSub, newMonocClust == curCluster)
  }
  
  transcripts <- StringToGRanges(curName)
  
  frags <- Fragments(object = subObj[["combined_peaks"]])
  cells <- colnames(x = subObj[["combined_peaks"]])
  counts <- FeatureMatrix(fragments = frags, features = transcripts, 
                          process_n = 2000, cells = cells, verbose = T)

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

filtContr = curChrDfContr[curChrDfContr$Merged_CellName%in%subObj@meta.data$Merge_CellID,]
filtStress = curChrDfStress[curChrDfStress$Merged_CellName%in%subObj@meta.data$Merge_CellID,]

# check using other cell ID
contrCells = subObj@meta.data$Orig_ID[subObj@meta.data$Group == "Control"]
stressCells = subObj@meta.data$Orig_ID[subObj@meta.data$Group == "Stress"]
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
curTest