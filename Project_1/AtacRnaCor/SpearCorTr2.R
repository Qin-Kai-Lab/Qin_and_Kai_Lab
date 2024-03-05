library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(MAST)
library(ggplot2)
library(rtracklayer)


setwd("~/Wang/output")

source('../programs/renameClusters.R')

allMarkers = read.csv("allMarkersRenameClust0.2Res_2022-10-21.csv")

curDate<-Sys.Date()

atacInt<-readRDS('atacIntegrated_macs2_2')
DefaultAssay(atacInt)

# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atacInt) <- annotations

#combine rna and atac
RNA.combined.norm$Merged_CellName<-colnames(RNA.combined.norm)
atacInt$Merged_CellName<-colnames(atacInt)
rnaFilt<-subset(RNA.combined.norm, subset = Merged_CellName %in% atacInt$Merged_CellName )
atacFilt = subset(atacInt, subset = Merged_CellName %in% rnaFilt$Merged_CellName )

identical(colnames(rnaFilt), colnames(atacFilt))
identical(rownames(rnaFilt@meta.data), rownames(atacFilt@meta.data))
identical(rnaFilt$Merged_CellName, atacFilt$Merged_CellName)

atacFilt[['RNA']]<-rnaFilt@assays$RNA
DefaultAssay(atacInt)<-'Combined_peaks'
atacFilt$Annotations = rnaFilt$Annotations

rm(RNA.combined.norm, rnaFilt, atacInt)
gc()

Idents(atacFilt) = atacFilt$Annotations

#
cond<-'Control_Stress'
targetDir<-paste0('Atac_Rna_coverPlots/macs2_Integ/', cond, '/')
dir.create(targetDir, recursive = T)

annotGenes<-data.frame(Annotation(atacFilt))

atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)

table(atacFilt$Annotations)


# get cell names 
# workw with atacMarix
subset_seurat_object <- subset(atacFilt, Annotations == "OPC")
targCellNames = subset_seurat_object$Merged_CellName
# get gene expression
gene1_expr = subset_seurat_object@assays$RNA@counts["Sox6" ,]

# find gene coordinates
extraRange = 20000

subset_gene <- annotations[annotations$gene_name == "Sox6"]
gene1 = data.frame(subset_gene@ranges)
startPos = gene1$start[1] - extraRange
endPos = gene1$end[nrow(gene1)] + extraRange
#endPos = gene1$end[nrow(gene1)]
chr = as.character(subset_gene@seqnames@values)

# get fragmetns
findOverlappingRanges = function(refStart, refEnd, rangesDf) {
  selRanges = data.frame()
  for ( i in 1:nrow(rangesDf) ) {
    testStart = rangesDf$V2[i]
    testEnd = rangesDf$V3[i]
    if ((refStart <= testEnd) && (testStart <= refEnd)) {
      curRow = rangesDf[i,]
      selRanges = rbind(selRanges,  curRow)
    }
  }
  return(selRanges)
}

controlFile = "/home/flyhunter/Wang/data/control/atac_fragments.tsv.gz"
stressFile = "/home/flyhunter/Wang/data/stress/atac_fragments.tsv.gz"

# use tabix to filter ranges
# getGeneRows = function(startPos, endPos, chr, curFile) {
#   curRange = paste0(chr, ":", startPos, "-", endPos)
#   cmd_str = paste0("tabix ", curFile, " ", curRange)
#   curList = system(cmd_str, intern = T)
#   combDf = data.frame()
#   for ( i in 1:length(curList) ) {
#     splitLine = strsplit(curList[i], "\t")[[1]]
#     curDf = data.frame(t(splitLine), stringsAsFactors = F)
#     combDf = rbind(combDf,curDf)
#   } 
#   return(combDf)
# }

# writing on disk with tabix is faster than combining lines in R memory
getGeneRows = function(startPos, endPos, chr, curFile) {
  curRange = paste0(chr, ":", startPos, "-", endPos)
  cmd_str = paste0("tabix ", curFile, " ", curRange, " > curFragment.tsv")
  curList = system(cmd_str, intern = F)
  combDf = read.table("curFragment.tsv", F, sep = "\t")
  file.remove("curFragment.tsv")
  return(combDf)
}

#42738
curChrDfContr = getGeneRows(startPos=startPos, endPos=endPos, chr=chr, curFile = controlFile)

startTime = Sys.time()
curChrDfStress = getGeneRows(startPos=startPos, endPos=endPos, chr=chr, curFile = stressFile)
endTime = Sys.time()
print(endTime - startTime)



#combineRows = function(linesList, refStart, refEnd) {
  #   combDf = data.frame()
  #   for ( i in 1:length(linesList) ) {
  #     if (!startsWith(linesList[i], "#")) {
  #       splitLine = strsplit(linesList[i], "\t")[[1]]
  #       testStart = as.numeric(splitLine[2])
  #       testEnd = as.numeric(splitLine[3])
  #       if ((refStart <= testEnd) && (testStart <= refEnd)) {
  #         curDf = data.frame(t(splitLine), stringsAsFactors = F)
  #         combDf = rbind(combDf,curDf)
  #       }
  #     }
  #   }
  #   return(combDf)
  # }

# cmd_str = paste0("grep ", chr, " ", controlFile)
# cmd_str_stress = paste0("grep ", chr, " ", stressFile)
# 
# # get rows for current chromosome
# linesListControl = system(cmd_str, intern = T)
# linesListStress = system(cmd_str_stress, intern = T)
# 
# # combine rows, this one can be paralel for all genes in the list
# combineRows = function(linesList, refStart, refEnd) {
#   combDf = data.frame()
#   for ( i in 1:length(linesList) ) {
#     if (!startsWith(linesList[i], "#")) {
#       splitLine = strsplit(linesList[i], "\t")[[1]]
#       testStart = as.numeric(splitLine[2])
#       testEnd = as.numeric(splitLine[3])
#       if ((refStart <= testEnd) && (testStart <= refEnd)) {
#         curDf = data.frame(t(splitLine), stringsAsFactors = F)
#         combDf = rbind(combDf,curDf)
#       }
#     }
#   }
#   return(combDf)
# }
# startTime = Sys.time()
# curChrDfContr = combineRows(linesList=linesListControl, refStart = startPos, refEnd = endPos)
# 
# curChrDfStress = combineRows(linesList=linesListStress, refStart = startPos, refEnd = endPos)
# endTime = Sys.time()
# print(endTime - startTime)

#parallel for reads
# combineRows = function(i) {
#   if (!startsWith(linesList[i], "#")) {
#     splitLine = strsplit(i, "\t")[[1]]
#     curDf = data.frame(t(splitLine), stringsAsFactors = F)
#     return(curDf)
#   }
# }
# runPar = function(linesList) {
#   useCores = 14
#   cl <- makeCluster(useCores, type = "FORK")
#   registerDoParallel(cl)
# 
#   combDat = foreach(i = linesList, .combine  = 'rbind') %dopar% {
#     combineRows(i = i)
#   }
#   parallel::stopCluster(cl = cl)
#   return(combDat)
# }
# 
# startTime =  Sys.time()
# curChrDfContrFull = runPar(linesList=linesListStress)
# endTime = Sys.time()
# print(endTime - startTime)
# 


sumFrags = function(contrDf, stressDf, targCellNames) {
  curNames = c("Chr", "Start", "End", "Barcode", "Fragments")
  colnames(contrDf) = curNames
  colnames(stressDf) = curNames
  contrDf$Barcode = paste0(contrDf$Barcode, "_1")
  contrDfSel = contrDf[contrDf$Barcode %in% targCellNames,]
  stressDf$Barcode = paste0(stressDf$Barcode, "_2")
  stressDfSel = stressDf[stressDf$Barcode %in% targCellNames,]
  combDat = rbind(contrDfSel, stressDfSel)
  combDat$Fragments = as.numeric(combDat$Fragments)
  fragmentsSum = aggregate(Fragments~Barcode, data = combDat, FUN = sum)
  return(fragmentsSum)
}

summedFragments = sumFrags(contrDf=curChrDfContr, stressDf=curChrDfStress, targCellNames=targCellNames)

addZeoroVals = function(df, targCellNames) {
  selCells = data.frame(targCellNames[!targCellNames%in%df$Barcode])
  colnames(selCells) = "Barcode"
  selCells$Fragments = 0
  combDf = rbind(df, selCells)
  return(combDf)
}

completeAtac = addZeoroVals(df=summedFragments, targCellNames=targCellNames)
# finish RNA seq data
completeRna = data.frame(cbind(names(gene1_expr),gene1_expr))
colnames(completeRna) = c('Barcode', "RNA_counts")

combRnaAtac = plyr::join(completeRna, completeAtac, by = "Barcode", type = "left", match = "all")
combRnaAtac$RNA_counts = as.numeric(combRnaAtac$RNA_counts)
combRnaAtac$Fragments = as.numeric(combRnaAtac$Fragments)

cor.test(combRnaAtac$RNA_counts, combRnaAtac$Fragments,  method = "spearman")

cor.test(combRnaAtac$RNA_counts, combRnaAtac$Fragments,  method = "pearson")



# combRnaAtacSel = combRnaAtac[combRnaAtac$RNA_counts > 0 | combRnaAtac$Fragments > 0 ,]
# cor.test(combRnaAtacSel$RNA_counts, combRnaAtacSel$Fragments,  method = "spearman")
# 
# library(MASS)
# summary(glm.nb(combRnaAtacSel$RNA_counts~combRnaAtacSel$Fragments))
summary(glm.nb(combRnaAtac$RNA_counts~combRnaAtac$Fragments))
