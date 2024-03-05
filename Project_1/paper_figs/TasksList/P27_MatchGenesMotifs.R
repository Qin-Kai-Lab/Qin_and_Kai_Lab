library(JASPAR2022)
library(TFBSTools)

targDir <- 'OPC_ODC/Monocle3/86PC/'
genes = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv'))
genesFilt = unique(genes$gene_id[genes$q_value < 0.05])

opts <- list()
opts[["species"]] <- 10090
#opts[["tax_group"]] <- "vertebrates"

PFMatrixList <- getMatrixSet(JASPAR2022, opts)

PFMatrixList[[90]]@name
PFMatrixList[[90]]@tags$symbol
PFMatrixList[[90]]@tags$remap_tf_name

allMotifs = character()

for ( i in 1:length(PFMatrixList ) ) {
  curMotif  = PFMatrixList[[i]]@name
  curAlias = PFMatrixList[[i]]@tags$alias
  curMap = PFMatrixList[[i]]@tags$remap_tf_name
  if (length(curAlias) > 0 ) {
    curAlias = strsplit(curAlias, ",")[[1]]
  }
  allMotifs = c(allMotifs, curMotif, curAlias,  curMap)
  
}

allMotifs = unique(allMotifs)

signMotifs = genesFilt[genesFilt%in%allMotifs]


## more organisms

pfm1 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(10090), all_versions = FALSE))

pfm2 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(9606), all_versions = FALSE))

pfm3 <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", species = c(10116), all_versions = FALSE))

curNames = unique(c(names(pfm1), names(pfm2), names(pfm3)))

opts = list()
opts[["ID"]] <- curNames

PFMatrixList <- getMatrixSet(JASPAR2022, opts)

allMotifs = character()

for ( i in 1:length(PFMatrixList ) ) {
  curMotif  = PFMatrixList[[i]]@name
  curAlias = PFMatrixList[[i]]@tags$alias
  curMap = PFMatrixList[[i]]@tags$remap_tf_name
  if (length(curAlias) > 0 ) {
    curAlias = strsplit(curAlias, ",")[[1]]
  }
  allMotifs = c(allMotifs, curMotif, curAlias,  curMap)
  
}

allMotifs = unique(allMotifs)

overlapGenes = character()
for ( i in genesFilt ) {
  if (any(grepl(i, allMotifs, ignore.case = T))) {
    overlapGenes = c(overlapGenes, i)
  }
}

## all vertebrates
opts <- list()
#opts[["species"]] <- 10090
opts[["tax_group"]] <- "vertebrates"

PFMatrixList <- getMatrixSet(JASPAR2022, opts)

PFMatrixList[[90]]@name
PFMatrixList[[90]]@tags$symbol
PFMatrixList[[90]]@tags$remap_tf_name

allMotifs = character()

for ( i in 1:length(PFMatrixList ) ) {
  curMotif  = PFMatrixList[[i]]@name
  curAlias = PFMatrixList[[i]]@tags$alias
  curMap = PFMatrixList[[i]]@tags$remap_tf_name
  if (length(curAlias) > 0 ) {
    curAlias = strsplit(curAlias, ",")[[1]]
  }
  allMotifs = c(allMotifs, curMotif, curAlias,  curMap)
  
}

allMotifs = unique(allMotifs)

overlapGenes = character()
for ( i in genesFilt ) {
  if (any(grepl(i, allMotifs, ignore.case = T))) {
    overlapGenes = c(overlapGenes, i)
  }
}