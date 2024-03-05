library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

source('../programs/renameClusters.R')

editClusterNames = function(meta) {
  meta$Annotations = as.character(meta$Annotations)
  meta$Cluster = meta$Annotations
  for ( i in 1:nrow(meta) ) {
    if ( meta$Annotations[i] == "OPC" || meta$Annotations[i] == "ODC"  ) {
      meta$Cluster[i] = "OPC_ODC"
    }
  }
  return(meta)
}

data.input = RNA.combined.norm@assays$RNA$data 
meta = RNA.combined.norm@meta.data 
cell.use = rownames(meta)[meta$Annotations != "C-R"] # extract the cell names from disease data
print(length(cell.use))

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
meta = editClusterNames(meta)

# create object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cluster")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multisession", workers = 8) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat,  "CellChat_RNA_Combined")