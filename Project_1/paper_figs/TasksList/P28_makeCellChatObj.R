library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggplot2)
library(foreach)
library(doParallel)
library(future)

setwd("/home/flyhunter/Wang/output")

source('../programs/renameClusters.R')

opcOdc = readRDS('integ_OPC_ODC')

opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

opc_inter = rownames(opcOdc@meta.data)[opcOdc@meta.data$newMonocClust == "OPC" | opcOdc@meta.data$newMonocClust == "Intermideate"]
all(opc_inter%in%rownames(RNA.combined.norm@meta.data))

odc = rownames(opcOdc@meta.data)[opcOdc@meta.data$newMonocClust == "ODC" ]
all(odc%in%rownames(RNA.combined.norm@meta.data))
table(RNA.combined.norm@meta.data$Annotations)
table(opcOdc$newMonocClust)

editClusterNames = function(meta) {
  meta$Annotations = as.character(meta$Annotations)
  meta$Cluster = meta$Annotations
  for ( i in 1:nrow(meta) ) {
    if ( meta$Cell_ID[i] %in%opc_inter  ) {
      meta$Cluster[i] = "OPC_Inter"
    } else if (meta$Cell_ID[i]%in%odc) {
      meta$Cluster[i] = "ODC"
    } else {
      meta$Cluster[i] = meta$Annotations[i]
    }
  }
  return(meta)
}

data.input = RNA.combined.norm@assays$RNA$data 
meta = RNA.combined.norm@meta.data 
meta$Cell_ID = rownames(meta)
meta = editClusterNames(meta)

allGroups = unique(meta$group)
cores = 8 
curCluster  = makeForkCluster(cores)
registerDoParallel(curCluster)

rm(opcOdc, RNA.combined.norm)

gc()

getCellChatObj = function(data.input, i, meta) {
  cell.use = rownames(meta)[meta$group == i & meta$Annotations != "C-R"]
  # Prepare input data for CelChat analysis
  curData = data.input[, cell.use]
  curMeta  = meta[cell.use, ]
  #curMeta = editClusterNames(curMeta)
  
  # create object
  cellchat <- createCellChat(object = curData, meta = curMeta, group.by = "Cluster")
  print(i)
  print("ObjectCreated")
  cellchat <- addMeta(cellchat, meta = curMeta)
  cellchat <- setIdent(cellchat, ident.use = "Cluster") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  #showDatabaseCategory(CellChatDB)
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  cat("Database created \n")
  flush.console()
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 4) # do parallel
  print("Start identifyOverExpressedGenes(cellchat)")
  flush.console()
  cellchat <- identifyOverExpressedGenes(cellchat)
  print("Start identifyOverExpressedInteractions(cellchat)")
  flush.console()
  cellchat <- identifyOverExpressedInteractions(cellchat)
  print('Start computeCommunProb(cellchat, type = "triMean")')
  flush.console()
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  fileName = paste0("CellChat_RNA_", i, "_2.RDS")
  saveRDS(cellchat, fileName)
  
}

res = foreach(i = allGroups, .combine = 'c')%dopar% {
  getCellChatObj(data.input=data.input, i=i, meta=meta)
}

stopCluster(curCluster)
