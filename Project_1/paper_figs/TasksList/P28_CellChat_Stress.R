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
cell.use = rownames(meta)[meta$group == "Stress" & meta$Annotations != "C-R"] # extract the cell names from disease data
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

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#saveRDS(cellchat,  "CellChat_RNA_Stress")

cellchat = readRDS("CellChat_RNA_Stress")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


targDir = "Paper_figs/TasksList/P28/Stress/"
dir.create(targDir, recursive = T, showWarnings = F)
# ggsave(paste0(targDir, "NetPlot_Number_Of_Interactions.png"), plot = p1,  height = 14, width = 14, units = 'in', dpi = 300)
# ggsave(paste0(targDir, "NetPlot_Interaction_strength.png"), plot = p2,  height = 14, width = 14, units = 'in', dpi = 300)

png(filename = paste0(targDir,"NetPlot_Number_Of_Interactions.png"), height = 6, width = 6, units = "in", res = 300)
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

png(filename = paste0(targDir,"NetPlot_Interaction_strength.png"), height = 6, width = 6, units = "in", res = 300)
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


###
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
for (i in 1:length(pathways.show.all)) {
  pathways.show = pathways.show.all[i]
  png(filename = paste0(targDir, pathways.show.all[i], "_2024-01-17.png"), height = 10, width = 10, units = "in", res = 300)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
}