library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

targDir = "Paper_figs/TasksList/P28/Combined/All/"
dir.create(targDir, recursive = T, showWarnings = F)

stress = readRDS("CellChat_RNA_Stress")
control = readRDS("CellChat_RNA_Control")

stress <- netAnalysis_computeCentrality(stress, slot.name = "netP")
control <- netAnalysis_computeCentrality(control, slot.name = "netP")

object.list <- list(Control = control, Stress = stress)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),) + theme(text = element_text(size = 20))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") + theme(text = element_text(size = 20))
gg3 = gg1 + gg2

ggsave(paste0(targDir, "BarPlot.png"), plot = gg3, height = 10, width = 10, units = "in", dpi = 300)

targDir = "Paper_figs/TasksList/P28/Combined/OPC_ODC/All/"
dir.create(targDir, recursive = T, showWarnings = F)


# Weights
maxWeight = max(c(max(control@net$weight["OPC_ODC", ]), max(stress@net$weight["OPC_ODC", ])))

png(paste0(targDir, "Network_Weights.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2["OPC_ODC", ] <- mat["OPC_ODC", ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = maxWeight, title.name = title0)
}
dev.off()

#Counts

maxCount = max(c(max(control@net$count["OPC_ODC", ]), max(stress@net$count["OPC_ODC", ])))

png(paste0(targDir, "Network_Interraction_Counts.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2["OPC_ODC", ] <- mat["OPC_ODC", ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max =maxCount, title.name = title0)
}
dev.off()

# make bar plot

weightDf = data.frame()
for (i in 1:length(object.list)) {
    #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    title0 =  names(object.list)[i]
    mat <- object.list[[i]]@net$weight
    curVal <- sum(mat["OPC_ODC", ])
    curDf = data.frame(Group = title0, Weight = curVal)
    weightDf = rbind(weightDf, curDf)
}

countDf = data.frame()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  curVal <- sum(mat["OPC_ODC", ])
  curDf = data.frame(Group = title0, Count =curVal)
  countDf = rbind(countDf, curDf)
}

p1 = ggplot(weightDf, aes(x = Group, y = Weight, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p2 = ggplot(countDf, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p3 = p1 + p2 

ggsave(paste0(targDir, "BarPlot.png"), plot = p3, height = 10, width = 10, units = "in", dpi = 300)

#
levels(control@idents)
levels(stress@idents)

controlPathways <- control@netP$pathways
stressPathways = stress@netP$pathways
identical(controlPathways, stressPathways)

combPaths = controlPathways[controlPathways%in%stressPathways]



library(ComplexHeatmap)
targDir = "Paper_figs/TasksList/P28/Combined/All/"

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

png(paste0(targDir, "Outgoing_Signal.png"), width = 18, height = 18, units = "in", res = 300)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 18, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 18, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

png(paste0(targDir, "Incoming_Signal.png"), width = 18, height = 18, units = "in", res = 300)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 18, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 18, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

targDir = "Paper_figs/TasksList/P28/Combined/OPC_ODC/All/"
#netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8),  comparison = c(1, 2), angle.x = 45, font.size = 16)

png(paste0(targDir, "Communication_probabilites.png"), width = 18, height = 18, units = "in", res = 300)
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8),  comparison = c(1, 2), angle.x = 45, font.size = 16)
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg3 = gg1 + gg2

png(paste0(targDir, "Communication_probabilites_StrContr.png"), width = 18, height = 18, units = "in", res = 300)
print(gg3)
dev.off()


targDir = "Paper_figs/TasksList/P28/Combined/All/Pathwayas/"
dir.create(targDir)
for ( curPathway in combPaths ) {
  png(paste0(targDir, curPathway, ".png"), width = 10, height = 10, units = "in", res = 300)
  pathways.show <- curPathway
  #weight.max <- max(c( getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show), getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathways.show) )) # control the edge weights across different datasets
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  }
  dev.off()
}

targDir = "Paper_figs/TasksList/P28/Combined/OPC_ODC/Pathwayas/"
dir.create(targDir)
for ( curPathway in combPaths ) {
  png(paste0(targDir, curPathway, ".png"), width = 10, height = 10, units = "in", res = 300)
  pathways.show <- curPathway
  #weight.max <- max(c( getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show), getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathways.show) )) # control the edge weights across different datasets
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                        signaling.name = paste(pathways.show, names(object.list)[i]), sources.use = 7, targets.use = c(1:6,8))
  }
  dev.off()
}
