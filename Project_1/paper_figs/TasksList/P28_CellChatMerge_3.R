library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ggplot2)

targDir = "Paper_figs/TasksList/P28/Merged_3/All/"
dir.create(targDir, recursive = T, showWarnings = F)

stress = readRDS("CellChat_RNA_Stress_3.RDS")
control = readRDS("CellChat_RNA_Control_3.RDS")

#stress <- netAnalysis_computeCentrality(stress, slot.name = "netP")
#control <- netAnalysis_computeCentrality(control, slot.name = "netP")

object.list <- list(Control = control, Stress = stress)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),) + theme(text = element_text(size = 20))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") + theme(text = element_text(size = 20))
gg3 = gg1 + gg2

ggsave(paste0(targDir, "BarPlot.png"), plot = gg3, height = 10, width = 10, units = "in", dpi = 300)

targDir = "Paper_figs/TasksList/P28/Merged_3/OPC_ODC/"
dir.create(targDir, recursive = T, showWarnings = F)

levels(control@idents)

custom_colors <- c("chocolate4", "blue", "forestgreen", "purple", 'magenta1',  "orange",  "aquamarine3", "red", "black", "gold")

## Weights

maxWeight =  max(c(max(control@net$weight), max(stress@net$weight)))
# Output
png(paste0(targDir, "Network_Weights_Output.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[c("OPC", "Intermideate", "ODC"), ] <- mat[c("OPC", "Intermideate", "ODC"), ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = maxWeight, title.name = title0, color.use = custom_colors)
}
dev.off()

png(paste0(targDir, "Network_Weights_Input.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, c("OPC", "Intermideate", "ODC")] <- mat[, c("OPC", "Intermideate", "ODC")]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = maxWeight, title.name = title0, color.use = custom_colors)
}
dev.off()

## Counts

maxcount = max(c(max(control@net$count), max(stress@net$count)))
png(paste0(targDir, "Network_Counts_Output.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, count.scale = T, label.edge= F, edge.count.max = count.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[c("OPC", "Intermideate", "ODC"), ] <- mat[c("OPC", "Intermideate", "ODC"), ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = maxcount, title.name = title0, color.use = custom_colors)
}
dev.off()

png(paste0(targDir, "Network_Counts_Input.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, count.scale = T, label.edge= F, edge.count.max = count.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, c("OPC", "Intermideate", "ODC")] <- mat[, c("OPC", "Intermideate", "ODC")]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = maxcount, title.name = title0, color.use = custom_colors)
}
dev.off()

## make bar plot

# Ouput

weightDf = data.frame()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  curVal <- c(sum(mat["OPC", ]), sum(mat["Intermideate", ]), sum(mat["ODC", ]))
  curDf = data.frame(Group = title0, Weight = curVal)
  curDf$Cluster = c("OPC", "Intermideate", "ODC")
  weightDf = rbind(weightDf, curDf)
}

countDf = data.frame()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  curVal <-  c(sum(mat["OPC", ]), sum(mat["Intermideate", ]), sum(mat["ODC", ]))
  curDf = data.frame(Group = title0, Count =curVal)
  curDf$Cluster = c("OPC", "Intermideate", "ODC")
  countDf = rbind(countDf, curDf)
}

weightDf$Cluster = factor(weightDf$Cluster, levels = c("OPC", "Intermideate", "ODC"))
countDf$Cluster = factor(countDf$Cluster, levels = c("OPC", "Intermideate", "ODC"))

p1 = ggplot(weightDf, aes(x = Group, y = Weight, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p1 = p1 + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

p2 = ggplot(countDf, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p2 = p2 + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

p3 = p1 + p2 

ggsave(paste0(targDir, "BarPlot_Output.png"), plot = p3, height = 10, width = 16, units = "in", dpi = 300)

# input

weightDf = data.frame()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  curVal <- c(sum(mat[, "OPC"]), sum(mat[, "Intermideate"]),  sum(mat[, "ODC"]))
  curDf = data.frame(Group = title0, Weight = curVal)
  curDf$Cluster = c("OPC", "Intermideate", "ODC")
  weightDf = rbind(weightDf, curDf)
}

countDf = data.frame()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  curVal <- c(sum(mat[, "OPC"]), sum(mat[, "Intermideate"]),  sum(mat[, "ODC"]))
  curDf = data.frame(Group = title0, Count =curVal)
  curDf$Cluster = c("OPC", "Intermideate", "ODC")
  countDf = rbind(countDf, curDf)
}

weightDf$Cluster = factor(weightDf$Cluster, levels = c("OPC", "Intermideate", "ODC"))
countDf$Cluster = factor(countDf$Cluster, levels = c("OPC", "Intermideate", "ODC"))

p1 = ggplot(weightDf, aes(x = Group, y = Weight, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p1 = p1 + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

p2 = ggplot(countDf, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

p2 = p2 + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

p3 = p1 + p2 

ggsave(paste0(targDir, "BarPlot_Input.png"), plot = p3, height = 10, width = 16, units = "in", dpi = 300)

#
levels(control@idents)
levels(stress@idents)

controlPathways <- control@netP$pathways
stressPathways = stress@netP$pathways
identical(controlPathways, stressPathways)

combPaths = controlPathways[controlPathways%in%stressPathways]



library(ComplexHeatmap)
targDir = "Paper_figs/TasksList/P28/Merged_3/All/"

#pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

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

targDir = "Paper_figs/TasksList/P28/Merged_3/OPC_ODC/"
#dir.create(targDir, recursive = T, showWarnings = F)
#netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6,8),  comparison = c(1, 2), angle.x = 45, font.size = 16)

levels(control@idents)


png(paste0(targDir, "Communication_Probabilites_Output.png"), width = 24, height = 28, units = "in", res = 300)
netVisual_bubble(cellchat, sources.use = c(6,8:9), targets.use = c(1:5,7,10),  comparison = c(1, 2), angle.x = 45, font.size = 16)
dev.off()

png(paste0(targDir, "Communication_Probabilites_Input.png"), width = 24, height = 38, units = "in", res = 300)
netVisual_bubble(cellchat, sources.use = c(1:5,7,10), targets.use = c(6,8:9),  comparison = c(1, 2), angle.x = 45, font.size = 16)
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = c(6,8:9), targets.use = c(1:5,7,10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(6,8:9), targets.use = c(1:5,7,10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg3 = gg1 + gg2

png(paste0(targDir, "Communication_Probabilites_StrContr_Output.png"), width = 34, height = 32, units = "in", res = 300)
print(gg3)
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = c(1:5,7,10), targets.use = c(6,8:9),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:5,7,10), targets.use = c(6,8:9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Stress", angle.x = 45, remove.isolate = T, font.size = 16)
#> Comparing communications on a merged object
gg3 = gg1 + gg2

png(paste0(targDir, "Communication_Probabilites_StrContr_Input.png"), width = 34, height = 36, units = "in", res = 300)
print(gg3)
dev.off()


targDir = "Paper_figs/TasksList/P28/Merged_3/All/Pathwayas/"
dir.create(targDir)
for ( curPathway in combPaths ) {
  png(paste0(targDir, curPathway, ".png"), width = 10, height = 10, units = "in", res = 300)
  pathways.show <- curPathway
  #weight.max <- max(c( getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show), getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathways.show) )) # control the edge weights across different datasets
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                        signaling.name = paste(pathways.show, names(object.list)[i]), color.use = custom_colors)
  }
  dev.off()
}

targDir = "Paper_figs/TasksList/P28/Merged_3/OPC_ODC/Pathwayas_Output/"
dir.create(targDir)
for ( curPathway in combPaths ) {
  png(paste0(targDir, curPathway, ".png"), width = 10, height = 10, units = "in", res = 300)
  pathways.show <- curPathway
  #weight.max <- max(c( getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show), getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathways.show) )) # control the edge weights across different datasets
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                        signaling.name = paste(pathways.show, names(object.list)[i]), sources.use = c(6,8:9), targets.use = c(1:5,7,10), color.use = custom_colors)
  }
  dev.off()
}

targDir = "Paper_figs/TasksList/P28/Merged_3/OPC_ODC/Pathwayas_Input/"
dir.create(targDir)
for ( curPathway in combPaths ) {
  png(paste0(targDir, curPathway, ".png"), width = 10, height = 10, units = "in", res = 300)
  pathways.show <- curPathway
  #weight.max <- max(c( getMaxWeight(object.list[1], slot.name = c("netP"), attribute = pathways.show), getMaxWeight(object.list[2], slot.name = c("netP"), attribute = pathways.show) )) # control the edge weights across different datasets
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, 
                        signaling.name = paste(pathways.show, names(object.list)[i]), sources.use = c(1:5,7,10), targets.use = c(6,8:9), color.use = custom_colors)
  }
  dev.off()
}
