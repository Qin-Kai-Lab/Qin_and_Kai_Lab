## Weights

maxWeight =  max(c(max(control@net$weight), max(stress@net$weight)))
# Output
png(paste0(targDir, "Network_Weights_Output_Combined.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,1), xpd=TRUE)
combGr = list()
combMat = list()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[c("OPC", "Intermideate", "ODC"), ] <- mat[c("OPC", "Intermideate", "ODC"), ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  combMat[[i]] = mat2
  combGr[[i]] = groupSize
}
avGr = (combGr[[1]] + combGr[[2]]) /2 
avMat = (combMat[[1]]  +  combMat[[2]]) /2
netVisual_circle(avMat, vertex.weight = avGr, weight.scale = T, edge.weight.max = maxWeight, color.use = custom_colors)
dev.off()

png(paste0(targDir, "Network_Weights_Input_Combined.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,1), xpd=TRUE)
combGr = list()
combMat = list()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$weight
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, c("OPC", "Intermideate", "ODC")] <- mat[, c("OPC", "Intermideate", "ODC")]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  combMat[[i]] = mat2
  combGr[[i]] = groupSize
}
avGr = (combGr[[1]] + combGr[[2]]) /2 
avMat = (combMat[[1]]  +  combMat[[2]]) /2
netVisual_circle(avMat, vertex.weight = avGr, weight.scale = T, edge.weight.max = maxWeight, color.use = custom_colors)
dev.off()

## Counts

maxcount = max(c(max(control@net$count), max(stress@net$count)))
png(paste0(targDir, "Network_Counts_Output_Combined.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,1), xpd=TRUE)
combGr = list()
combMat = list()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, count.scale = T, label.edge= F, edge.count.max = count.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[c("OPC", "Intermideate", "ODC"), ] <- mat[c("OPC", "Intermideate", "ODC"), ]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  combMat[[i]] = mat2
  combGr[[i]] = groupSize
}
avGr = (combGr[[1]] + combGr[[2]]) /2 
avMat = (combMat[[1]]  +  combMat[[2]]) /2
netVisual_circle(avMat, vertex.weight = avGr, weight.scale = T, edge.weight.max = maxcount, color.use = custom_colors)
dev.off()

png(paste0(targDir, "Network_Counts_Input_Combined.png"), width = 14, height = 14, units = "in", res = 300)
par(mfrow = c(1,1), xpd=TRUE)
combGr = list()
combMat = list()
for (i in 1:length(object.list)) {
  #netVisual_circle(object.list[[i]]@net$count, count.scale = T, label.edge= F, edge.count.max = count.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  title0 =  names(object.list)[i]
  mat <- object.list[[i]]@net$count
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[, c("OPC", "Intermideate", "ODC")] <- mat[, c("OPC", "Intermideate", "ODC")]
  groupSize <- as.numeric(table(object.list[[i]]@idents))
  combMat[[i]] = mat2
  combGr[[i]] = groupSize
}
avGr = (combGr[[1]] + combGr[[2]]) /2 
avMat = (combMat[[1]]  +  combMat[[2]]) /2
netVisual_circle(avMat, vertex.weight = avGr, weight.scale = T, edge.weight.max = maxcount, color.use = custom_colors)
dev.off()
