library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ggplot2)

targDir = "Paper_figs/TasksList/P28/Merged_2/All/"
dir.create(targDir, recursive = T, showWarnings = F)

stress = readRDS("CellChat_RNA_Stress_2.RDS")
control = readRDS("CellChat_RNA_Control_2.RDS")

#stress <- netAnalysis_computeCentrality(stress, slot.name = "netP")
#control <- netAnalysis_computeCentrality(control, slot.name = "netP")

object.list <- list(Control = control, Stress = stress)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Stress"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Stress",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Control",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(7,8), targets.use = c(1:6,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(7,8), targets.use = c(1:6,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


write.csv(net.up , paste0(targDir, "DiffExpr_StressUp.csv"), row.names = F)
write.csv(net.down , paste0(targDir, "DiffExpr_StressDown.csv"), row.names = F)

# Chord diagram

# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_chord_gene(object.list[[2]], sources.use = c(7,8), targets.use = c(1:6,9), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# netVisual_chord_gene(object.list[[1]], sources.use = c(7,8), targets.use = c(1:6,9), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

computeEnrichmentScore(net.down, species = 'mouse', variable.both = TRUE)

targDir = "Paper_figs/TasksList/P28/Merged_2/OPC_ODC/"

# OPC_ODC
#Up
upFilt = net.up[net.up$source == "OPC_Inter" | net.up$source == "ODC",]

png(paste0(targDir, "Ligand_Enrichment_StressUp_OPC_ODC_Output.png"), height = 14, width = 14, units = 'in', res = 300)
computeEnrichmentScore(upFilt, species = 'mouse', variable.both = TRUE, measure  = "ligand")
dev.off()

png(paste0(targDir, "Signaling_Enrichment_StressUp_OPC_ODC_Output.png"), height = 14, width = 14, units = 'in', res = 300)
computeEnrichmentScore(upFilt, species = 'mouse', variable.both = TRUE, measure  = "signaling")
dev.off()

# Down
downFilt = net.down[net.down$source == "OPC_Inter" | net.down$source == "ODC",]

png(paste0(targDir, "Ligand_Enrichment_StressDown_OPC_ODC_Output.png"), height = 14, width = 14, units = 'in', res = 300)
computeEnrichmentScore(downFilt, species = 'mouse', variable.both = TRUE, measure  = "ligand")
dev.off()

png(paste0(targDir, "Signaling_Enrichment_StressDown_OPC_ODC_Output.png"), height = 14, width = 14, units = 'in', res = 300)
computeEnrichmentScore(downFilt, species = 'mouse', variable.both = TRUE, measure  = "signaling")
dev.off()