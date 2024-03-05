library(ArchR)
library(ggplot2)



addArchRGenome("mm10")

addArchRThreads(threads = 14) 

proj = loadArchRProject(path = "./archContrStress")

curDate<-Sys.Date()

addArchNames = function(x) {
  x$archName <- NA
  for (i in 1:nrow(x)) {
    if (x$group[i] == 'Control') {
      x$archName[i] = paste0('control#', x$CellName[i])
    } else if (x$group[i] == 'Stress') {
      x$archName[i] = paste0('stress#', x$CellName[i])
    }
  }
  return(x)
}


rnaMetadat = read.csv('RnaSeqMetadat.csv')
archMeta = data.frame(proj$cellNames)
colnames(archMeta)[1] = 'archName'

rnaMetadat = addArchNames(rnaMetadat)

atacRnaPres = merge(archMeta, rnaMetadat, by = 'archName', sort = F)
nrow(atacRnaPres)
length(archMeta$archName[archMeta$archName %in% rnaMetadat$archName])

projFilt = subsetCells(ArchRProj = proj, cellNames = atacRnaPres$archName)

identical(projFilt$cellNames, atacRnaPres$archName)

projFilt$Annotations = atacRnaPres$Annotations

projFilt <- addIterativeLSI(ArchRProj = projFilt, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
projFilt <- addClusters(input = projFilt, reducedDims = "IterativeLSI", force = T)

rm(proj)
gc()
##
projFilt <- addUMAP(ArchRProj = projFilt, reducedDims = "IterativeLSI", force = TRUE)
condGr <- gsub("\\#.*", "", projFilt$cellNames)
projFilt$Group <- condGr

targDir <- './archContrStress/atacRna/UMAPs/'
dir.create(targDir, recursive = T)


plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "Group", embedding = "UMAP")+
  theme_ArchR(
    color = "black",
    textFamily = "sans",
    baseSize = 20,
    baseLineSize = 0.5,
    baseRectSize = 0.5,
    plotMarginCm = 1,
    legendPosition = "bottom",
    legendTextSize = 20,
    axisTickCm = 0.1,
    xText90 = FALSE,
    yText90 = FALSE
  )

plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "Annotations", embedding = "UMAP")
annotUmap = getEmbedding(ArchRProj = projFilt, embedding = "UMAP", returnDF = TRUE)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#levels(x =projFilt) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

#projFilt$Annotations <- factor(projFilt$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))

archClusters = character()

for (i in projFilt$Annotations) {
  if (i == 'SUB') {
    j = '1-SUB'
  } else if (i == 'CA1' ) { 
      j = '2-CA1'
  } else if (i == 'CA2' ) {
      j = '3-CA2'
  } else if (i == 'CA3' ) {
      j = '4-CA3'
  } else if (i == 'DG') {
      j = '5-DG'
  } else if ( i == 'GABA') {
      j = '6-GABA'
  } else if (i == 'C-R') {
    j = '7-C-R'
  } else if (i == 'OPC') {
    j = '8-OPC'
  } else if (i == 'ODC') {
    j = '9-ODC'
  }  else if (i == 'MG') {
    j = '10-MG'
  }
  archClusters = c(archClusters, j)
}

projFilt$newClust = archClusters

# save Umaps not finished
fName <- paste0(targDir, 'atacFiltRna_Umap_annot', '_', curDate, '.png')
png(file= fName, width=14, height=14, units = 'in', res = 300)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "newClust", embedding = "UMAP")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)
dev.off()

fName <- paste0(targDir, 'atacFiltRna_Umap_Group', '_', curDate, '.png')
png(file= fName, width=14, height=14, units = 'in', res = 300)
plotEmbedding(ArchRProj = projFilt, colorBy = "cellColData", name = "Group", embedding = "UMAP")+
  theme_ArchR(baseSize = 20,legendTextSize = 20)
dev.off()

cellsPerClust = data.frame(table(projFilt$Group, projFilt$Annotations))
colnames(cellsPerClust)[1:2] = c('Group', 'Cluster')
write.csv(cellsPerClust, paste0(targDir, 'CellsPer_GroupCluster_', curDate, '.csv'), row.names = F)

# control vs stress

projFilt$Annotations = as.character(projFilt$Annotations)
projFilt$Group = as.character(projFilt$Group)

targDir <- './archContrStress/atacRna/StressVsControl/'
dir.create(targDir, recursive = T)

atacMark = function() {
  allClustSum = data.frame(matrix(nrow = 0, ncol = 0))
  clusters = unique(projFilt$Annotations)
  for ( cluster in clusters) {
    cellSubs = projFilt$cellNames[projFilt$Annotations == cluster]
    cellsClust = subsetCells(ArchRProj = projFilt, cellNames = cellSubs)
    atacMark <- getMarkerFeatures(
      ArchRProj = cellsClust, 
      groupBy = "Group",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      useGroups = "stress",
      bgdGroups = "control",
      maxCells = 20000,
      useMatrix = "GeneScoreMatrix"
    )
    #markerList <- getMarkers(atacMark, cutOff = "Pval <= 1")
    #clust1Mark <- data.frame(markerList$stress)
    #clust1Mark$fdr_p <- atacMark@assays@data$FDR
    manSum = data.frame(cbind(atacMark@elementMetadata$name, atacMark@assays@data$Pval, atacMark@assays@data$FDR, atacMark@assays@data$Log2FC)) # checked manual sumary is identical with package
    colnames(manSum) = c('Genes', 'Pval', 'FDR_def', 'Log2FC')
    #sumSel = manSum[(manSum$Log2FC > 0),]
    manSum$bonferonni_p = p.adjust(p=manSum$Pval, method = 'bonferroni')
    manSum$fdr_p = p.adjust(p=manSum$Pval, method = 'fdr')
    manSum$Cluster = cluster
    outName = paste0(targDir, 'StressVsControl_GeneScoremat_', cluster, '_', curDate, '.csv')
    write.csv(manSum, outName, row.names = F)
    allClustSum = rbind(allClustSum, manSum)
  }
  return(allClustSum)
}

allClustMark = atacMark()
write.csv(allClustMark, paste0(targDir, 'StressVsControl_GeneScoremat_allClusters', curDate, '.csv'), row.names = F)

# check than manual summary is the same as archR
#q1 <- manSum[order(manSum$x.2),]
#q2 <- clust1Mark[order(clust1Mark$Log2FC),]
#identical(q1$atacMark.elementMetadata.name, q2$name)

#saveArchRProject(ArchRProj = projFilt, outputDirectory = "./archContrStress/atacRna/projFilt")

saveRDS(projFilt, 'atacArchRnaFilt')