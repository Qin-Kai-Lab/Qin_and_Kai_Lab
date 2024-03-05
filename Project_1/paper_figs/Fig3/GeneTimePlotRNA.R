library(Seurat)
library(ggplot2)
library(monocle3)

targDir <- 'OPC_ODC/Monocle3/86PC/'

rnaDat = readRDS('integ_OPC_ODC')

# get time 
cds<-readRDS('OpcOdcInt_MonocClust_86PC')
timeDf = data.frame(cds$monocle3_pseudotime)
timeDf$Cells = row.names(timeDf)
colnames(timeDf)[1] = "PsedTime"

# get expression 
getExpr = function(Obj, timeDf, genesList) {
  combDat = data.frame()
  for (curGene in genesList) {
    curExpr = Obj@assays$RNA@data[curGene ,]
    q1 = colnames(Obj@assays$RNA$data)
    exprDat = data.frame(cbind(q1, curExpr))
    colnames(exprDat) = c("Cells", "Expression")
    timeExpr = plyr::join(exprDat, timeDf, by = "Cells", type = "left", match = "all")
    timeExpr$Gene = curGene
    combDat = rbind(combDat, timeExpr)
  }
  return(combDat)
}

# 
genes = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv'))

findTop = function(genesDf, topN) {
  genesDf = genesDf[order(genesDf$p_value, decreasing = F),]
  genesDf$Abs_effect = abs(genesDf$estimate)
  genesSel = head(genesDf, topN*2)
  genesSel = genesSel[order(genesSel$Abs_effect, decreasing = T),]
  for (i in topN:nrow(genesSel) ) {
    curGenes = unique(head(genesSel$gene_id, i))
    if (length(curGenes == topN)) {
      break
    }
  }
  return(curGenes)
}

topGenes = findTop(genesDf=genes, topN = 10)


combDat = getExpr(Obj=rnaDat, timeDf=timeDf, genesList=topGenes)
combDat$Expression = as.numeric(as.character(combDat$Expression))


curPlot = ggplot(combDat, aes(x = PsedTime, y = Expression, color = Gene)) +
  theme_classic() +
  geom_smooth(method = "lm", formula  = y ~ x , se = FALSE) +
  theme(text = element_text(size = 28)) 

curPlot 

topN = "10"
typeGenes = "RNA_Linear"
png(filename = paste0(targDir,"TimeExpression_", topN, "_", typeGenes, "_2023-12-17.png"), width = 22, height = 16, units = "in", res = 300)
print(curPlot )
dev.off()

## quadratics
genes = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_TimeCor_RNA_Poly_2023-12-15.csv'))

findTop = function(genesDf, topN) {
  genesDf = genesDf[order(genesDf$p_value, decreasing = F),]
  genesDf$Abs_effect = abs(genesDf$estimate)
  genesSel = head(genesDf, topN*3)
  genesSel = genesSel[order(genesSel$Abs_effect, decreasing = T),]
  for (i in topN:nrow(genesSel) ) {
    curGenes = unique(head(genesSel$gene_id, i))
    if (length(curGenes == topN)) {
      break
    }
  }
  return(curGenes)
}

topGenes = findTop(genesDf=genes, topN = 10)


combDat = getExpr(Obj=rnaDat, timeDf=timeDf, genesList=topGenes)
combDat$Expression = as.numeric(as.character(combDat$Expression))


curPlot = ggplot(combDat, aes(x = PsedTime, y = Expression, color = Gene)) +
  theme_classic() +
  geom_smooth(method = "lm", formula  = y ~ x + I(x^2), se = FALSE) +
  theme(text = element_text(size = 28)) 

curPlot 

topN = "10"
typeGenes = "RNA_Quadratic"
png(filename = paste0(targDir,"TimeExpression_", topN, "_", typeGenes, "_2023-12-17.png"), width = 22, height = 16, units = "in", res = 300)
print(curPlot )
dev.off()