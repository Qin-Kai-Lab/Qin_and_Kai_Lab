
setwd("~/Wang/output")


targDir <- 'OPC_ODC/Monocle3/86PC/'
genes = read.csv(paste0(targDir,'OPC_ODC_MonocClust_86PC_top100G_TimeCor_2023-02-27.csv'))



# get motifs 
getMotifs = function(genes, curGenes, targDir) {
  dir.create(targDir, recursive = T, showWarnings = F)
  if (curGenes == "All") {
    genesFilt = unique(genes$gene_id[genes$q_value < 0.05])
  } else if (curGenes == "Positive") {
    genesFilt = unique(genes$gene_id[genes$q_value < 0.05 & genes$estimate > 0 ])
  } else if (curGenes == "Negative") {
    genesFilt = unique(genes$gene_id[genes$q_value < 0.05 & genes$estimate < 0 ])
  }
  curTable = paste0(targDir, curGenes, "_genes.txt")
  outDir = paste0(targDir, curGenes, "_genes/")
  write.table(genesFilt, curTable, row.names = F, col.names = F, quote = F)
  # broader region gave more enriched motifs
  cmd_str = paste0("findMotifs.pl ", curTable, " mouse ", outDir, " -start -400 -end 100 -len 8,10,12 -p 8")
  #cmd_str = paste0("findMotifs.pl ", curTable, " mouse ", outDir, " -p 8")
  system(cmd_str)
}

getMotifs(genes = genes, curGenes = "All", targDir = "Paper_figs/TasksList/P27/")
getMotifs(genes = genes, curGenes = "Positive", targDir = "Paper_figs/TasksList/P27/")
getMotifs(genes = genes, curGenes = "Negative", targDir = "Paper_figs/TasksList/P27/")