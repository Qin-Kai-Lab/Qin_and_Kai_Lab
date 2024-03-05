filesList = list.files("Paper_figs/Fig1", pattern = "DotPlot", full.names = T)
basename(filesList)


dirList = list.files("Paper_figs", full.names = T)

for (curDir in dirList) {
  try({
    
    newDir = paste0(curDir, "/DotPlot")
    dir.create(newDir, recursive = T)
    filesList = list.files(curDir, pattern = "DotPlot", full.names = T)
    for (curFile in filesList) {
      newFileName = file.path(newDir, basename(curFile))
      file.rename(curFile, newFileName)
    }
    
  })
  
}



library(KEGGREST)

pathways <- keggList("pathway", organism = "mmu")
pathway_id <- "mmu00010"
pathway_info <- keggGet(pathway_id)

q1 = data.frame(pathway_info[[1]]$GENE)