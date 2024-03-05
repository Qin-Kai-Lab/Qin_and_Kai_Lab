clusters = c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

groups = c("Control", "Stress")


combGroups = function(gr1, gr2) {
  newGr = character()
  for ( cond in gr1) {
    for (cluster in gr2) {
      curCond = paste(cond, cluster, sep = "_")
      newGr = c(newGr, curCond)
    }
  }
  return(newGr)
}

