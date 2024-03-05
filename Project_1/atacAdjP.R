setwd('/home/flyhunter/Wang/output/atacMergeDE_origPeaks')

curDate<-Sys.Date()

dir.create('adj_p')

filesList<-list.files(pattern = '.csv')

for ( i in filesList){
  df<-read.csv(i)
  df$bonferroni_p<-p.adjust(df$p_val, method = "bonferroni")
  df$fdr_p<-p.adjust(df$p_val, method = "fdr")
  outFile<-sub("_[^_]+$", "_", i)
  write.csv(df, paste0('adj_p/', outFile, curDate, '.csv'))
}