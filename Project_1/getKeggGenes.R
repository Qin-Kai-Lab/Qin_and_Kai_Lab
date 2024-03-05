library(stringr)

setwd('C:/Users/abomb/Downloads')

df<-read.delim('KEGG_2019_Mouse', F, sep = '_')

dat <- data.frame(do.call('rbind', strsplit(as.character(df$V1),'\t\t',fixed=TRUE)))

#snps5<-strsplit(foo$X2[2], "\t")

#group2<-str_to_title(snps5[[1]])

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

colnames(dat)<-c('Pathway', 'Genes')

write.csv(dat, 'kegg_mouse_2019.csv', row.names = F)

write.table(dat, 'kegg_mouse_2019.tsv', sep='\t', col.names = T, row.names = F)
