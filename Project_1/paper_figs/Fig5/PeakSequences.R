library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
setwd("~/Wang/output")

peak1 = StringToGRanges("chr18-32147000-32148243")

genome_sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peak1)
peak1Str = as.character(genome_sequence)
write.table(peak1Str, "Paper_figs/Fig5/chr18-32147000-32148243_seq.txt", row.names = F, col.names = F, quote = F)


peak1 = StringToGRanges("chr18-31788558-31789864")

genome_sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peak1)
peak1Str = as.character(genome_sequence)
write.table(peak1Str, "Paper_figs/Fig5/chr18-31788558-31789864_seq.txt", row.names = F, col.names = F, quote = F)