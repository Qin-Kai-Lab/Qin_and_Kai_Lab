library(Seurat)
library(ggplot2)

setwd("~/Projects/Project_8/output")

objSubFilt = readRDS("RegNormInteg_selFov.rds")

#ElbowPlot(objSubFilt, ndims = 50, reduction = "integrated.dr")

library(future)
plan("multicore", workers = 14)

options(future.globals.maxSize = 50 * 1024 ^ 3)

objSubFilt <- JackStraw(objSubFilt, num.replicate = 1000, dims = 50, reduction = "pca")
objSubFilt <- ScoreJackStraw(objSubFilt, dims = 1:50, reduction = "pca")
contrPlot = JackStrawPlot(objSubFilt, dims = 1:50, reduction = "pca")

ggsave(paste0("IntegRegNormJackStraw.png"), plot=contrPlot, height = 12, width = 14, units = 'in', dpi = 300)