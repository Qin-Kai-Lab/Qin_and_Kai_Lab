library(Seurat)
library(Signac)
library(ggplot2)
library(gridExtra)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

atacInt<-readRDS('atacIntegrated_macs2_stressFirst')

DefaultAssay(atacInt)

source('../programs/renameClusters.R')

# renmae rna 

cellType<- RNA.combined.norm
newNames<-gsub('_1', '_control', colnames(cellType))
newNames<-gsub('_2', '_stress', newNames)


cellType<-RenameCells(cellType, new.names = newNames)

# rename atac
nextNames<-gsub('_2', '_control', colnames(atacInt))
nextNames<-gsub('_1', '_stress', nextNames)


atacInt<-RenameCells(atacInt, new.names = nextNames)


# check cell numbers between RNA and ATAC
RNA.combined.norm<-cellType

atacCells<-colnames(atacInt)
rnaCells<-colnames(RNA.combined.norm)

length(atacCells)
length(rnaCells)

length(atacCells[(atacCells%in%rnaCells)])

RNA.combined.norm@meta.data$IntCellId<-colnames(RNA.combined.norm)

rnaAtac<-subset(RNA.combined.norm, subset = IntCellId %in% atacCells)

# check that the oreder of cells is the same

identical(atacCells, colnames(rnaAtac))

identical(atacInt@meta.data$gex_barcode, rnaAtac@meta.data$CellName)


rnaAtac@meta.data$Cluster<-rnaAtac@meta.data$Annotations
ind = match(colnames(rnaAtac), colnames(atacInt))
atacInt@meta.data$Annotations<-NA
atacInt@meta.data[ind, "Cluster"] = rnaAtac@meta.data[, "Cluster"]


# change order to doublecheck

#atacCells_sort<-sort(atacCells)

#identical(atacCells_sort, atacCells)

# add cell IDs from RNA to atac

#atacInt@meta.data$Annotations<-rnaAtac@meta.data$Annotations

atacInt$Annotations<-atacInt$Cluster
Idents(atacInt)<-atacInt$Annotations

#RunUMAP(atacInt, reduction = "integrated_lsi", dims = 2:30)

DimPlot(atacInt)

p1 <- DimPlot(atacInt, reduction = "umap", group.by = "dataset", label = F)
p2 <- DimPlot(atacInt, reduction = "umap", label = T, repel=T)
p3<-grid.arrange(p1,p2, ncol=2, nrow=1)

ggsave(paste0('AtacClustRnaLab_macs2_',curDate, ".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)

saveRDS(atacInt, file = 'atacIntegrated_macs2_stressFirst')