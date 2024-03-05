library(Seurat)
library(Signac)
library(ggplot2)
library(gridExtra)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

atacInt<-readRDS('atacIntegrated_macs2')

DefaultAssay(atacInt)

source('../programs/renameClusters.R')

# check cell numbers between RNA and ATAC
atacCells<-rownames(atacInt@meta.data)
rnaCells<-rownames(RNA.combined.norm@meta.data)

length(atacCells)
length(rnaCells)

length(atacCells[(atacCells%in%rnaCells)])

RNA.combined.norm@meta.data$IntCellId<-rownames(RNA.combined.norm@meta.data)

rnaAtac<-subset(RNA.combined.norm, subset = IntCellId %in% atacCells)

# check that the oreder of cells is the same

identical(atacCells, rownames(rnaAtac@meta.data))

identical(atacInt@meta.data$gex_barcode, rnaAtac@meta.data$CellName)

# change order to doublecheck

atacCells_sort<-sort(atacCells)

identical(atacCells_sort, atacCells)

# add cell IDs from RNA to atac

atacInt@meta.data$Annotations<-rnaAtac@meta.data$Annotations

Idents(atacInt)<-atacInt$Annotations

#RunUMAP(atacInt, reduction = "integrated_lsi", dims = 2:30)

DimPlot(atacInt)

p1 <- DimPlot(atacInt, reduction = "umap", group.by = "dataset", label = F)
p2 <- DimPlot(atacInt, reduction = "umap", label = T, repel=T)
p3<-grid.arrange(p1,p2, ncol=2, nrow=1)

ggsave(paste0('AtacClustRnaLab_macs2_',curDate, ".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)

saveRDS(atacInt, file = 'atacIntegrated_macs2')