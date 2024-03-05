library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(clustree)
library(DoubletFinder)

## processing control

#please set the working directory after cleaning the memory
setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

control = Read10X(data.dir = "../data/stress")

RNA_control = CreateSeuratObject(counts = control$`Gene Expression`)
RNA_control$group = "Stress"

rm(control)

RNA_control[["percent.mt"]] <- PercentageFeatureSet(RNA_control, pattern = "^mt-")

VlnPlot(RNA_control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# calculate median and 2 standartd deviations
control_median<-median(RNA_control$nFeature_RNA)
control_sd<-sd(RNA_control$nFeature_RNA)
control_median+(control_sd*2)

contrMedRna<-median(RNA_control$nCount_RNA)
contrSdRna<-sd(RNA_control$nCount_RNA)
contrMedRna-(contrSdRna*2)

rnaContrFiltr<- subset(RNA_control, subset = nFeature_RNA > 200 & percent.mt < 5 & nCount_RNA > 400)
rm(RNA_control)
rnaContrFiltr<-NormalizeData(object = rnaContrFiltr)
rnaContrFiltr<-FindVariableFeatures(object = rnaContrFiltr)
rnaContrFiltr<-ScaleData(object = rnaContrFiltr)
rnaContrFiltr<-RunPCA(object = rnaContrFiltr)
ElbowPlot(rnaContrFiltr, ndims = 50)
rnaContrFiltr<-FindNeighbors(object = rnaContrFiltr, dims = 1:40)
rnaContrFiltr<-FindClusters(object = rnaContrFiltr)
rnaContrFiltr<-RunUMAP(object = rnaContrFiltr, dims = 1:40)

# pk
sweep.res.list_nsclc<- paramSweep_v3(rnaContrFiltr, PCs = 1:40, sct = F)
sweep.stats_nsclc<-summarizeSweep(sweep.res.list_nsclc, GT = F)
bcmvn_nsclc<- find.pK(sweep.stats_nsclc)

# plotting is not nessesary
ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1))+
  geom_point() +
  geom_line()

pK<- bcmvn_nsclc %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pk<-as.numeric(as.character(pK[[1]]))

nExp_poi<-round(0.078*nrow(rnaContrFiltr@meta.data))
annotations<-rnaContrFiltr@meta.data$seurat_clusters
homotypic.pop<-modelHomotypic(annotations)
nExp_poi.adj<- round(nExp_poi*(1-homotypic.pop))

# run doublet finder

rnaContrFiltr<-doubletFinder_v3(rnaContrFiltr, PCs = 1:40,
                                pN=0.25,
                                pK=pk,
                                nExp = nExp_poi.adj,
                                reuse.pANN = F, sct = F)

table(rnaContrFiltr@meta.data$DF.classifications_0.25_0.28_628)

rnaContrFiltr@meta.data$CellID<-rownames(rnaContrFiltr@meta.data)

doubletFiltr<-subset(x = rnaContrFiltr, subset = DF.classifications_0.25_0.28_628 == "Singlet")

singletCellId<-doubletFiltr@meta.data$CellID

write.table(singletCellId, 'cellNames_doubletFinderFiltered_stress.txt', row.names = F, col.names = F, quote = F)