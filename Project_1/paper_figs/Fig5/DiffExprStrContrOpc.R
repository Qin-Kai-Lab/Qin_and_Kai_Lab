library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)

editClusterNames = function(meta) {
  meta$newMonocClust = as.character(meta$newMonocClust)
  meta$Cluster = meta$newMonocClust
  for ( i in 1:nrow(meta) ) {
    if ( meta$Cell_ID[i] %in%opc_inter  ) {
      meta$Cluster[i] = "OPC_Inter"
    } else if (meta$Cell_ID[i]%in%odc) {
      meta$Cluster[i] = "ODC"
    } else {
      meta$Cluster[i] = meta$newMonocClust[i]
    }
  }
  return(meta)
}

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

opcOdc =  readRDS('integ_OPC_ODC')
opcOdc$newMonocClust = opcOdc$MonocClust
opcOdc$newMonocClust[opcOdc$newMonocClust == 4] = 3
opcOdc$newMonocClust[opcOdc$newMonocClust == 1] = "ODC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 2] = "OPC"
opcOdc$newMonocClust[opcOdc$newMonocClust == 3] = "Intermideate"

opcOdc$Cell_ID = rownames(opcOdc@meta.data)

atacFiltSub = subset(atacFilt, Merged_CellName%in%rownames(opcOdc@meta.data))
opcOdcSub = subset(opcOdc, Cell_ID%in%rownames(atacFiltSub@meta.data))

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

opc_inter = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "OPC" | opcOdcSub@meta.data$newMonocClust == "Intermideate"]
odc = rownames(opcOdcSub@meta.data)[opcOdcSub@meta.data$newMonocClust == "ODC" ]

opcOdcSub@meta.data = editClusterNames(opcOdcSub@meta.data)

identical(atacFiltSub$Merged_CellName, opcOdcSub$Cell_ID)
identical(colnames(atacFiltSub), colnames(opcOdcSub ))

atacFiltSub$newMonocClust = opcOdcSub$newMonocClust
atacFiltSub$Cluster = opcOdcSub$Cluster

subObj = subset(atacFiltSub, newMonocClust == "OPC")

Idents(subObj) = subObj$dataset

da_peaks <- FindMarkers(
  object = subObj,
  ident.1 = "Stress",
  ident.2 = "Control",
  test.use = 'LR',
  latent.vars = 'nCount_Combined_peaks',
  logfc.threshold = 0, only.pos = F)

da_peaks$Peak = rownames(da_peaks)

genes = ClosestFeature(subObj, regions = da_peaks$Peak)

combDat = cbind(da_peaks, genes)

targDir = './Paper_figs/Fig5/'

write.csv(combDat, 'Paper_figs/Fig5/DiffExprStrContrOpc_2024-02-23.csv', row.names = F)

gpr = combDat[combDat$gene_name == "Gpr17",]

annotation <- Annotation(object = subObj[["Combined_peaks"]])
curAnnot = annotation[annotation$gene_name == "Gpr17",]
transcripts <- CollapseToLongestTranscript(ranges = curAnnot)
transcripts <- Extend(x = transcripts, upstream = 2000, 
                      downstream = 2000)

gr <- GRanges(
  seqnames = sapply(strsplit(combDat$Peak, "-"), `[`, 1),
  ranges = IRanges(
    start = as.numeric(sapply(strsplit(combDat$Peak, "-"), `[`, 2)),
    end = as.numeric(sapply(strsplit(combDat$Peak, "-"), `[`, 3))
  )
)


#
# Extract overlapping regions from the original GRanges object
overlapping_regions <- subsetByOverlaps(gr,transcripts)

gprOverlapPeaks = c("chr18-31943368-31943733", "chr18-31950032-31950341")

gprOverlap = combDat[combDat$Peak%in%gprOverlapPeaks,]


