library(Seurat)
library(Signac)
library(MAST)
library(ggplot2)
library(data.table)

setwd("~/Wang/output")

#my_vector <- seq(from = 0, to = 20000, by = 500)
#outFile = "Paper_figs/Fig5/Pvals_Gpr17_20KB_by500_Mean_2024-02-15.csv"

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

CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
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

#
rna_count = atacFiltSub@assays$RNA@counts["Gpr17" ,]
predict_count =  atacFiltSub@assays$PredictActivity@counts["Gpr17" ,]

rna_data = atacFiltSub@assays$RNA@data["Gpr17" ,]
predict_data =  atacFiltSub@assays$PredictActivity@data["Gpr17" ,]

spTest = cor.test(rna_count, predict_count, method = "spearman")
prTest = cor.test(rna_data, predict_data, method = "pearson")
prTest_2 = cor.test(rna_data, predict_count, method = "pearson")

##

curDf = data.frame(Cell_ID = names(rna_data), RNA_EXPR=rna_data, ATAC_COUNT = predict_count, ATAC_DATA = predict_data)

model1 = lm(RNA_EXPR~ATAC_COUNT, data = curDf)
summary(model1)

model2 = lm(RNA_EXPR~ATAC_DATA, data = curDf)
summary(model2)

curDf$RNA_EXPR = as.numeric(as.character(curDf$RNA_EXPR))
curDf$ATAC_COUNT = as.numeric(as.character(curDf$ATAC_COUNT))

jitter_amount <- 0.1

# Create the plot with jitter
curPlot <- ggplot(curDf, aes(x = ATAC_COUNT, y = RNA_EXPR)) +
  geom_point(position = position_jitter(width = jitter_amount, height = jitter_amount)) +  # Scatter plot with jitter
  geom_smooth(method = "lm", se = F, linewidth = 3) +  # Add linear regression line
  labs(x = "Predicted Expression (ATAC)", y = "RNA Expression") + 
  theme_classic() +
  theme(text = element_text(size = 36))

# boxplot

curDf$ATAC_COUNT = as.character(curDf$ATAC_COUNT)

boxP= ggplot(curDf, aes(x=ATAC_COUNT, y=RNA_EXPR)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  )

ggplot(curDf, aes(x=ATAC_COUNT, y=RNA_EXPR)) + 
  geom_violin(trim=T)+
  geom_jitter(size = 0.5) + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 4,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

# Calculate sum and mean for each group
summary_data <- curDf %>%
  group_by(ATAC_COUNT) %>%
  summarise(sum_value = sum(RNA_EXPR),
            mean_value = mean(RNA_EXPR), median_value = median(RNA_EXPR), 
            se_value = sd(RNA_EXPR) / sqrt(n()), 
            sd_value = sd(RNA_EXPR))

# Create bar plot
ggplot(summary_data, aes(x = ATAC_COUNT)) +
  geom_bar(aes(y = mean_value), stat = "identity", fill = "black", alpha = 0.8) +
  geom_text(aes(y = mean_value, label = round(mean_value, 2)), vjust = -0.1, hjust = -0.08 ,color = "red", size = 8) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2, color = "blue") +
  labs(x = "ATAC COUNT", y = "RNA EXPRESSION") + 
  theme_classic() +
  theme(text = element_text(size = 36))

ggplot(summary_data, aes(x = ATAC_COUNT)) +
  geom_bar(aes(y = mean_value), stat = "identity", fill = "black", alpha = 0.8) +
  geom_text(aes(y = mean_value, label = round(mean_value, 2)), vjust = -0.5,color = "red", size = 8) +
  geom_point(aes(y = mean_value), color = "red", size = 4) +
  #geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2, color = "blue") +
  labs(x = "ATAC COUNT", y = "RNA EXPRESSION") + 
  theme_classic() +
  theme(text = element_text(size = 36))
