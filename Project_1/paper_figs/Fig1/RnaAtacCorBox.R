library(ggplot2)

dfDef = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/Default/PredictGene_RNA_spearman_AllClusters_2023-10-23.csv")
#dfDef = dfDef[dfDef$fdr_p < 0.05,]
df20 = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/20kbUp_20kbDown/PredictGene_RNA_spearman_AllClusters_2023-10-23.csv")
#df20 = df20[df20$fdr_p < 0.05,]

dfDef$ClustMeth = paste0(dfDef$Cluster, "_2KB")
dfDef$Method = "2KB"
df20$Method = "20KB"

df20$ClustMeth = paste0(df20$Cluster, "_20KB")

combDat = rbind(dfDef, df20)
min(combDat$fdr_p)
combDat = combDat[!is.na(combDat$fdr_p),]

combDat$`-Log10P.adj` = log10(combDat$fdr_p) * -1

combDat = combDat[combDat$Cluster != "C-R",]

combDat$Method  = factor(combDat$Method , levels =  c('2KB', '20KB'))

p= ggplot(combDat, aes(x=Method, y=`-Log10P.adj`, fill = Method)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  ) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue", linewidth = 1.1) &
  scale_y_continuous(limits = c(0, 1300), breaks = c(1, 10, 100, 1000) , labels = c(1, 10, 100, 1000))

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

ggsave(file = "Paper_figs/Fig1/Boxp_RNA_ATAC_Cor_fdrp_2023-12-12.png", plot = p, height = 10, width = 22, units = "in", dpi = 300)

# anudjasted pvaluse
dfDef = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/Default/PredictGene_RNA_spearman_AllClusters_2023-10-23.csv")
df20 = read.csv("/home/flyhunter/Wang/output/atacRna/PredictGeneActCor/20kbUp_20kbDown/PredictGene_RNA_spearman_AllClusters_2023-10-23.csv")

dfDef$ClustMeth = paste0(dfDef$Cluster, "_2KB")
dfDef$Method = "2KB"
df20$Method = "20KB"

df20$ClustMeth = paste0(df20$Cluster, "_20KB")

combDat = rbind(dfDef, df20)
combDat$`-Log10P` = log10(combDat$Pval) * -1

combDat = combDat[combDat$Cluster != "C-R",]
combDat$Method  = factor(combDat$Method , levels =  c('2KB', '20KB'))
combDat$LogP = log(combDat$`-Log10P`)
combDat$ScaleLogP = scale(combDat$`-Log10P`)

p= ggplot(combDat, aes(x=Method, y=`-Log10P`, fill = Method)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  ) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue", linewidth = 1.1) &
  scale_y_continuous(limits = c(0, 1300), breaks = c(0, 1, 10, 100, 1000) , labels = c(0, 1, 10, 100, 1000))

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

ggsave(file = "Paper_figs/Fig1/Boxp_RNA_ATAC_Cor_Pval_2023-12-12.png", plot = p, height = 10, width = 22, units = "in", dpi = 300)

p= ggplot(combDat, aes(x=Method, y=`-Log10P`, fill = Method)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  ) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue", linewidth = 1.1) &
  scale_y_continuous(limits = c(0, 300), breaks = c(1, 10, 100, 200, 300) , labels = c(1, 10, 100, 200, 300))

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

ggsave(file = "Paper_figs/Fig1/Boxp_RNA_ATAC_Cor_Pval_2023-12-12_cust.png", plot = p, height = 10, width = 22, units = "in", dpi = 300)



p= ggplot(combDat, aes(x=Method, y=LogP, fill = Method)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  ) +
  geom_hline(yintercept = 0.263, linetype = "dashed", color = "blue", linewidth = 1.1)

dir.create("Paper_figs/TasksList/P26/", recursive = T, showWarnings = F)
ggsave(file = "Paper_figs/TasksList/P26/Boxp_RNA_ATAC_Cor_Pval_2024-01-02_Log.png", plot = p, height = 10, width = 12, units = "in", dpi = 300)

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))
ggsave(file = "Paper_figs/TasksList/P26/Boxp_RNA_ATAC_Cor_Pval_Group_2024-01-03_Log.png", plot = p, height = 10, width = 12, units = "in", dpi = 300)

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

ggsave(file = "Paper_figs/Fig1/Boxp_RNA_ATAC_Cor_Pval_2023-12-12_Log.png", plot = p, height = 10, width = 22, units = "in", dpi = 300)


p= ggplot(combDat, aes(x=Method, y=ScaleLogP, fill = Method)) + 
  geom_boxplot() + theme_classic() +  theme(
    text = element_text(size = 26),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
  ) +
  geom_hline(yintercept = -0.38, linetype = "dashed", color = "blue", linewidth = 1.1)

p = p + facet_grid(cols = vars(Cluster))  + theme(axis.text.y = element_text(size = 12))

ggsave(file = "Paper_figs/Fig1/Boxp_RNA_ATAC_Cor_Pval_2023-12-12_Scale.png", plot = p, height = 10, width = 22, units = "in", dpi = 300)


q1 = combDat[combDat$Pval > 0.05,]
table(q1$Cluster)
q2 = combDat[combDat$Pval <0.05,]
table(q2$Cluster)