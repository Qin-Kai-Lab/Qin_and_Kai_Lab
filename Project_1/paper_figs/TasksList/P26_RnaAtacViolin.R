library(ggplot2)

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


# p= ggplot(combDat, aes(x=Method, y=LogP, fill = Method)) + 
#   geom_boxplot() + theme_classic() +  theme(
#     text = element_text(size = 26),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 18)  # Adjust the angle and size
#   ) +
#   geom_hline(yintercept = 0.263, linetype = "dashed", color = "blue", linewidth = 1.1)
# 
# dir.create("Paper_figs/TasksList/P26/", recursive = T, showWarnings = F)
# ggsave(file = "Paper_figs/TasksList/P26/Boxp_RNA_ATAC_Cor_Pval_2024-01-02_Log.png", plot = p, height = 10, width = 12, units = "in", dpi = 300)

p = ggplot(combDat, aes(x=Method, y=LogP, fill = Method)) + 
  geom_violin(trim=FALSE)+
  geom_jitter(size = 0.5) + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 4,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26)) +
  geom_hline(yintercept = 0.263, linetype = "dashed", color = "blue", linewidth = 1.1)

ggsave(file = "Paper_figs/TasksList/P26/Violin_RNA_ATAC_Cor_Pval_2024-01-02_Log.png", plot = p, height = 10, width = 12, units = "in", dpi = 300)
