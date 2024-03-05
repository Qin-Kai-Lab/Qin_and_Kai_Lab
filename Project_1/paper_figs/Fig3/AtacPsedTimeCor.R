library(Seurat)
library(monocle3)
library(MASS)
library(SeuratWrappers)

library(limma)

atac = readRDS("atacIntegrated_macs2_2_RNA")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

cds<-readRDS('OpcOdcInt_MonocClust_86PC')

timeDf = data.frame(cds$monocle3_pseudotime)
timeDf$Cells = row.names(timeDf)
colnames(timeDf)[1] = "PsedTime"

DimPlot(atac)

any(colnames(RNA.combined.norm) == rownames(timeDf))
identical(colnames(RNA.combined.norm), rownames(timeDf))
identical(rownames(RNA.combined.norm@meta.data), rownames(timeDf))


RNA.combined.norm$PesdTime = timeDf$cds.monocle3_pseudotime

q1 = colnames(RNA.combined.norm@assays$RNA$counts)
identical(q1, rownames(RNA.combined.norm@meta.data))

curGene = "Xkr4"
curExpr = RNA.combined.norm@assays$RNA@data[curGene ,]

exprDat = data.frame(cbind(q1, curExpr))
colnames(exprDat) = c("Cells", "Expression")

combDat = plyr::join(exprDat, timeDf, by = "Cells", type = "left", match = "all")

model <- glm.nb(Expression ~ PsedTime, data = combDat)

# Display model summary
summary(model)

## try with brining atac expression to monocle object
atacFilt = subset(atac, subset = Merged_CellName %in% colnames(RNA.combined.norm))
DefaultAssay(atacFilt) = "PredictActivity"
atacFilt@assays$Combined_peaks = NULL
atacFilt@assays$RNA = NULL

cds_atac <- as.cell_data_set(atacFilt)
# cds_atac <- preprocess_cds(cds_atac, num_dim = 86)
# cds_atac <- align_cds(cds_atac, num_dim = 86, alignment_group = "group")

metadat<-data.frame(cds_atac@colData)


#cds <- cluster_cells(cds, resolution=2e-3)
#cds <- cluster_cells(cds)

timeAtac = timeDf[timeDf$Cells%in%colnames(cds_atac),]
identical(timeAtac$Cells, colnames(cds_atac))
identical(timeAtac$Cells, rownames(metadat))

saveRDS(atacFilt, "OpcOdc_PredExpr")

#
rm(atac, atacFilt, RNA.combined.norm)
gc()

#
cds_atac$monocle3_pseudotime<- timeAtac$PsedTime
cds$monocle3_pseudotime<- pseudotime(cds)
print("started models")

saveRDS(cds_atac, "OpcOdc_PredExpr_Monoc3")

# gene_fits <- fit_models(cds_atac, model_formula_str = "~monocle3_pseudotime", expression_family="negbinomial")
# gene_Poly = fit_models(cds_atac, model_formula_str = "~monocle3_pseudotime+I(monocle3_pseudotime^2)", expression_family="negbinomial")
#rna_Poly = fit_models(cds, model_formula_str = "~monocle3_pseudotime+I(monocle3_pseudotime^2)", expression_family="negbinomial")

# fit_coefs <- coefficient_table(gene_fits)
# coefPoly = coefficient_table(gene_Poly)
#coefPolyRna = coefficient_table(rna_Poly)

# 
# coefFiltr<-data.frame(fit_coefs[(fit_coefs$status == 'OK') & (fit_coefs$term == 'monocle3_pseudotime'),])
# coefFiltr<-coefFiltr[, c(1:2, 5:ncol(coefFiltr))]
# 
# coefPolyFiltr = data.frame(coefPoly[(coefPoly$status == 'OK'),])
# coefPolyFiltr = coefPolyFiltr[coefPolyFiltr$term == 'monocle3_pseudotime' | coefPolyFiltr$term == 'I(monocle3_pseudotime^2)',]
# coefPolyFiltr<-coefPolyFiltr[, c(1:2, 5:ncol(coefPolyFiltr))]

# coefPolyFiltrRna = data.frame(coefPolyRna[(coefPolyRna$status == 'OK'),])
# coefPolyFiltrRna = coefPolyFiltrRna[coefPolyFiltrRna$term == 'monocle3_pseudotime' | coefPolyFiltrRna$term == 'I(monocle3_pseudotime^2)',]
# coefPolyFiltrRna<-coefPolyFiltrRna[, c(1:2, 5:ncol(coefPolyFiltrRna))]

# targDir <- 'OPC_ODC/Monocle3/86PC/'
# 
# write.csv(coefFiltr, paste0(targDir,'OPC_ODC_MonocClust_86PC_TimeCor_ATAC_2023-12-15.csv'), row.names = F)
# write.csv(coefPolyFiltr, paste0(targDir,'OPC_ODC_MonocClust_86PC_TimeCor_ATAC_Poly_2023-12-15.csv'), row.names = F)
#write.csv(coefPolyFiltrRna, paste0(targDir,'OPC_ODC_MonocClust_86PC_TimeCor_RNA_Poly_2023-12-15.csv'), row.names = F)

print("ALL complete!")