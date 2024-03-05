library(readxl)

list.files("~/Downloads", pattern = "xlsx")

excel_file <- "~/Downloads/DiffExprRNAMonoc3ClustAdjP_2023-06-08, kj.xlsx"

# Get the sheet names in the Excel document
sheet_names <- excel_sheets(excel_file)

# Create an empty list to store the data frames for each sheet
data_list <- list()

# Loop through each sheet and read the data into a data frame
for (sheet in sheet_names) {
  data <- read_excel(excel_file, sheet = sheet)
  data_list[[sheet]] <- data
}

# Access the data frames for each sheet using the sheet names as keys


curOpc = sort(data_list[[2]]$Genes)

## compare genes
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

targDir = './Paper_figs/Fig5/'

markers = read.csv(paste0(targDir, "DiffExprRNAMonoc3ClustAdjP_2023-06-08.csv"))
clusters = unique(markers$Cluster)
curDate = Sys.Date()
dfFilt<-markers[(markers$fdr_p < 0.05),]

dfOpc = dfFilt[dfFilt$Cluster == "OPC" ,]
opc = sort(dfFilt$Genes[dfFilt$Cluster == "OPC"])
identical(opc, curOpc)
opc[!(opc%in%curOpc)]
curOpc[!(curOpc%in%opc)]
identical(length(opc), length(curOpc))

# Intermediate
interm = sort(dfFilt$Genes[dfFilt$Cluster == "Intermideate"])
curInterm = sort(data_list[[3]]$Genes)
identical(interm, curInterm)
interm[!(interm%in%curInterm)]
curInterm[!(curInterm%in%interm)]
identical(length(interm), length(curInterm))

#odc

odc = sort(dfFilt$Genes[dfFilt$Cluster == "ODC"])
curOdc = sort(data_list[[4]]$Genes)
identical(odc, curOdc)
odc[!(odc%in%curOdc)]
curOdc[!(curOdc%in%odc)]
identical(length(odc), length(curOdc))