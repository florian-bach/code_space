
data <- read.csv("/Users/s1249052/PhD/cytof/vac69a/T_cells_only/csv/Volunteer_03_T_cells_T_6.csv")

wee_data <- head(data, n=100)

transposed_csv <- t(data)

write.csv(transposed_csv, "transposed_csv_03_t6.csv")

BiocManager::install("scater")

library(SC3)
library(scater)



#create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(transposed_csv),
    logcounts = log2(as.matrix(transposed_csv) + 1)
    ))

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]



sce <- sc3(sce, ks = 10:20, biology = TRUE)
Sys.time()
