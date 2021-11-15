library(Seurat)

hello <- Seurat::Read10X('~/postdoc/scRNAseq/debarcode_try/read_count/', gene.column=1)

hello_df <- data.frame(t(hello[1:6, 1:350]))

hist(hello_df[,6])
