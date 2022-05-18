library(cellhashR)

# Example 1: parse CITE-seq-Count output, printing QC
barcodeData <- ProcessCountMatrix(rawCountData = '~/postdoc/scRNAseq/CITE_seq_count_output/primary/umi_count', minCountPerCell = 5)


# Create QC plots of barcode normalization
#PlotNormalizationQC(barcodeData)

# Generate the final cell hashing calls
calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = c('bff_cluster'))

# Inspect negative cells:
SummarizeCellsByClassification(calls = calls, barcodeMatrix = barcodeData)
