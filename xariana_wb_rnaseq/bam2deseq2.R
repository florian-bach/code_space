library(GenomicFeatures)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)

hse <- GenomicFeatures::makeTxDbFromGFF("/Volumes/fbach/ncbi-genomes-2023-04-27/GCF_000001405.40_GRCh38.p14_genomic.gtf", format="gtf" )
exonsByGene <- exonsBy(hse, by="gene")


fls <- "/Volumes/fbach/STAR_output/Aligned.out.bam"


bamLst <- BamFileList( fls, yieldSize=100000 )

se <- summarizeOverlaps(exonsByGene, bamLst,
                         mode="Union",
                         singleEnd=FALSE,
                         ignore.strand=TRUE,
                         fragments=TRUE )


countdata <- assay(se)
head( countdata )

countdata[grep("IFNG", rownames(countdata))]
