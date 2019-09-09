library(ImpulseDE2)
library(splines)
library(DESeq2)
library(BiocParallel)
library(ggplot2)
## simulate data

lsSim <- simulateDataSetImpulseDE2(
  vecTimePointsA=rep(1:6, 3),
  vecTimePointsB=NULL,
  vecBatchesA=rep("B_NULL", 6*3),
  vecBatchesB=NULL,
  scaNConst=100,
  scaNImp=50,
  scaNLin=50,
  scaNSig=50,
  scaNRand=0,
  scaMumax=1000,
  scaSeedInit = 1)

## lsSim[1] is $dfAnnotation, a small dataframe containing sample name, condition (cas/control), Time and batch identifiers; SAMPLE TAGS
## lsSim[2] is $ matObservedCounts, a matrix where every row is a gene and every column is a sample


## run DESeq2
dfAnnotationDESeq2 <- lsSim$dfAnnotation

# create a natural cubic spline basis with 4 degrees of freedom; you should have less DoF than timepoints
matTimeSplineBasis <- ns(
  dfAnnotationDESeq2$Time, df=4)

colnames(matTimeSplineBasis) <- paste0("spline", seq(1, dim(matTimeSplineBasis)[2]))



# add spline basis into annotation data frame;
# now the annotation df has 4 more columns with an integer between -1 and 1 in every row, according to the splines from above
dfAnnotationDESeq2 <- cbind(dfAnnotationDESeq2, matTimeSplineBasis)


# create DESeq2 object with spline basis vectors as individual linear predictors
dds <- DESeqDataSetFromMatrix(
  countData = lsSim$matObservedCounts, # matrix of genes x samples
  colData = dfAnnotationDESeq2, # sample tags
  design = ~spline1 + spline2 + spline3 + spline4) # why do we add all splines? these need to be  colnames of coldata)


#plot those splines

broad_splines <- data.frame(cbind(dds$spline1, dds$spline2, dds$spline3, dds$spline4))
broad_splines$time <- 1:18
colnames(broad_splines)[1:4] <- paste0("spline", seq(1, dim(matTimeSplineBasis)[2]))

long_splines <- tidyr::gather(broad_splines, spline, y, colnames(broad_splines)[1:4])

ggplot(long_splines, aes(x=time, y=y))+
  geom_line()+
  facet_wrap(~spline)







## converting counts to integer mode

# estimate dispersion with less stringent dispersion outlier handling to avoid missing differentially expressed genes
# with observations that are zero; the outlier calling is controlled via outlierSD of estimateDispersionsMAP()

dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)

## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.

dds <- estimateDispersionsMAP(dds, outlierSD = 10)
# perform a log-likelihood ratio test
dds <- nbinomLRT(dds, full = ~spline1 + spline2 + spline3 + spline4,
                 reduced = ~1)
# extract results
print(head(results(dds), n=250))






###########       with batch correction         ############








# With batch correction:
  
suppressPackageStartupMessages(library(ImpulseDE2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
## simulate data
lsSim <- simulateDataSetImpulseDE2(
  vecTimePointsA=rep(1:6, 3),
  vecTimePointsB=NULL,
  vecBatchesA=do.call(c, lapply(
    seq(1,3), function(i) rep(paste0("B",i), 6) )),
  vecBatchesB=NULL,
  scaMuBatchEffect=1,
  scaSDBatchEffect=0.2,
  scaNConst=100,
  scaNImp=50,
  scaNLin=50,
  scaNSig=50,
  scaNRand=0,
  scaMumax=1000,
  scaSeedInit = 1)
## run DESeq2
dfAnnotationDESeq2 <- lsSim$dfAnnotation
# create a natural cubic spline basis with 4 degrees of freedom
matTimeSplineBasis <- ns(
  dfAnnotationDESeq2$Time, df=4)
colnames(matTimeSplineBasis) <- 
  paste0("spline", seq(1, dim(matTimeSplineBasis)[2]))
# add spline basis into annotation data frame
dfAnnotationDESeq2 <- cbind(dfAnnotationDESeq2, matTimeSplineBasis)
# create DESeq2 object with spline basis vectors as individual linear predictors
dds <- DESeqDataSetFromMatrix(
  countData = lsSim$matObservedCounts,
  colData = dfAnnotationDESeq2,
  design = ~Batch + spline1 + spline2 + spline3 + spline4)

## converting counts to integer mode

## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors

# estimate dispersion with less stringent dispersion outlier handling
# the outlier calling is controlled via outlierSD of estimateDispersionsMAP()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)

## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.

dds <- estimateDispersionsMAP(dds, outlierSD = 10)
# perform a log-likelihood ratio test
dds <- nbinomLRT(dds, full = ~Batch + spline1 + spline2 + spline3 + spline4,
                 reduced = ~Batch)
# extract results
print(head(results(dds)))

## log2 fold change (MLE): spline4 
## LRT p-value: '~ Batch + spline1 + spline2 + spline3 + spline4' vs '~ Batch' 
## DataFrame with 6 rows and 6 columns
##        baseMean log2FoldChange     lfcSE       stat     pvalue       padj
##       <numeric>      <numeric> <numeric>  <numeric>  <numeric>  <numeric>
## gene1  332.2238   0.0879040309 0.2924181  0.9434907 0.91824014 0.96535590
## gene2  422.1556   0.9170551773 0.3212048 13.1827163 0.01041667 0.03364355
## gene3  780.6804   0.2937970080 0.3761408  0.9062017 0.92366941 0.96535590
## gene4 1211.7727  -0.4093646633 0.3156516  5.8675786 0.20925825 0.35368283
## gene5  239.4681  -0.1719123830 0.3563992  1.8938956 0.75526629 0.89486528
## gene6 1102.1868   0.0006391691 0.2693686 12.9481841 0.01153177 0.03649294




#############################################







data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
data[,1:2] <- NULL

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:ncol(data)])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)
long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)

# make unique identifier of each sample
long_data$file_name <- paste(long_data$Volunteer, long_data$timepoint, sep = '_')

data2 <- spread(long_data, Analyte, Concentration)
data2[,1:2] <- NULL

# transpose data2 to make matrix with analyte rows and sample columns
data3 <- t(data2[,2:ncol(data2)])
colnames(data3) <- data2$file_name

plasma_data_spline <- data.frame(data3)
hello <- apply(plasma_data_spline, 2, FUN=function(x){as.integer(x)})
rownames(hello) <- rownames(plasma_data_spline)

hello_01 <- apply(hello, 1, FUN=function(x){scales::rescale(x, t=c(100, 1000))})
hello_01 <- apply(hello_01, 1, FUN=function(x){as.integer(x)})
rownames(hello_01) <- rownames(hello)



sample_tags <- data.frame(Sample=colnames(data3))


time  <- data.frame(regx=as.character(c("C-1", "DoD", "T", "45")), time= 1:4)
time$regx <- as.character(time$regx)


# stringr::str_extract(colnames(data3), "C-1|DoD|T|45")
Timer <- regmatches(colnames(data3),regexpr("C-1|DoD|T|45",colnames(data3)))

# this is fragile, you should fix it; struggle to pattern match substring with list entry
sample_tags$Time <- match(Timer, time$regx)


plasma_splines <- ns(
  sample_tags$Time, df=3)

colnames(plasma_splines) <- paste0("spline", seq(1, dim(plasma_splines)[2]))


sample_tags <- cbind(sample_tags, plasma_splines)
sample_tags$Volunteer <- substr(sample_tags$Sample, 1, 4)


plasma_dds <- DESeqDataSetFromMatrix(
  countData = hello_01, # matrix of genes x samples
  colData = sample_tags, # sample tags
  design = ~ Volunteer +spline3)



plasma_dds <- estimateSizeFactors(plasma_dds)
plasma_dds <- estimateDispersionsGeneEst(plasma_dds)
plasma_dds <- estimateDispersionsFit(plasma_dds)

## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.

plasma_dds <- estimateDispersionsMAP(plasma_dds, outlierSD = 10)
# perform a log-likelihood ratio test
plasma_dds <- nbinomLRT(plasma_dds, full = ~spline1 + spline2 + spline3,
                 reduced = ~1)
# extract results
View(as.data.frame(results(resOrdered)))
