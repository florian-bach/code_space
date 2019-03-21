source("https://bioconductor.org/biocLite.R")
biocLite("cytofkit")
library(cytofkit) 
cytofkit_GUI()  



set.seed(100) #set as whatever
dir <- system.file('extdata',package='cytofkit')
file <- list.files(dir ,pattern='.fcs$', full=TRUE)
parameters <- list.files(dir, pattern='.txt$', full=TRUE)
res <- cytofkit(fcsFiles = file, 
                markers = parameters, 
                projectName = 'cytofkit_test',
                transformMethod = "cytofAsinh", 
                mergeMethod = "ceil",
                fixedNum = 500,                                    ## set at 500 for faster run
                dimReductionMethod = "tsne",
                clusterMethods = c("Rphenograph", "ClusterX"),    ## accept multiple methods
                visualizationMethods = c("tsne", "pca"),          ## accept multiple methods
                progressionMethod = "isomap",
                clusterSampleSize = 500,
                resultDir = getwd(),
                saveResults = TRUE, 
                saveObject = TRUE)





## Loading the FCS data:  
dir <- system.file('extdata',package='cytofkit')
file <- list.files(dir ,pattern='.fcs$', full=TRUE)
paraFile <- list.files(dir, pattern='.txt$', full=TRUE)
parameters <- as.character(read.table(paraFile, header = TRUE)[,1])

## Extract the expression matrix with transformation

data_transformed <- cytof_exprsExtract(fcsFile = file, 
                                       comp = FALSE, 
                                       transformMethod = "cytofAsinh")

# for multiple files use

combined_data_transformed <- cytof_exprsMerge(fcsFiles = file, comp=FALSE,
                                              transformMethod = "cytofAsinh",
                                              mergeMethod = "all")

## use clustering algorithm to detect cell subsets
## to speed up our test here, we only use 100 cells
data_transformed_1k <- data_transformed[1:100, ]

## run PhenoGraph
cluster_PhenoGraph <- cytof_cluster(xdata = data_transformed_1k, method = "Rphenograph")

#tsne
data_transformed_1k_tsne <- cytof_dimReduction(data=data_transformed_1k, method = "tsne")

#clusterX
cluster_ClusterX <- cytof_cluster(ydata = data_transformed_1k_tsne,  method="ClusterX")

#run flowsom with cluster 15 

cluster_FlowSOM <- cytof_cluster(xdata = data_transformed_1k, method = "FlowSOM", FlowSOM_k = 12)


data_1k_all <- cbind(data_transformed_1k, data_transformed_1k_tsne, 
                     PhenoGraph = cluster_PhenoGraph, ClusterX=cluster_ClusterX, 
                     FlowSOM=cluster_FlowSOM)
data_1k_all <- as.data.frame(data_1k_all)

#vis

## PhenoGraph plot on tsne
cytof_clusterPlot(data=data_1k_all, xlab="tsne_1", ylab="tsne_2", 
                  cluster="PhenoGraph", sampleLabel = FALSE)



## PhenoGraph cluster heatmap
PhenoGraph_cluster_median <- aggregate(. ~ PhenoGraph, data = data_1k_all, median)
cytof_heatmap(PhenoGraph_cluster_median[, 2:37], baseName = "PhenoGraph Cluster Median")