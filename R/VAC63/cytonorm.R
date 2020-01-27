library(readxl) 
library(CATALYST)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(SingleCellExperiment)
library(CytoNorm)
library(FlowSOM)
library(devtools)


library(BiocManager)
# This should pull all dependencies.
BiocManager::install("FlowSOM") 

# Then install latest dependencies from github, using devtools.
install.packages("devtools") 
library(devtools) #load it

# install_github("RGLab/flowWorkspace")
# install_github("RGLab/openCyto")
BiocManager::install("CytoML") 
install_github('saeyslab/CytoNorm')


`%!in%` = Negate(`%in%`)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

red_palette <- c(myPalette(100), rep(myPalette(100)[100], 200))

sc <- scale_colour_gradientn(colours = red_palette, limits=c(0, 8))




#setwd("C:/Users/Florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")

setwd("/home/florian/PhD/cytof/vac63c/t_cells/experiment_279034_files")
files <- list.files(pattern = "*.fcs")
dir <- getwd()

data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Type = stringr::str_match(files, "concat")[, 1],
                   Batch = stringr::str_match(files, "b[0-9]*")[, 1],
                   stringsAsFactors = FALSE)

data$Batch <- ifelse(is.na(data$Batch)==T, "b3", data$Batch)
data$Type <- ifelse(is.na(data$Type)==T, "Validation", "Train")






#### cytonorm vignette ####
# 
# dir <- system.file("extdata", package = "CytoNorm")
# files <- list.files(dir, pattern = "fcs$")
# data <- data.frame(File = files,
#                    Path = file.path(dir, files),
#                    Type = stringr::str_match(files, "_([12]).fcs")[, 2],
#                    Batch = stringr::str_match(files, "PTLG[0-9]*")[, 1],
#                    stringsAsFactors = FALSE)
# data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
# 
# train_data <- dplyr::filter(data, Type == "Train")
# validation_data <- dplyr::filter(data, Type == "Validation")
# 
# ff <- flowCore::read.FCS(data$Path[1])
# channels <- flowCore::colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
#                                      40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
#                                      49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30, 
#                                      28, 29, 25, 35)]
# transformList <- flowCore::transformList(channels,
#                                          cytofTransform)
# transformList.reverse <- flowCore::transformList(channels,
#                                                  cytofTransform.reverse)


#### my stuff ####
ff <- flowCore::read.FCS(data$Path[1])

proper_channels <- colnames(ff)[c(3,14:15,23:57,63,65)]


transformList <- flowCore::transformList(proper_channels,
                                         CytoNorm::cytofTransform)

transformList_reverse <- flowCore::transformList(proper_channels,
                                                 CytoNorm::cytofTransform.reverse)

fsom <- prepareFlowSOM(data$Path,
                       proper_channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

cvs <- CytoNorm::testCV(fsom, cluster_values = c(5, 10, 15)) 


Training the model

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 6000, 
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 10,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)





