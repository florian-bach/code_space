library(dplyr)
library(tidyr)
library(ggplot2)
library(MOFA) # <3

# read in data (laptop)

#CD4+ Flowsom
data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_abundances.csv")
colnames(data)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')

mof <- data.frame("Feature"=paste("CD4_Cluster_", data$ClusterID, sep=''), data[3:20])

#CD8 flowsom
data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_abundances.csv")
colnames(data)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')

data$Feature <- paste("CD8_Cluster_", data$ClusterID, sep='')
cd8 <-data.frame("Feature"=data[,21], data[,3:20])

mof <- rbind(cd8, mof) 


### plasma

data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:ncol(data)])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)
long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)

long_data$Volunteer <- gsub("v00", "V0", long_data$Volunteer, fixed = T)
long_data$timepoint <- gsub("+", ".", long_data$timepoint, fixed = T)
long_data$timepoint <- gsub("C-1", "Baseline", long_data$timepoint, fixed = T)


long_data$col <- paste(long_data$Volunteer, long_data$timepoint, sep = "_")

spread_data <- data.frame(t(spread(long_data, Analyte, Concentration)))
spread_data$Feature <- rownames(spread_data)

colnames(spread_data) <- unname(as.matrix(spread_data)[3,])
spread_data$col <- NULL
spread_data$Feature <- rownames(spread_data)

spread_data <- spread_data[4:nrow(spread_data),] 
spread_data <- data.frame(spread_data$Feature, spread_data[,1:ncol(spread_data)-1])
spread_data <- as.matrix(spread_data)
colnames(spread_data)[1] <- "Feature"

try <- as.data.frame(spread_data, stringsAsFactors = FALSE)
feature <- try$Feature
try$Feature <- NULL
try2 <- data.frame(sapply(try, FUN=function(x){as.numeric(x)}))
try2$Feature <- feature
colnames(try2)[22] <- colnames(mof)[10]

mof2 <- bind_rows(try2, mof)



library(MOFAdata)
data("CLL_data")
