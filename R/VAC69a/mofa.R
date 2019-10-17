library(dplyr)
library(tidyr)
library(ggplot2)

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

long_data$

spread_data <- spread(long_data, Volunteer, timepoint)





