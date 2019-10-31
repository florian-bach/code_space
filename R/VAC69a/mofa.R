library(dplyr)
library(tidyr)
library(ggplot2)
library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)

library(reticulate)

# config <- py_config()
# config$numpy


# setup python environment in anaconda first
use_condaenv("r-reticulate", required = TRUE)

# py_install("mofapy", envname = "r-reticulate", method="auto")


# data clean up: remove features with little/no variance



#####     CYTOF DATA 


#CD4+ Flowsom
datacd4 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_results/results/cluster_abundances.csv")
colnames(datacd4)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')
datacd4 <- data.frame("Feature"=paste("CD4_Cluster_", datacd4$ClusterID, sep=''), datacd4[3:20])
rownames(datacd4) <- datacd4$Feature
datacd4$Feature <- NULL
datacd4 <- as.matrix(datacd4)

#CD8 flowsom
#data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_abundances.csv")
datacd8 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_results/results/cluster_abundances.csv")
colnames(datacd8)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')
datacd8 <- data.frame("Feature"=paste("CD8_Cluster_", datacd8$ClusterID, sep=''), datacd8[3:20])
rownames(datacd8) <- datacd8$Feature
datacd8$Feature <- NULL
datacd8 <- as.matrix(datacd8)

cytof_data <- rbind(datacd4, datacd8)

cytof_data <- cytof_data*100
#, substr(colnames(cytof_data), 1,3), substr(colnames(cytof_data), 5,nchar(colnames(cytof_data))))
# rownames(cytof_data)[(nrow(cytof_data)-1):nrow(cytof_data)] <- c("Volunteer", "Timepoint")


### PLASMA DATA


# data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")
dataplasma <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
t_dataplasma <- t(dataplasma)#

colnames(t_dataplasma) <- paste(t_dataplasma[1,], t_dataplasma[2,], sep='_')

dataplasma <- data.frame(t_dataplasma[3:nrow(t_dataplasma),])

colnames(dataplasma) <- gsub("C.1", "Baseline", colnames(dataplasma))
colnames(dataplasma) <- gsub("v00", "V0", colnames(dataplasma))
colnames(dataplasma) <- gsub("T", "DoD", colnames(dataplasma))

#convert to numerical matri
plasma_data <- as.matrix(select(dataplasma, colnames(cytof_data))) # <3
class(plasma_data) <- "double" # <3

#put it together, names() makes view names
#log transform plasma data

log_plasma_data <- log10(plasma_data)

# log_plasma_data <- rbind(log_plasma_data, substr(colnames(log_plasma_data), 1,3), substr(colnames(log_plasma_data), 5,nchar(colnames(log_plasma_data))))
# rownames(log_plasma_data)[(nrow(log_plasma_data)-1):nrow(log_plasma_data)] <- c("Volunteer", "Timepoint")


###### clinical & biochem data

### cell counts

data <- read.csv("/Users/s1249052/PhD/clinical data/vac69a/haem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("V", substr(data$Volunteer, 6, 7), sep='')

short <- filter(data, timepoint %in% c("_C_1", "_T6", "_EP"))

short2 <- select(short, Volunteer, timepoint, wbc, platelets, neutrophils, lymphocytes, monocytes, eosinophils)
short2$timepoint <- as.character(short2$timepoint)

###  NAUGHTY because dod and ep1 arent the same, but let's jsut assume theyre roughly the same
short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "DoD.6"
short2$timepoint[short2$timepoint=="_EP"] <- "DoD"

short2 <- filter(short2, timepoint %in% c("Baseline", "DoD.6", "DoD"))


mofa_haem <- t(short2)
colnames(mofa_haem) <- paste(mofa_haem[1,], mofa_haem[2,], sep="_") 
mofa_haem <- mofa_haem[3:nrow(mofa_haem),]
mofa_haem <- as.matrix(select(data.frame(mofa_haem), colnames(cytof_data)))
class(mofa_haem) <- "double"

# haem_data <- log10(mofa_haem)
haem_data <- mofa_haem
haem_data <- ifelse(is.infinite(haem_data), NA, haem_data)
haem_data <- select(haem_data, colnames(cytof_data))


## biochem

data <- read.csv("/Users/s1249052/PhD/clinical data/vac69a/biochem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("V", substr(data$Volunteer, 6, 7), sep='')

short <- filter(data, timepoint %in% c("_C_1", "_T6", "_EP"))
short2 <- select(short,  Volunteer, timepoint, sodium, potassium, creatinine, bilirubin, alt, alkphos, albumin, ast, ggt)

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "DoD.6"
short2$timepoint[short2$timepoint=="_EP"] <- "DoD"

short2 <- filter(short2, timepoint %in% c("Baseline", "DoD.6", "DoD"))


mofa_biochem <- t(short2)
colnames(mofa_biochem) <- paste(mofa_biochem[1,], mofa_biochem[2,], sep="_") 
mofa_biochem <- mofa_biochem[3:nrow(mofa_biochem),]
mofa_biochem <- as.matrix(select(data.frame(mofa_biochem), colnames(cytof_data)))
class(mofa_biochem) <- "double"

# biochem_data <- log10(mofa_biochem)
biochem_data <- mofa_biochem
biochem_data <- ifelse(is.infinite(biochem_data), NA, biochem_data)
######     MetaData, rownames of this should be the names of the columns in the data matrices

MetaData <- rbind(
  substr(colnames(cytof_data), 1,3),
  substr(colnames(cytof_data), 5, nchar(colnames(cytof_data)))
         )
colnames(MetaData) <- colnames(cytof_data)                                            

MetaData <- t(MetaData)
colnames(MetaData) <- c("Volunteer", "Timepoint")
# 

######   time??

# time <- gsub("Baseline", 1, MetaData)
# time <- gsub("DoD.6", 20, time, fixed = T)
# time <- gsub("DoD", 14, time, fixed=T)
# 
# time <- t(data.frame(time[2,]))
# class(time) <- "numeric"
# rownames(time) <- "Time"
# time <- as.matrix(time)

######   MOFA BABYYYYYYY


moby <- list(cytof_data, log_plasma_data, haem_data, biochem_data)
names(moby) <- c("CyTOF", "Plasma", "Haem", "Biochem")

mae_vivax <- MultiAssayExperiment(
  experiments = moby, 
  colData = MetaData
)

#mofa baby
MOFAobject <- createMOFAobject(mae_vivax) #<333333 <33333





# getting default model training options from the vignette

TrainOptions <- getDefaultTrainOptions()
ModelOptions <- getDefaultModelOptions(MOFAobject)

DataOptions <- getDefaultDataOptions()
DataOptions$scaleViews <- TRUE

TrainOptions$DropFactorThreshold <- 0.01
TrainOptions$tolerance <- 0.01


MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

MOFAobject <- regressCovariates(
  object = MOFAobject,
  views = c("CyTOF","Plasma", "Haem", "Biochem"),
  covariates = MOFAobject@InputData@colData$Volunteer
)


n_inits <- 50



MOFAlist <- lapply(seq_len(n_inits), function(it) {
  
  TrainOptions$seed <- 2018 + it
  
  MOFAobject <- prepareMOFA(
    MOFAobject, 
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )
  
  runMOFA(MOFAobject)
})


compareModels(MOFAlist)
#higher elbo values preferred

compareFactors(MOFAlist)
# visualises robsustness between runs


MOFAobject <- selectModel(MOFAlist, plotit = FALSE)
MOFAobject

plotVarianceExplained(MOFAobject)


plotWeightsHeatmap(
  MOFAobject, 
  view = "Plasma", 
  factors = 1:3,
  show_colnames = T
)

plotFactorScatter(
  MOFAobject,
  factors = c(1,3),
  color_by = "Volunteer",      # color by the IGHV values that are part of the training data
  shape_by = "Timepoint"  # shape by the trisomy12 values that are part of the training data
)


plotWeights(
  MOFAobject, 
  view = "Haem", 
  factor = 1, 
  nfeatures = 5
)

plotDataHeatmap(
  MOFAobject, 
  view = "Haem", 
  factor = 1, 
  features = 20, 
  show_rownames = T
)


plotFactorScatters(
  MOFAobject,
  factors = 1:5,
  color_by = "Volunteer",
  shape_by = "Timepoint"
)


plotFactorBeeswarm(
  MOFAobject,
  factors = 1,
  color_by = "Timepoint"
)


r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total

# Variance explained by each factor in each view
head(r2$R2PerFactor)

# Plot it
plotVarianceExplained(MOFAobject)




