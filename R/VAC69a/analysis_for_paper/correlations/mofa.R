library(dplyr)
library(tidyr)
library(ggplot2)
<<<<<<< HEAD
library(MOFA) # <3
=======
library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
>>>>>>> 5992d616634ea1357ef19c85c4d7b4f781cc9ca1

library(reticulate)

# config <- py_config()
# config$numpy


# setup python environment in anaconda first
#use_condaenv("r-reticulate", required = TRUE)
use_python("/home/flobuntu/anaconda3/bin/python", required = TRUE)


#py_install("mofapy", envname = "r-reticulate", method="auto")


# data clean up: remove features with little/no variance



#####     CYTOF DATA ####

#consider doing a square root asinh transform, like the complex heatmap

cytof_data <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header = T, stringsAsFactors = F)
cytof_data <- subset(cytof_data, cytof_data$timepoint!="C10")

#cytof_data$trans_freq=scale(asin(sqrt(cytof_data$frequency/100)), center = TRUE, scale = TRUE)

cytof_data$trans_freq=asin(sqrt(cytof_data$frequency/100))


# wide_cytof <- data.frame(cytof_data %>%
#   select(cluster_id, sample_id, frequency) %>%
#   pivot_wider(names_from = sample_id, values_from = frequency))

wide_cytof <- data.frame(cytof_data %>%
                           select(cluster_id, sample_id, trans_freq) %>%
                           pivot_wider(names_from = sample_id, values_from = trans_freq))


rownames(wide_cytof) <- gsub("ï", "i", wide_cytof$cluster_id)
wide_cytof <- as.matrix(select(wide_cytof, -cluster_id))
class(wide_cytof) <- "double"


### haem ####

#data <- read.csv("/Users/s1249052/PhD/clinical data/vac69a/haem.csv")
data <- read.csv("~/PhD/clinical_data/vac69a/haem.csv", header=T, stringsAsFactors = F)


colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("V", substr(data$Volunteer, 6, 7), sep='')

short <- filter(data, timepoint %in% c("_C_1", "_T6", "_EP"))

short2 <- select(short, Volunteer, timepoint, wbc, platelets, neutrophils, lymphocytes, monocytes, eosinophils)
short2$timepoint <- as.character(short2$timepoint)

###  NAUGHTY because dod and ep1 arent the same, but let's jsut assume theyre roughly the same
short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "T6"
short2$timepoint[short2$timepoint=="_EP"] <- "DoD"

short2 <- filter(short2, timepoint %in% c("Baseline", "T6", "DoD"))


mofa_haem <- t(short2)
colnames(mofa_haem) <- paste(mofa_haem[1,], mofa_haem[2,], sep="_") 
mofa_haem <- mofa_haem[3:nrow(mofa_haem),]
mofa_haem <- as.matrix(mofa_haem[,colnames(mofa_haem)%in%colnames(wide_cytof)[2:ncol(wide_cytof)]]) # <3

class(mofa_haem) <- "double"
haem_data <- mofa_haem

log_haem_data <- log10(mofa_haem)
log_haem_data <- ifelse(is.infinite(log_haem_data), -10, log_haem_data)
#haem_data <- na.omit(log_haem_data)

## biochem ####

#data <- read.csv("/Users/s1249052/PhD/clinical data/vac69a/biochem.csv")
data <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")

colnames(data)[1] <- "Volunteer"

data$Volunteer <- paste("V", substr(data$Volunteer, 6, 7), sep='')

short <- filter(data, timepoint %in% c("_C_1", "_T6", "_EP"))
short2 <- select(short,  Volunteer, timepoint, sodium, potassium, creatinine, bilirubin, alt, alkphos, albumin, ast, ggt)

short2$timepoint <- as.character(short2$timepoint)

short2$timepoint[short2$timepoint=="_C_1"] <- "Baseline"
short2$timepoint[short2$timepoint=="_T6"] <- "T6"
short2$timepoint[short2$timepoint=="_EP"] <- "DoD"

short2 <- filter(short2, timepoint %in% c("Baseline", "T6", "DoD"))


mofa_biochem <- t(short2)
colnames(mofa_biochem) <- paste(mofa_biochem[1,], mofa_biochem[2,], sep="_") 
mofa_biochem <- mofa_biochem[3:nrow(mofa_biochem),]

mofa_biochem <- as.matrix(mofa_biochem[,colnames(mofa_biochem)%in%colnames(wide_cytof)]) # <3
class(mofa_biochem) <- "double"

#impuyte v05 NA at Baseline by assuming value at DoD
log_biochem_data <- mofa_biochem[1:7,]
log_biochem_data[5,1] <- log_biochem_data[5,17]
log_biochem_data <- log10(log_biochem_data)
#biochem_data <- ifelse(is.infinite(biochem_data), NA, biochem_data)

log_biochem_data <- log_biochem_data[c(3,5:7), ]




alt_timecourse <- mofa_biochem[5,order(colnames(mofa_biochem))]
alt_timecourse[7] <- 13
cor(t(wide_cytof), alt_timecourse)


cytof_data$alt <- alt_timecourse[match(cytof_data$sample_id, names(alt_timecourse))]

correlation_dfs <- split(cytof_data, cytof_data$cluster_id)
names(correlation_dfs) <- lapply(correlation_dfs, function(x)unique(x$cluster_id))
correlation_results <- lapply(correlation_dfs, function(x)cor.test(x$frequency, x$alt, method="pearson"))

corr_res <- lapply(correlation_results, function(x) cbind("p_value"=x$p.value, "rho"=x$estimate))
  
corr_res <- data.frame(do.call(rbind, corr_res))
rownames(corr_res) <- names(correlation_results)
corr_res$p_adj <- p.adjust(corr_res$p_value, "BH")

subset(corr_res, corr_res$p_adj<=0.05)

#> subset(corr_res, corr_res$p_adj<=0.05)
#                                       p_value        rho        p_adj
# activated  CD4 CM                 5.362225e-04 0.7331917 5.719707e-03
# activated  MAIT                   1.384127e-03 0.6944552 8.858416e-03
# activated  Vd2+                   1.134533e-03 0.7030937 8.858416e-03
# activated PD1+CD27-HLADR+  CD4 EM 1.139985e-04 0.7849951 1.823976e-03
# activated PD1+HLADR+ CD4 EM       1.851518e-06 0.8763470 5.924857e-05
# activated PD1+HLADR+ CD8 EM       7.048489e-03 0.6111414 3.759194e-02


######   time??

# time <- gsub("Baseline", 1, MetaData)
# time <- gsub("DoD.6", 20, time, fixed = T)
# time <- gsub("DoD", 14, time, fixed=T)
# 
# time <- t(data.frame(time[2,]))
# class(time) <- "numeric"
# rownames(time) <- "Time"
# time <- as.matrix(time)

### PLASMA DATA ####


# data <- read.csv("C:/Users/Florian/PhD/oxford/vac69/plasma_analytes_vivax_no_inequalitites.csv")
#dataplasma <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
dataplasma <- read.csv("~/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")
dataplasma <- subset(dataplasma, dataplasma$Volunteer!="v009")
t_dataplasma <- t(dataplasma)#

colnames(t_dataplasma) <- paste(t_dataplasma[1,], t_dataplasma[2,], sep='_')

dataplasma <- data.frame(t_dataplasma[3:nrow(t_dataplasma),])

colnames(dataplasma) <- gsub("C.1", "Baseline", colnames(dataplasma))
colnames(dataplasma) <- gsub("v00", "V0", colnames(dataplasma))
colnames(dataplasma) <- gsub("T", "DoD", colnames(dataplasma))
colnames(dataplasma) <- gsub("DoD.6", "T6", colnames(dataplasma))

#convert to numerical matri
plasma_data <- as.matrix(dataplasma[,colnames(dataplasma)%in%colnames(wide_cytof)]) # <3

changing_analytes <- sig_analytes <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt", what="", skip = 1)
changing_analytes <- changing_analytes[1:18]

rownames(plasma_data) <- gsub("pg.ml.", "", rownames(plasma_data), fixed = T)
rownames(plasma_data) <- gsub(".", "",rownames (plasma_data), fixed=T)

plasma_data <- subset(plasma_data, rownames(plasma_data)%in%changing_analytes)

class(plasma_data) <- "double" # <3

#put it together, names() makes view names
#log transform plasma data

log_plasma_data <- log2(plasma_data)

log_plasma_data <- rbind(log_plasma_data, "alt"=biochem_data["alt",!grepl("V05", colnames(biochem_data))])


# log_plasma_data <- rbind(log_plasma_data, substr(colnames(log_plasma_data), 1,3), substr(colnames(log_plasma_data), 5,nchar(colnames(log_plasma_data))))
# rownames(log_plasma_data)[(nrow(log_plasma_data)-1):nrow(log_plasma_data)] <- c("Volunteer", "Timepoint")



# metadata ####

######     MetaData, rownames of this should be the names of the columns in the data matrices

MetaData <- rbind(
  substr(colnames(wide_cytof), 1,3),
  substr(colnames(wide_cytof), 5, nchar(colnames(wide_cytof)))
)

colnames(MetaData) <- colnames(wide_cytof)

MetaData <- t(MetaData)
colnames(MetaData) <- c("Volunteer", "Timepoint")
# 


######   MOFA BABYYYYYYY #####



moby <- list(wide_cytof, log_plasma_data, log_haem_data, log_biochem_data)  #, haem_data, biochem_data
names(moby) <- c("CyTOF", "Plasma", "Haem", "Biochem") #, "Haem", "Biochem"

mae_vivax <- MultiAssayExperiment::MultiAssayExperiment(
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


#regress out personal identity

# MOFAobject <- regressCovariates(
#   object = MOFAobject,
#   views = c("CyTOF","Plasma", "Biochem"), #, #, "Haem", "Biochem"
#   covariates = MOFAobject@InputData@colData$Volunteer
# )


n_inits <- 5


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

# picks the model with highest elbo
MOFAobject <- selectModel(MOFAlist, plotit = FALSE)

plotVarianceExplained(MOFAobject)


plotFactorScatters(
  MOFAobject,
  factors = 1:5,
  color_by = "Volunteer",
  shape_by = "Timepoint"
)


plotWeightsHeatmap(
  MOFAobject, 
  view = "CyTOF", 
  factors = 1:5,
  show_colnames = T
)

plotWeightsHeatmap(
  MOFAobject, 
  view = "Haem", 
  factors = c(1,2, 3,4),
  show_colnames = T
)


plotFactorScatter(
  MOFAobject,
  factors = c(1,4),
  color_by = "Volunteer",      # color by the IGHV values that are part of the training data
  shape_by = "Timepoint"  # shape by the trisomy12 values that are part of the training data
)


plotWeights(
  MOFAobject, 
  view = "Haem", 
  factor = 1, 
  nfeatures = 12
)

plotDataHeatmap(
  MOFAobject, 
  view = "Plasma", 
  factor = 2, 
  features = 20, 
  show_rownames = T
)


plotFactorScatters(
  MOFAobject,
  factors = 1:4,
  color_by = "Volunteer",
  shape_by = "Timepoint"
)






###regression ####

MOFAobject <- regressCovariates(
  object = MOFAobject,
  views = c("CyTOF","Plasma")#, #, "Haem", "Biochem"
  #covariates = MOFAobject@InputData@colData$Volunteer
)


n_inits <- 50



regressed_MOFAlist <- lapply(seq_len(n_inits), function(it) {
  
  TrainOptions$seed <- 2018 + it
  
  MOFAobject <- prepareMOFA(
    MOFAobject, 
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )
  
  runMOFA(MOFAobject)
})


compareModels(regressed_MOFAlist)
#higher elbo values preferred

compareFactors(regressed_MOFAlist)
# visualises robsustness between runs


regressed_MOFAobject <- selectModel(regressed_MOFAlist, plotit = FALSE)
regressed_MOFAobject

plotVarianceExplained(MOFAobject)


plotWeightsHeatmap(
  MOFAobject, 
  view = "Plasma", 
  factors = 1:3,
  show_colnames = T
)

plotFactorScatter(
  regressed_MOFAobject,
  factors = c(1,2),
  color_by = "Volunteer",      # color by the IGHV values that are part of the training data
  shape_by = "Timepoint"  # shape by the trisomy12 values that are part of the training data
)


plotWeights(
  regressed_MOFAobject, 
  view = "CyTOF", 
  factor = 2, 
  nfeatures = 8
)

plotDataHeatmap(
  regressed_MOFAobject, 
  view = "CyTOF", 
  factor = 1, 
  features = 20, 
  show_rownames = T
)


plotFactorScatters(
  regressed_MOFAobject,
  factors = 1:3,
  color_by = "Volunteer",
  shape_by = "Timepoint"
)


<<<<<<< HEAD
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
=======
plotFactorBeeswarm(
  MOFAobject,
  factors = 1,
  color_by = "Timepoint"
)


r2 <- calculateVarianceExplained(regressed_MOFAobject)
r2$R2Total

# Variance explained by each factor in each view
head(r2$R2PerFactor)

# Plot it
plotVarianceExplained(MOFAobject)
>>>>>>> 5992d616634ea1357ef19c85c4d7b4f781cc9ca1

mof2 <- bind_rows(try2, mof)



<<<<<<< HEAD
library(MOFAdata)
data("CLL_data")
=======





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
  factors = c(1,2),
  color_by = "Volunteer",      # color by the IGHV values that are part of the training data
  shape_by = "Timepoint"  # shape by the trisomy12 values that are part of the training data
)


plotWeights(
  MOFAobject, 
  view = "Plasma", 
  factor = 3, 
  nfeatures = 8
)

plotDataHeatmap(
  MOFAobject, 
  view = "CyTOF", 
  factor = 1, 
  features = 20, 
  show_rownames = T
)


plotFactorScatters(
  MOFAobject,
  factors = 1:3,
  color_by = "Volunteer",
  shape_by = "Timepoint"
)









>>>>>>> f1e7fbc5cd72344b7848b02ff174a0b613b12f18




# trash code you'll proabbly never need ####

# 
# #CD4+ Flowsom
# datacd4 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_results/results/cluster_abundances.csv")
# colnames(datacd4)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')
# datacd4 <- data.frame("Feature"=paste("CD4_Cluster_", datacd4$ClusterID, sep=''), datacd4[3:20])
# rownames(datacd4) <- datacd4$Feature
# datacd4$Feature <- NULL
# datacd4 <- as.matrix(datacd4)
# 
# #CD8 flowsom
# #data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_abundances.csv")
# datacd8 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_results/results/cluster_abundances.csv")
# colnames(datacd8)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline", "DoD", "DoD+6"), times=6), sep='_')
# datacd8 <- data.frame("Feature"=paste("CD8_Cluster_", datacd8$ClusterID, sep=''), datacd8[3:20])
# rownames(datacd8) <- datacd8$Feature
# datacd8$Feature <- NULL
# datacd8 <- as.matrix(datacd8)
# 
# cytof_data <- rbind(datacd4, datacd8)
# 
# cytof_data <- cytof_data*100
# #, substr(colnames(cytof_data), 1,3), substr(colnames(cytof_data), 5,nchar(colnames(cytof_data))))
# # rownames(cytof_data)[(nrow(cytof_data)-1):nrow(cytof_data)] <- c("Volunteer", "Timepoint")

