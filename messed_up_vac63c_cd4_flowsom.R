library(edgeR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)


# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# read in data (laptop)
# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_abundances.csv")
# data2 <-read.csv("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_results/results/cluster_abundances.csv")

# read in data (iMac)

data <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_results/results/cluster_abundances.csv")
data2 <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_2_results/results/cluster_abundances.csv")


colnames(data)[3:49] <- substr(colnames(data)[3:49], nchar(colnames(data)[3:49])-10, nchar(colnames(data)[3:49])-4)
colnames(data2)[3:49] <- substr(colnames(data2)[3:49], nchar(colnames(data2)[3:49])-10, nchar(colnames(data2)[3:49])-4)


# get rid of uneccessary columns, convert fractions to counts

#cd8 3,659

short<-data
short2<-data2

# get rid of control files
short <- short[,-grep("tr", colnames(short), fixed = T)]
short2 <- short2[,-grep("tr", colnames(short2), fixed = T)]


number_of_cells <- 2306

short[3:46] <- short[3:46]*number_of_cells # this is the number of cells from each fcs file
short2[3:46] <- short2[3:46]*number_of_cells # this is the number of cells from each fcs file

# combine
short <- cbind(short, short2[3:46])


#clean up column
colnames(short) <- gsub("_", "", colnames(short), fixed=T)
colnames(short) <- gsub(".", "", colnames(short), fixed=T)

# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
groups <- colnames(short)[3:ncol(short)]

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)

colnames(design) <- substr(colnames(design), 7, 20)
colnames(design) <- make.names(colnames(design))


# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short[,3:ncol(short)], group=groups, lib.size = colSums(short[,3:ncol(short)]))
fit <- estimateDisp(y, design)
fit <- glmQLFit(fit, design, robust=TRUE) 


# make a contrast (matrix calculating differential expression between to or more samples) for each person

# volunteers <- as.character(c(302, 303, 304, 305, 306, 307, 308, 310, 313, 315, 318, 320))

# c-1 to c+45
pre_post_301 <- makeContrasts(X301C45  - X301C1, levels=design)
pre_post_302 <- makeContrasts(X302C45  - X302C1, levels=design)
pre_post_304 <- makeContrasts(X304C45  - X304C1, levels=design)
pre_post_305 <- makeContrasts(X305C45  - X305C1, levels=design)
pre_post_306 <- makeContrasts(X306C45  - X306C1, levels=design)
pre_post_307 <- makeContrasts(X307C45  - X307C1, levels=design)
pre_post_308 <- makeContrasts(X308C45  - X308C1, levels=design)
pre_post_310 <- makeContrasts(X310C45  - X310C1, levels=design)
pre_post_313 <- makeContrasts(X313C45  - X313C1, levels=design)
pre_post_315 <- makeContrasts(X315C45  - X315C1, levels=design)
pre_post_320 <- makeContrasts(X320C45  - X320C1, levels=design)

# c-1 until dod
pre_dod_301 <- makeContrasts(X301DoD  - X301C1, levels=design)
pre_dod_302 <- makeContrasts(X302DoD  - X302C1, levels=design)
pre_dod_304 <- makeContrasts(X304DoD  - X304C1, levels=design)
pre_dod_305 <- makeContrasts(X305DoD  - X305C1, levels=design)
pre_dod_306 <- makeContrasts(X306DoD  - X306C1, levels=design)
pre_dod_307 <- makeContrasts(X307DoD  - X307C1, levels=design)
pre_dod_308 <- makeContrasts(X308DoD  - X308C1, levels=design)
pre_dod_310 <- makeContrasts(X310DoD  - X310C1, levels=design)
pre_dod_313 <- makeContrasts(X313DoD  - X313C1, levels=design)
pre_dod_315 <- makeContrasts(X315Dod  - X315C1, levels=design)
pre_dod_320 <- makeContrasts(X320DoD  - X320C1, levels=design)

# dod until t+6
dod_t6_301 <- makeContrasts(X301T6  - X301DoD, levels=design)
dod_t6_302 <- makeContrasts(X302T6  - X302DoD, levels=design)
dod_t6_304 <- makeContrasts(X304T6  - X304DoD, levels=design)
dod_t6_305 <- makeContrasts(X305T6  - X305DoD, levels=design)
dod_t6_306 <- makeContrasts(X306T6  - X306DoD, levels=design)
dod_t6_307 <- makeContrasts(X307T6  - X307DoD, levels=design)
dod_t6_308 <- makeContrasts(X308T6  - X308DoD, levels=design)
dod_t6_310 <- makeContrasts(X310T6  - X310DoD, levels=design)
dod_t6_313 <- makeContrasts(X313T6  - X313DoD, levels=design)
dod_t6_315 <- makeContrasts(X315T6  - X315Dod, levels=design)
dod_t6_320 <- makeContrasts(X320T6  - X320DoD, levels=design)


# make deg list for each person & comparison                                                                            
deg_301 <- glmQLFTest(fit, contrast=dod_t6_301) 
deg_302 <- glmQLFTest(fit, contrast=dod_t6_302) 
deg_304 <- glmQLFTest(fit, contrast=dod_t6_304) 
deg_305 <- glmQLFTest(fit, contrast=dod_t6_305) 
deg_306 <- glmQLFTest(fit, contrast=dod_t6_306) 
deg_307 <- glmQLFTest(fit, contrast=dod_t6_307) 
deg_308 <- glmQLFTest(fit, contrast=dod_t6_308) 
deg_310 <- glmQLFTest(fit, contrast=dod_t6_310) 
deg_313 <- glmQLFTest(fit, contrast=dod_t6_313) 
deg_315 <- glmQLFTest(fit, contrast=dod_t6_315) 
deg_320 <- glmQLFTest(fit, contrast=dod_t6_320) 

# save as dataframes

deg_301 <- topTags(deg_301, n=100)$table
deg_302 <- topTags(deg_302, n=100)$table
deg_304 <- topTags(deg_304, n=100)$table
deg_305 <- topTags(deg_305, n=100)$table
deg_306 <- topTags(deg_306, n=100)$table
deg_307 <- topTags(deg_307, n=100)$table
deg_308 <- topTags(deg_308, n=100)$table
deg_310 <- topTags(deg_310, n=100)$table
deg_313 <- topTags(deg_313, n=100)$table
deg_315 <- topTags(deg_315, n=100)$table
deg_320 <- topTags(deg_320, n=100)$table

#add cluser as column 

deg_301$Cluster <- rownames(deg_301)
deg_302$Cluster <- rownames(deg_302)
deg_304$Cluster <- rownames(deg_304)
deg_305$Cluster <- rownames(deg_305)
deg_306$Cluster <- rownames(deg_306)
deg_307$Cluster <- rownames(deg_307)
deg_308$Cluster <- rownames(deg_308)
deg_310$Cluster <- rownames(deg_310)
deg_313$Cluster <- rownames(deg_313)
deg_315$Cluster <- rownames(deg_315)
deg_320$Cluster <- rownames(deg_320)


# deg_301 <- deg_301[order(as.numeric(deg_301$Cluster)),]; deg_301$Baseline <- short[,'301C1']; deg_301$Treatment <- short[,'301C45']
# deg_302 <- deg_302[order(as.numeric(deg_302$Cluster)),]; deg_302$Baseline <- short[,'302C1']; deg_302$Treatment <- short[,'302C45']
# deg_304 <- deg_304[order(as.numeric(deg_304$Cluster)),]; deg_304$Baseline <- short[,'304C1']; deg_304$Treatment <- short[,'304C45']
# deg_305 <- deg_305[order(as.numeric(deg_305$Cluster)),]; deg_305$Baseline <- short[,'305C1']; deg_305$Treatment <- short[,'305C45']
# deg_306 <- deg_306[order(as.numeric(deg_306$Cluster)),]; deg_306$Baseline <- short[,'306C1']; deg_306$Treatment <- short[,'306C45']
# deg_307 <- deg_307[order(as.numeric(deg_307$Cluster)),]; deg_307$Baseline <- short[,'307C1']; deg_307$Treatment <- short[,'307C45']
# deg_308 <- deg_308[order(as.numeric(deg_308$Cluster)),]; deg_308$Baseline <- short[,'308C1']; deg_308$Treatment <- short[,'308C45']
# deg_310 <- deg_310[order(as.numeric(deg_310$Cluster)),]; deg_310$Baseline <- short[,'310C1']; deg_310$Treatment <- short[,'310C45']
# deg_313 <- deg_313[order(as.numeric(deg_313$Cluster)),]; deg_313$Baseline <- short[,'313C1']; deg_313$Treatment <- short[,'313C45']
# deg_315 <- deg_315[order(as.numeric(deg_315$Cluster)),]; deg_315$Baseline <- short[,'315C1']; deg_315$Treatment <- short[,'315C45']
# deg_320 <- deg_320[order(as.numeric(deg_320$Cluster)),]; deg_320$Baseline <- short[,'320C1']; deg_320$Treatment <- short[,'320C45']



deg_301 <- deg_301[order(as.numeric(deg_301$Cluster)),]; deg_301$pre <- short[,'301DoD']; deg_301$post <- short[,'301T6']
deg_302 <- deg_302[order(as.numeric(deg_302$Cluster)),]; deg_302$pre <- short[,'302DoD']; deg_302$post <- short[,'302T6']
deg_304 <- deg_304[order(as.numeric(deg_304$Cluster)),]; deg_304$pre <- short[,'304DoD']; deg_304$post <- short[,'304T6']
deg_305 <- deg_305[order(as.numeric(deg_305$Cluster)),]; deg_305$pre <- short[,'305DoD']; deg_305$post <- short[,'305T6']
deg_306 <- deg_306[order(as.numeric(deg_306$Cluster)),]; deg_306$pre <- short[,'306DoD']; deg_306$post <- short[,'306T6']
deg_307 <- deg_307[order(as.numeric(deg_307$Cluster)),]; deg_307$pre <- short[,'307DoD']; deg_307$post <- short[,'307T6']
deg_308 <- deg_308[order(as.numeric(deg_308$Cluster)),]; deg_308$pre <- short[,'308DoD']; deg_308$post <- short[,'308T6']
deg_310 <- deg_310[order(as.numeric(deg_310$Cluster)),]; deg_310$pre <- short[,'310DoD']; deg_310$post <- short[,'310T6']
deg_313 <- deg_313[order(as.numeric(deg_313$Cluster)),]; deg_313$pre <- short[,'313DoD']; deg_313$post <- short[,'313T6']
deg_315 <- deg_315[order(as.numeric(deg_315$Cluster)),]; deg_315$pre <- short[,'315Dod']; deg_315$post <- short[,'315T6']
deg_320 <- deg_320[order(as.numeric(deg_320$Cluster)),]; deg_320$pre <- short[,'320DoD']; deg_320$post <- short[,'320T6']





#mark stuff that never surpasses 1% frequency
deg_301$matters <- ifelse(deg_301$pre < 0.01*number_of_cells, ifelse(deg_301$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_302$matters <- ifelse(deg_302$pre < 0.01*number_of_cells, ifelse(deg_302$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_304$matters <- ifelse(deg_304$pre < 0.01*number_of_cells, ifelse(deg_304$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_305$matters <- ifelse(deg_305$pre < 0.01*number_of_cells, ifelse(deg_305$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_306$matters <- ifelse(deg_306$pre < 0.01*number_of_cells, ifelse(deg_306$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_307$matters <- ifelse(deg_307$pre < 0.01*number_of_cells, ifelse(deg_307$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_308$matters <- ifelse(deg_308$pre < 0.01*number_of_cells, ifelse(deg_308$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_310$matters <- ifelse(deg_310$pre < 0.01*number_of_cells, ifelse(deg_310$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_313$matters <- ifelse(deg_313$pre < 0.01*number_of_cells, ifelse(deg_313$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_315$matters <- ifelse(deg_315$pre < 0.01*number_of_cells, ifelse(deg_315$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_320$matters <- ifelse(deg_320$pre < 0.01*number_of_cells, ifelse(deg_320$post < 0.01*number_of_cells, "matters_not", "matters"), "matters")



# add volunteer, cluster column & collate deglists 
deg_301$Volunteer <- "V_301"
deg_302$Volunteer <- "V_302"
deg_304$Volunteer <- "V_304"
deg_305$Volunteer <- "V_305"
deg_306$Volunteer <- "V_306"
deg_307$Volunteer <- "V_307"
deg_308$Volunteer <- "V_308"
deg_310$Volunteer <- "V_310"
deg_313$Volunteer <- "V_313"
deg_315$Volunteer <- "V_315"
deg_320$Volunteer <- "V_320"

# deg_02$Cluster <- rownames(deg_02)
# deg_03$Cluster <- rownames(deg_03)
# deg_05$Cluster <- rownames(deg_05)
# deg_06$Cluster <- rownames(deg_06)
# deg_07$Cluster <- rownames(deg_07)
# deg_09$Cluster <- rownames(deg_09)

# combine them in one big file
individual_from_all <- rbind(deg_301, deg_302, deg_304, deg_305, deg_306, deg_307, deg_308, deg_310, deg_313, deg_315, deg_320)

# subset dataframe so only fold changes over 2 and less than 0.5 are included
upper_cut_off <- dplyr::filter(individual_from_all, logFC > 1)
lower_cut_off <- dplyr::filter(individual_from_all, logFC < -1)
cut_off <- rbind(upper_cut_off, lower_cut_off)
cut_off <- dplyr::filter(cut_off, matters == "matters")
cut_off$Direction <- ifelse(cut_off$logFC>1, "up", "down")
nrow(cut_off)
#131

# disassemble big dataframe for making figures


list_of_degs <- split(cut_off, cut_off$Volunteer)



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change

# laptop
# setwd("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_medians")

setwd("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_results/results/cluster_medians")

list_of_files <- list.files()

median_301 <- read.csv(list_of_files[max(grep("301", list_of_files))])
median_302 <- read.csv(list_of_files[max(grep("302", list_of_files))])
median_304 <- read.csv(list_of_files[max(grep("304", list_of_files))])
median_305 <- read.csv(list_of_files[max(grep("305", list_of_files))])
median_306 <- read.csv(list_of_files[max(grep("306", list_of_files))])
median_307 <- read.csv(list_of_files[max(grep("307", list_of_files))])
median_308 <- read.csv(list_of_files[max(grep("308", list_of_files))])
median_310 <- read.csv(list_of_files[max(grep("310", list_of_files))])
median_313 <- read.csv(list_of_files[max(grep("313", list_of_files))])
median_315 <- read.csv(list_of_files[max(grep("315", list_of_files))])
median_320 <- read.csv(list_of_files[max(grep("320", list_of_files))])


# super convoluted way of doing it but it works: restrict cluster medians to clusters that can be found in deg list with cutoff of log2 1 and -1


deg_medians_301 <- median_301 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "301") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)

deg_medians_302 <- median_302 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[2]][order(as.numeric(list_of_degs[[2]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "302") %>%
  mutate(Direction = list_of_degs[[2]][order(as.numeric(list_of_degs[[2]]$Cluster)),]$Direction)

deg_medians_304 <- median_304 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[3]][order(as.numeric(list_of_degs[[3]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "304") %>%
  mutate(Direction = list_of_degs[[3]][order(as.numeric(list_of_degs[[3]]$Cluster)),]$Direction)

deg_medians_305 <- median_305 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[4]][order(as.numeric(list_of_degs[[4]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "305") %>%
  mutate(Direction = list_of_degs[[4]][order(as.numeric(list_of_degs[[4]]$Cluster)),]$Direction)

deg_medians_306 <- median_306 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[5]][order(as.numeric(list_of_degs[[5]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "306") %>%
  mutate(Direction = list_of_degs[[5]][order(as.numeric(list_of_degs[[5]]$Cluster)),]$Direction)

deg_medians_307 <- median_307 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[6]][order(as.numeric(list_of_degs[[6]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "307") %>%
  mutate(Direction = list_of_degs[[6]][order(as.numeric(list_of_degs[[6]]$Cluster)),]$Direction)

deg_medians_308 <- median_308 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[7]][order(as.numeric(list_of_degs[[7]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "308") %>%
  mutate(Direction = list_of_degs[[7]][order(as.numeric(list_of_degs[[7]]$Cluster)),]$Direction)

deg_medians_310 <- median_310 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[8]][order(as.numeric(list_of_degs[[8]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "310") %>%
  mutate(Direction = list_of_degs[[8]][order(as.numeric(list_of_degs[[8]]$Cluster)),]$Direction)

deg_medians_313 <- median_313 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[9]][order(as.numeric(list_of_degs[[9]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "313") %>%
  mutate(Direction = list_of_degs[[9]][order(as.numeric(list_of_degs[[9]]$Cluster)),]$Direction)

deg_medians_315 <- median_315 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[10]][order(as.numeric(list_of_degs[[10]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "315") %>%
  mutate(Direction = list_of_degs[[10]][order(as.numeric(list_of_degs[[10]]$Cluster)),]$Direction)

deg_medians_320 <- median_320 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[11]][order(as.numeric(list_of_degs[[11]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "320") %>%
  mutate(Direction = list_of_degs[[11]][order(as.numeric(list_of_degs[[11]]$Cluster)),]$Direction)


#put it all together to make ggplots; drop UMAP channels
deg_medians_all <- rbind(deg_medians_301, deg_medians_302, deg_medians_304, deg_medians_305, deg_medians_306, deg_medians_307, deg_medians_308, deg_medians_310, deg_medians_313, deg_medians_315, deg_medians_320)
deg_medians_all$ClusterID <- as.character(deg_medians_all$ClusterID)

data_medians_all <- select(deg_medians_all, colnames(deg_medians_all)[c(1, 5, 17, 25:59, 65, 67, 72, 73)])

# this bit spikes in the lowest and hightest value for each channel taken from a flowsom run on all T cells
# in order to adapt the notion of positiviy away from a z score specific to cd4s or cd8s

spike <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_results/results/cluster_medians/aggregate_cluster_medians.csv")
spike <- select(spike, colnames(spike)[c(1, 5, 17, 25:59, 65, 67)])

spike_col_min <- unlist(lapply(spike, min))
spike_col_max <- unlist(lapply(spike, max))

spike_min_max <- data.frame(rbind(spike_col_min, spike_col_max))
spike_min_max$ClusterID <- 0
spike_min_max <- mutate(spike_min_max, Volunteer=0, Direction=0) 


data_medians_all <- rbind(data_medians_all, spike_min_max)


#make column names marker names
# colnames(data_medians_all) <- c('ClusterID', 'MetaclusterID', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
#                           '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
#                           '156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
#                           '165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
#                           '174Yb_CLA', '175Lu_Perforin', '198Pt_CD8', '209Bi_CD16', 'Volunteer', 'Cluster (Fold Change)')

colnames(data_medians_all) <- substr(colnames(data_medians_all), 8, nchar(colnames(data_medians_all))-10)

# convert to long format
colnames(data_medians_all)[c(1,2,41,42)] <- c("ClusterID", "CD45", "Volunteer", "Direction")

data_medians_all <- dplyr::select(data_medians_all, -CD45, -CD8, -CD4, -Va72, -Vd2, -CD3, -TIM3, -CXCR5, -CX3CR1, -TCRgd)

data_medians_all[,2:(ncol(data_medians_all)-2)] <- lapply(data_medians_all[,2:(ncol(data_medians_all)-2)], function(x){scales::rescale(x,to=c(0,1))})

colnames(data_medians_all) <- gsub(".", "-", colnames(data_medians_all), fixed=T)

long_deg_medians_all <- tidyr::gather(data_medians_all, Marker, Intensity, colnames(data_medians_all)[2:(ncol(data_medians_all)-2)])

# make beautiful iris color palette

#my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)
#my_palette <- rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))


# make levels to reorder markers in a meaningful way

marker_levels <- c("CD4",
                   "CD8",
                   "Vd2",
                   "CD69",
                   "CD38",
                   "HLADR",
                   "ICOS",
                   "CD28",
                   "PD1",
                   "TIM3",
                   "CD95",
                   "BCL2",
                   "CD27",
                   "Perforin",
                   "GZB",
                   "Tbet",
                   "GATA3",
                   "Eomes",
                   "CTLA4",
                   "Ki67",
                   "CD127",
                   "CD56",
                   "CD16",
                   "CD161",
                   "RORgt",
                   "CD49d",
                   "CD103",
                   "CD25",
                   "FoxP3",
                   "CD39",
                   "CLA",
                   "CXCR5",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")
# make some heatmaps bruv



my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")
#color_blind <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")




# figure only with only what's up from dod to dod6

# order clusters by correlation (starting point = handpicked for high activation level)
mat <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_results/results/aggregate_cluster_CVs.csv")
mat <- select(mat, colnames(mat)[c(1, 5, 17, 25:59, 65, 67)])


aggregate_medians <- select(deg_medians_all, colnames(deg_medians_all)[c(1, 5, 17, 25:59, 65, 67, 72)])

colnames(mat) <- substr(colnames(aggregate_medians)[1:40], 8, nchar(colnames(aggregate_medians)[1:40])-10)
colnames(mat)[c(1,2)] <- c("ClusterID", "CD45")


# convert to long format

mat2 <- dplyr::select(mat, -CD45, -CD8, -CD4, -Va72, -Vd2, -CD3, -TIM3, -CXCR5, -CX3CR1, -TCRgd)


rownames(mat2) <- mat2$ClusterID

mat2 <- as.matrix(mat2)
tmat <- t(mat2[,2:ncol(mat2)])

corr_mat=cor(tmat, method="s")

specific_levels <- rownames(corr_mat[order(corr_mat[,19], decreasing = T),])




#######         figures for cluster abundances

data <- read.csv("/Users/s1249052/PhD/cytof/vac63c/analysis/FlowSOM_all_cd4+_results/results/cluster_abundances.csv")
short <-  select(data, colnames(data[3:49]))
# 

clusters_301 <- as.integer(list_of_degs[[1]]$Cluster)
clusters_302 <- as.integer(list_of_degs[[2]]$Cluster)
clusters_304 <- as.integer(list_of_degs[[3]]$Cluster)
clusters_305 <- as.integer(list_of_degs[[4]]$Cluster)
clusters_306 <- as.integer(list_of_degs[[5]]$Cluster)
clusters_307 <- as.integer(list_of_degs[[6]]$Cluster)
clusters_308 <- as.integer(list_of_degs[[7]]$Cluster)
clusters_310 <- as.integer(list_of_degs[[8]]$Cluster)
clusters_313 <- as.integer(list_of_degs[[9]]$Cluster)
clusters_315 <- as.integer(list_of_degs[[10]]$Cluster)
clusters_320 <- as.integer(list_of_degs[[11]]$Cluster)


#########     change the numbers at the end to pick different timepoints: 1=baseline, 2=dod, 3=t6, 4=c+45


abun_clusters_301 <- short[c(clusters_301),c(grep(301, colnames(short))[c(2,3)])]
abun_clusters_301$ClusterID <- as.character(clusters_301)
abun_clusters_301$Volunteer <- "301"
colnames(abun_clusters_301) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_302 <- short[c(clusters_302),c(grep(302, colnames(short))[c(2,3)])]
abun_clusters_302$ClusterID <- as.character(clusters_302)
abun_clusters_302$Volunteer <- "302"
colnames(abun_clusters_302) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_304 <- short[c(clusters_304),c(grep(304, colnames(short))[c(2,3)])]
abun_clusters_304$ClusterID <- as.character(clusters_304)
abun_clusters_304$Volunteer <- "304"
colnames(abun_clusters_304) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_305 <- short[c(clusters_305),c(grep(305, colnames(short))[c(2,3)])]
abun_clusters_305$ClusterID <- as.character(clusters_305)
abun_clusters_305$Volunteer <- "305"
colnames(abun_clusters_305) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_306 <- short[c(clusters_306),c(grep(306, colnames(short))[c(2,3)])]
abun_clusters_306$ClusterID <- as.character(clusters_306)
abun_clusters_306$Volunteer <- "306"
colnames(abun_clusters_306) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_307 <- short[c(clusters_307),c(grep(307, colnames(short))[c(2,3)])]
abun_clusters_307$ClusterID <- as.character(clusters_307)
abun_clusters_307$Volunteer <- "307"
colnames(abun_clusters_307) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_308 <- short[c(clusters_308),c(grep(308, colnames(short))[c(2,3)])]
abun_clusters_308$ClusterID <- as.character(clusters_308)
abun_clusters_308$Volunteer <- "308"
colnames(abun_clusters_308) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_310 <- short[c(clusters_310),c(grep(310, colnames(short))[c(2,3)])]
abun_clusters_310$ClusterID <- as.character(clusters_310)
abun_clusters_310$Volunteer <- "310"
colnames(abun_clusters_310) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_313 <- short[c(clusters_313),c(grep(313, colnames(short))[c(2,3)])]
abun_clusters_313$ClusterID <- as.character(clusters_313)
abun_clusters_313$Volunteer <- "313"
colnames(abun_clusters_313) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_315 <- short[c(clusters_315),c(grep(315, colnames(short))[c(2,3)])]
abun_clusters_315$ClusterID <- as.character(clusters_315)
abun_clusters_315$Volunteer <- "315"
colnames(abun_clusters_315) <- c("pre", "post", "ClusterID", "Volunteer")

abun_clusters_320 <- short[c(clusters_320),c(grep(320, colnames(short))[c(2,3)])]
abun_clusters_320$ClusterID <- as.character(clusters_320)
abun_clusters_320$Volunteer <- "320"
colnames(abun_clusters_320) <- c("pre", "post", "ClusterID", "Volunteer")



abun_clusters <- rbind(abun_clusters_301, abun_clusters_302, abun_clusters_304, abun_clusters_305, abun_clusters_306, abun_clusters_307, abun_clusters_308, abun_clusters_310, abun_clusters_313, abun_clusters_315, abun_clusters_320)
long_abun_clusters <- tidyr::gather(abun_clusters, Timepoint, Frequency, c("pre", "post"))



##############          working figure


for(i in unique(long_deg_medians_all$Volunteer)){
  # ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  # ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
  # 
  sub_set <- dplyr::filter(long_deg_medians_all, Volunteer == i)
  sub_set <- filter(sub_set, Direction=="up")
  
  assign(paste("Volunteer_", unique(sub_set$Volunteer), sep=''),
         # ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Volunteer))+
         ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Volunteer))+
           geom_tile(aes(fill=Intensity), color="white")+
           scale_fill_gradientn(colors=rev(my_palette))+
           # scale_y_discrete(position = result1)+
           xlab(NULL)+
           theme(panel.border = element_blank(),
                 # axis.text.y.left = result,
                 axis.line.y.left = element_blank(),
                 axis.line.y.right = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_text(size = 33),
                 axis.text.y.right = element_text(size = 35),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.title = element_blank(),
                 legend.position = "none",
                 plot.margin = unit(c(1,0,1,0), "cm"))
  )
  
  
  sub_set2 <- dplyr::filter(long_abun_clusters, Volunteer == i)
  sub_set2 <- dplyr::filter(sub_set2, ClusterID %in% sub_set$ClusterID)
  
  assign(paste("V", i,"_bar", sep=''), 
         
         ggplot(data = sub_set2,
                aes_(x=factor(sub_set2$ClusterID, levels = specific_levels), y=sub_set2$Frequency, fill=factor(sub_set2$Timepoint, levels=c("pre", "post")))
         )+
           geom_bar(stat="identity", position=position_dodge())+
           scale_fill_brewer(palette="Paired")+
           
           ylab("% of CD4+ T cells")+
           scale_y_continuous(position= "left", labels = scales::percent_format())+
           ggtitle(paste("Volunteer ", i, "\n", sep=''))+
           theme(legend.title = element_blank(),
                 legend.text = element_text(size = 20),
                 legend.position = "top", 
                 legend.justification = "center",
                 legend.direction = "horizontal",
                 axis.line = element_line(colour = "black"),
                 plot.title = element_text(size = 45, hjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 # axis.title.y = result,
                 #axis.title.y = element_text(size=20, color="black"),
                 axis.text.y = element_text(size=20, color="black")))
} 


# laptop
# setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
# iMac
setwd("/Users/s1249052/PhD/cytof/vac63c/figures")



plot_grid(V301_bar, Volunteer_301, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V302_bar, Volunteer_302, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V304_bar, Volunteer_304, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V305_bar, Volunteer_305, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V306_bar, Volunteer_306, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V307_bar, Volunteer_307, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V308_bar, Volunteer_308, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V310_bar, Volunteer_310, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V313_bar, Volunteer_313, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V315_bar, Volunteer_315, ncol=1, rel_heights = c(1,2), align="v")
plot_grid(V320_bar, Volunteer_320, ncol=1, rel_heights = c(1,2), align="v")

