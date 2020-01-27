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


number_of_cells <- 2306

short[3:49] <- short[3:49]*number_of_cells # this is the number of cells from each fcs file
short2[3:49] <- short2[3:49]*number_of_cells # this is the number of cells from each fcs file

# combine
short <- cbind(short, short2[3:49])
colnames(short) <- gsub("_", "", colnames(short), fixed=T)
colnames(short) <- gsub(".", "", colnames(short), fixed=T)

# get rid of control files
short[,grep("tr", colnames(short), fixed = T)] <- NULL

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
deg_301 <- glmQLFTest(fit, contrast=pre_post_301) 
deg_302 <- glmQLFTest(fit, contrast=pre_post_302) 
deg_304 <- glmQLFTest(fit, contrast=pre_post_304) 
deg_305 <- glmQLFTest(fit, contrast=pre_post_305) 
deg_306 <- glmQLFTest(fit, contrast=pre_post_306) 
deg_307 <- glmQLFTest(fit, contrast=pre_post_307) 
deg_308 <- glmQLFTest(fit, contrast=pre_post_308) 
deg_310 <- glmQLFTest(fit, contrast=pre_post_310) 
deg_313 <- glmQLFTest(fit, contrast=pre_post_313) 
deg_315 <- glmQLFTest(fit, contrast=pre_post_315) 
deg_320 <- glmQLFTest(fit, contrast=pre_post_320) 

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


deg_301 <- deg_301[order(as.numeric(deg_301$Cluster)),]; deg_301$Baseline <- short[,'302C1']; deg_301$Treatment <- short[,'302C45']
deg_302 <- deg_302[order(as.numeric(deg_302$Cluster)),]; deg_302$Baseline <- short[,'302C1']; deg_302$Treatment <- short[,'302C45']
deg_304 <- deg_304[order(as.numeric(deg_304$Cluster)),]; deg_304$Baseline <- short[,'302C1']; deg_304$Treatment <- short[,'302C45']
deg_305 <- deg_305[order(as.numeric(deg_305$Cluster)),]; deg_305$Baseline <- short[,'302C1']; deg_305$Treatment <- short[,'302C45']
deg_306 <- deg_306[order(as.numeric(deg_306$Cluster)),]; deg_306$Baseline <- short[,'302C1']; deg_306$Treatment <- short[,'302C45']
deg_307 <- deg_307[order(as.numeric(deg_307$Cluster)),]; deg_307$Baseline <- short[,'302C1']; deg_307$Treatment <- short[,'302C45']
deg_308 <- deg_308[order(as.numeric(deg_308$Cluster)),]; deg_308$Baseline <- short[,'302C1']; deg_308$Treatment <- short[,'302C45']
deg_310 <- deg_310[order(as.numeric(deg_310$Cluster)),]; deg_310$Baseline <- short[,'302C1']; deg_310$Treatment <- short[,'302C45']
deg_313 <- deg_313[order(as.numeric(deg_313$Cluster)),]; deg_313$Baseline <- short[,'302C1']; deg_313$Treatment <- short[,'302C45']
deg_315 <- deg_315[order(as.numeric(deg_315$Cluster)),]; deg_315$Baseline <- short[,'302C1']; deg_315$Treatment <- short[,'302C45']
deg_320 <- deg_320[order(as.numeric(deg_320$Cluster)),]; deg_320$Baseline <- short[,'302C1']; deg_320$Treatment <- short[,'302C45']


#mark stuff that never surpasses 1% frequency
deg_301$matters <- ifelse(deg_301$Baseline < 0.01*number_of_cells, ifelse(deg_301$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_302$matters <- ifelse(deg_302$Baseline < 0.01*number_of_cells, ifelse(deg_302$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_304$matters <- ifelse(deg_304$Baseline < 0.01*number_of_cells, ifelse(deg_304$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_305$matters <- ifelse(deg_305$Baseline < 0.01*number_of_cells, ifelse(deg_305$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_306$matters <- ifelse(deg_306$Baseline < 0.01*number_of_cells, ifelse(deg_306$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_307$matters <- ifelse(deg_307$Baseline < 0.01*number_of_cells, ifelse(deg_307$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_308$matters <- ifelse(deg_308$Baseline < 0.01*number_of_cells, ifelse(deg_308$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_310$matters <- ifelse(deg_310$Baseline < 0.01*number_of_cells, ifelse(deg_310$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_313$matters <- ifelse(deg_313$Baseline < 0.01*number_of_cells, ifelse(deg_313$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_315$matters <- ifelse(deg_315$Baseline < 0.01*number_of_cells, ifelse(deg_315$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_320$matters <- ifelse(deg_320$Baseline < 0.01*number_of_cells, ifelse(deg_320$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")



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
nrow(cut_off)
cut_off$Direction <- ifelse(cut_off$logFC>1, "up", "down")
#273

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
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "302") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_304 <- median_304 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "304") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_305 <- median_305 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "305") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_306 <- median_306 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "306") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_307 <- median_307 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "307") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_308 <- median_308 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "308") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_310 <- median_310 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "310") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_313 <- median_313 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "313") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_315 <- median_315 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "315") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)
deg_medians_320 <- median_320 %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "320") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)


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
                   "HLA-DR",
                   "ICOS",
                   "CD28",
                   "PD1",
                   "TIM-3",
                   "CD95",
                   "BCL-2",
                   "CD27",
                   "Perforin",
                   "GZB",
                   "Tbet",
                   "GATA3",
                   "Eomes",
                   "CTLA4",
                   "Ki-67",
                   "CD127",
                   "IntegrinB7",
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


abun_clusters_301 <- short[c(clusters_301),c(grep(301, colnames(short))[c(1,4)])]
abun_clusters_301$ClusterID <- as.character(clusters_301)
abun_clusters_301$Volunteer <- "301"
colnames(abun_clusters_301) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_302 <- short[c(clusters_302),c(grep(302, colnames(short))[c(1,4)])]
abun_clusters_302$ClusterID <- as.character(clusters_302)
abun_clusters_302$Volunteer <- "302"
colnames(abun_clusters_302) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_304 <- short[c(clusters_304),c(grep(304, colnames(short))[c(1,4)])]
abun_clusters_304$ClusterID <- as.character(clusters_304)
abun_clusters_304$Volunteer <- "304"
colnames(abun_clusters_304) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_305 <- short[c(clusters_305),c(grep(305, colnames(short))[c(1,4)])]
abun_clusters_305$ClusterID <- as.character(clusters_305)
abun_clusters_305$Volunteer <- "305"
colnames(abun_clusters_305) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_306 <- short[c(clusters_306),c(grep(306, colnames(short))[c(1,4)])]
abun_clusters_306$ClusterID <- as.character(clusters_306)
abun_clusters_306$Volunteer <- "306"
colnames(abun_clusters_306) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_307 <- short[c(clusters_307),c(grep(307, colnames(short))[c(1,4)])]
abun_clusters_307$ClusterID <- as.character(clusters_307)
abun_clusters_307$Volunteer <- "307"
colnames(abun_clusters_307) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_308 <- short[c(clusters_308),c(grep(308, colnames(short))[c(1,4)])]
abun_clusters_308$ClusterID <- as.character(clusters_308)
abun_clusters_308$Volunteer <- "308"
colnames(abun_clusters_308) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_310 <- short[c(clusters_310),c(grep(310, colnames(short))[c(1,4)])]
abun_clusters_310$ClusterID <- as.character(clusters_310)
abun_clusters_310$Volunteer <- "310"
colnames(abun_clusters_310) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_313 <- short[c(clusters_313),c(grep(313, colnames(short))[c(1,4)])]
abun_clusters_313$ClusterID <- as.character(clusters_313)
abun_clusters_313$Volunteer <- "313"
colnames(abun_clusters_313) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_315 <- short[c(clusters_315),c(grep(315, colnames(short))[c(1,4)])]
abun_clusters_315$ClusterID <- as.character(clusters_315)
abun_clusters_315$Volunteer <- "315"
colnames(abun_clusters_315) <- c("Baseline", "C+45", "ClusterID", "Volunteer")

abun_clusters_320 <- short[c(clusters_320),c(grep(320, colnames(short))[c(1,4)])]
abun_clusters_320$ClusterID <- as.character(clusters_320)
abun_clusters_320$Volunteer <- "320"
colnames(abun_clusters_320) <- c("Baseline", "C+45", "ClusterID", "Volunteer")



abun_clusters <- rbind(abun_clusters_301, abun_clusters_302, abun_clusters_304, abun_clusters_305, abun_clusters_306, abun_clusters_307, abun_clusters_308, abun_clusters_310, abun_clusters_313, abun_clusters_315, abun_clusters_320)
long_abun_clusters <- tidyr::gather(abun_clusters, Timepoint, Frequency, c("Baseline", "C+45"))



##############          working figure


for(i in unique(long_deg_medians_all$Volunteer)){
  # ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  # ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
  # 
  sub_set <- dplyr::filter(long_deg_medians_all, Volunteer == i)
  # sub_set <- filter(sub_set, Direction=="up")
  
  assign(paste("Volunteer_", unique(sub_set$Volunteer), sep=''),
         ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Volunteer))+
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
                aes_(x=factor(sub_set2$ClusterID, levels = specific_levels), y=sub_set2$Frequency, fill=factor(sub_set2$Timepoint, levels=c("Baseline", "C+45")))
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
setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("cd4_01_d6.png", grid.arrange(Volunteer_02, Volunteer_03, Volunteer_05, Volunteer_06, Volunteer_07, Volunteer_09, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                               c(4,5,6))
),  width = 40, height = 40, limitsize = F)

setwd("/Users/s1249052/PhD/cytof/vac63c/figures")
ggsave("cd4_baseline_c45.png", grid.arrange(Volunteer_301, Volunteer_302, Volunteer_304, Volunteer_305, Volunteer_306, Volunteer_307, Volunteer_308, Volunteer_310, Volunteer_313, Volunteer_315, Volunteer_320, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                               c(4,5,6),
                                                                                                                                                               c(7, 8, 9),
                                                                                                                                                               c(10, 11))
),  width = 40, height = 40, limitsize = F)

grid.arrange(V02_bar , V03_bar, V05_bar, V06_bar, V07_bar, V09_bar)

ggsave("cd4_01_d6.png", grid.arrange(Volunteer_02, V02_bar, Volunteer_03, V03_bar,Volunteer_05, V05_bar,Volunteer_06, V06_bar,Volunteer_07, V07_bar, Volunteer_09, V09_bar), ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                                                                                   c(4,5,6))
),  width = 40, height = 40, limitsize = F)




ggsave("v02_heatmaps_barplots.pdf", plot_grid(V02_bar, Volunteer_02, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v301_heatmaps_barplots.pdf", plot_grid(V301_bar, Volunteer_301, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v05_heatmaps_barplots.pdf", plot_grid(V05_bar, Volunteer_05, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v06_heatmaps_barplots.pdf", plot_grid(V06_bar, Volunteer_06, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v07_heatmaps_barplots.pdf", plot_grid(V07_bar, Volunteer_07, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v09_heatmaps_barplots.pdf", plot_grid(V09_bar, Volunteer_09, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)


ggsave("v02_heatmaps_barplots.png", plot_grid(V02_bar, Volunteer_02, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v03_heatmaps_barplots.png", plot_grid(V03_bar, Volunteer_03, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v05_heatmaps_barplots.png", plot_grid(V05_bar, Volunteer_05, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v06_heatmaps_barplots.png", plot_grid(V06_bar, Volunteer_06, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v07_heatmaps_barplots.png", plot_grid(V07_bar, Volunteer_07, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v09_heatmaps_barplots.png", plot_grid(V09_bar, Volunteer_09, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)


for(i in unique(long_abun_clusters$Volunteer)){
  
  (assign(paste(i, "_bar", sep=''), ggplot(data=filter(long_abun_clusters, Volunteer == i),
                                           aes(x=factor(ClusterID, levels = specific_levels), y=Frequency, fill=Timepoint)
  )+
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    xlab("Cluster ID")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.position = "top", 
          legend.justification = "center",
          legend.direction = "horizontal",
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=20, color="black"),
          axis.title.x = element_text(size=24, color="black"),
          axis.title.y = element_text(size=24, color="black"),
          axis.text.y = element_text(size=20, color="black"))))
  
  
}

ggsave("heatmap_plus_abundance_02.pdf", grid.arrange(Volunteer_02, Volunteer_02_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_03.pdf", grid.arrange(Volunteer_03, Volunteer_03_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_05.pdf", grid.arrange(Volunteer_05, Volunteer_05_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)

ggsave("heatmap_plus_abundance_06.pdf", grid.arrange(Volunteer_06, Volunteer_06_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_07.pdf", grid.arrange(Volunteer_07, Volunteer_07_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_09.pdf", grid.arrange(Volunteer_09, Volunteer_09_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)



ggsave("indie_up_d6.png", grid.arrange(up_v2, up_v3, up_v5, up_v6, up_v7, up_v9, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),c(4,5,6))),width = 40, height = 40, limitsize = F)

theme_boy <- theme_bw()+theme(panel.border = element_blank(),
                              axis.line.y.left = element_blank(),
                              axis.line.y.right = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text.x = element_text(size = 33),
                              axis.text.y.left = element_text(size = 35),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              legend.title = element_blank(),
                              legend.position = "none",
                              plot.title = element_text(size = 45, hjust = 0.5),
                              plot.margin = unit(c(1,0,1,0), "cm"))

ggsave("indie_up_d6.png", grid.arrange(up_v2+theme_boy, up_v3+theme_boy, up_v5+theme_boy, up_v6+theme_boy, up_v7+theme_boy, up_v9+theme_boy, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                                                   c(4,5,6))
),  width = 40, height = 40, limitsize = F)
