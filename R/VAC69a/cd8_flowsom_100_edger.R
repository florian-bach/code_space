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

data <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_results/results/cluster_abundances.csv")
data2 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_(copy)_(copy)_results/results/cluster_abundances.csv")

colnames(data)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline_01", "DoD_01", "DoD+6_01"), times=6), sep='_')
colnames(data2)[3:20] <- paste(rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3), rep(c("Baseline_02", "DoD_02", "DoD+6_02"), times=6), sep='_')


# get rid of uneccessary columns, convert fractions to counts

#cd8 3,659

short<-data
short2<-data2


number_of_cells <- 3659

short[3:20] <- short[3:20]*number_of_cells # this is the number of cells from each fcs file
short2[3:20] <- short2[3:20]*number_of_cells # this is the number of cells from each fcs file

# combine
short <- cbind(short, short2[3:20])

# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
groups <- rep(
  paste(
    rep(c("V06", "V07", "V09", "V02", "V03", "V05"), each=3),
    rep(c("Baseline", "DoD", "DoD_6"), times=6), sep='_'), times=2)

#groups <- gsub("T.6", "DoD.6", groups, fixed=T) 

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)

colnames(design) <- substr(colnames(design), 7, 20)

# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short[,3:38], group=groups, lib.size = colSums(short[,3:38]))
fit <- estimateDisp(y, design)
fit <- glmQLFit(fit, design, robust=TRUE) 

# make a contrast (matrix calculating differential expression between to or more samples) for each person
pre_post_02 <- makeContrasts(V02_DoD_6 - V02_Baseline,levels=design)
pre_post_03 <- makeContrasts(V03_DoD_6 - V03_Baseline,levels=design)
pre_post_05 <- makeContrasts(V05_DoD_6 - V05_Baseline,levels=design)
pre_post_06 <- makeContrasts(V06_DoD_6 - V06_Baseline,levels=design)
pre_post_07 <- makeContrasts(V07_DoD_6 - V07_Baseline,levels=design)
pre_post_09 <- makeContrasts(V09_DoD_6 - V09_Baseline,levels=design)

# make deg list for each person                                                                             
deg_02 <- glmQLFTest(fit, contrast=pre_post_02) 
deg_03 <- glmQLFTest(fit, contrast=pre_post_03) 
deg_05 <- glmQLFTest(fit, contrast=pre_post_05) 
deg_06 <- glmQLFTest(fit, contrast=pre_post_06) 
deg_07 <- glmQLFTest(fit, contrast=pre_post_07) 
deg_09 <- glmQLFTest(fit, contrast=pre_post_09) 

# save as dataframes
deg_02 <- topTags(deg_02, n=100)$table
deg_03 <- topTags(deg_03, n=100)$table
deg_05 <- topTags(deg_05, n=100)$table 
deg_06 <- topTags(deg_06, n=100)$table 
deg_07 <- topTags(deg_07, n=100)$table 
deg_09 <- topTags(deg_09, n=100)$table

#add cluser as column 
deg_02$Cluster <- rownames(deg_02)
deg_03$Cluster <- rownames(deg_03)
deg_05$Cluster <- rownames(deg_05)
deg_06$Cluster <- rownames(deg_06)
deg_07$Cluster <- rownames(deg_07)
deg_09$Cluster <- rownames(deg_09)

deg_02 <- deg_02[order(as.numeric(deg_02$Cluster)),]; deg_02$Baseline <- short[,'V02_Baseline_01']; deg_02$Treatment <- short[,'V02_DoD+6_01']
deg_03 <- deg_03[order(as.numeric(deg_03$Cluster)),]; deg_03$Baseline <- short[,'V03_Baseline_01']; deg_03$Treatment <- short[,'V03_DoD+6_01']
deg_05 <- deg_05[order(as.numeric(deg_05$Cluster)),]; deg_05$Baseline <- short[,'V05_Baseline_01']; deg_05$Treatment <- short[,'V05_DoD+6_01']
deg_06 <- deg_06[order(as.numeric(deg_06$Cluster)),]; deg_06$Baseline <- short[,'V06_Baseline_01']; deg_06$Treatment <- short[,'V06_DoD+6_01']
deg_07 <- deg_07[order(as.numeric(deg_07$Cluster)),]; deg_07$Baseline <- short[,'V07_Baseline_01']; deg_07$Treatment <- short[,'V07_DoD+6_01']
deg_09 <- deg_09[order(as.numeric(deg_09$Cluster)),]; deg_09$Baseline <- short[,'V09_Baseline_01']; deg_09$Treatment <- short[,'V09_DoD+6_01']

#mark stuff that never surpasses 1% frequency
deg_02$matters <- ifelse(deg_02$Baseline < 0.01*number_of_cells, ifelse(deg_02$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_03$matters <- ifelse(deg_03$Baseline < 0.01*number_of_cells, ifelse(deg_03$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_05$matters <- ifelse(deg_05$Baseline < 0.01*number_of_cells, ifelse(deg_05$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_06$matters <- ifelse(deg_06$Baseline < 0.01*number_of_cells, ifelse(deg_06$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_07$matters <- ifelse(deg_07$Baseline < 0.01*number_of_cells, ifelse(deg_07$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")
deg_09$matters <- ifelse(deg_09$Baseline < 0.01*number_of_cells, ifelse(deg_09$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")



# add volunteer, cluster column & collate deglists 
deg_02$Volunteer <- "V_02"
deg_03$Volunteer <- "V_03"
deg_05$Volunteer <- "V_05"
deg_06$Volunteer <- "V_06"
deg_07$Volunteer <- "V_07"
deg_09$Volunteer <- "V_09"

# deg_02$Cluster <- rownames(deg_02)
# deg_03$Cluster <- rownames(deg_03)
# deg_05$Cluster <- rownames(deg_05)
# deg_06$Cluster <- rownames(deg_06)
# deg_07$Cluster <- rownames(deg_07)
# deg_09$Cluster <- rownames(deg_09)

# combine them in one big file
individual_from_all <- rbind(deg_02, deg_03, deg_05, deg_06, deg_07, deg_09)

# subset dataframe so only fold changes over 2 and less than 0.5 are included
upper_cut_off <- filter(individual_from_all, logFC > 1)
lower_cut_off <- filter(individual_from_all, logFC < -1)
cut_off <- rbind(upper_cut_off, lower_cut_off)
cut_off <- filter(cut_off, matters == "matters")
nrow(cut_off)
# 61

cut_off$Direction <- ifelse(cut_off$logFC>1, "up", "down")


# disassemble big dataframe for making figures


list_of_degs <- split(cut_off, cut_off$Volunteer)



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change

# laptop
# setwd("C:/Users/Florian/PhD/cytof/vac69a/big_flowsoms/FlowSOM_all_cd4s_baseline_dod_t6_(copy)_(copy)_results/results/cluster_medians")

setwd("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_results/results/cluster_medians")

median_02 <- read.csv("12c05_sample1_01_0_0_cluster_medians.csv")
median_03 <- read.csv("12c10_sample1_01_0_0_cluster_medians.csv")
median_05 <- read.csv("12c15_sample1_01_0_0_cluster_medians.csv")
median_06 <- read.csv("14c05_sample1_01_0_0_cluster_medians.csv")
median_07 <- read.csv("14c10_sample1_01_0_0_cluster_medians.csv")
median_09 <- read.csv("14c15_sample1_01_0_0_cluster_medians.csv")

# super convoluted way of doing it but it works: restrict cluster medians to clusters that can be found in deg list with cutoff of log2 1 and -1
deg_medians_02 <- median_02 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "02") %>%
  mutate(Direction = list_of_degs[[1]][order(as.numeric(list_of_degs[[1]]$Cluster)),]$Direction)

deg_medians_03 <- median_03 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[2]][order(as.numeric(list_of_degs[[2]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "03") %>%
  mutate(Direction = list_of_degs[[2]][order(as.numeric(list_of_degs[[2]]$Cluster)),]$Direction)


deg_medians_05 <- median_05 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[3]][order(as.numeric(list_of_degs[[3]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "05") %>%
  mutate(Direction = list_of_degs[[3]][order(as.numeric(list_of_degs[[3]]$Cluster)),]$Direction)


deg_medians_06 <- median_06 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[4]][order(as.numeric(list_of_degs[[4]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "06") %>%
  mutate(Direction = list_of_degs[[4]][order(as.numeric(list_of_degs[[4]]$Cluster)),]$Direction)



deg_medians_07 <- median_07 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[5]][order(as.numeric(list_of_degs[[5]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "07") %>%
  mutate(Direction = list_of_degs[[5]][order(as.numeric(list_of_degs[[5]]$Cluster)),]$Direction)


deg_medians_09 <- median_09 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[6]][order(as.numeric(list_of_degs[[6]]$Cluster)),]$Cluster)) %>%
  mutate(Volunteer = "09") %>%
  mutate(Direction = list_of_degs[[6]][order(as.numeric(list_of_degs[[6]]$Cluster)),]$Direction)

#put it all together to make ggplots; drop UMAP channels
deg_medians_all <- rbind(deg_medians_02, deg_medians_03, deg_medians_05, deg_medians_06, deg_medians_07, deg_medians_09)
deg_medians_all$ClusterID <- as.character(deg_medians_all$ClusterID)

data_medians_all <- select(deg_medians_all, colnames(deg_medians_all)[c(1, 5, 17, 25:59, 65, 67, 72, 73)])


# this bit spikes in the lowest and hightest value for each channel taken from a flowsom run on all T cells
# in order to adapt the notion of positiviy away from a z score specific to cd4s or cd8s
spike <- read.csv("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_06b_results/results/cluster_medians/aggregate_cluster_medians.csv")
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

data_medians_all <- dplyr::select(data_medians_all, -CD45, -CD8, -CD4, -Vd2, -CD3, -CD20, -TIM.3, -CXCR5, -CX3CR1, -CD103,  -IntegrinB7, -ICOS, -TCRgd, -FoxP3, -CD25, -CD39, -CLA, -CD16)





data_medians_all[,2:(ncol(data_medians_all)-2)] <- lapply(data_medians_all[,2:(ncol(data_medians_all)-2)], function(x){scales::rescale(x,to=c(0,1))})







colnames(data_medians_all) <- gsub(".", "-", colnames(data_medians_all), fixed=T)
colnames(data_medians_all) <- gsub("Va7-2", "Va7.2", colnames(data_medians_all), fixed=T)

long_deg_medians_all <- gather(data_medians_all, Marker, Intensity, colnames(data_medians_all)[2:(ncol(data_medians_all)-2)])

# make beautiful iris color palette

#my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)
#my_palette <- rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))


# make levels to reorder markers in a meaningful way

marker_levels <- c("CD4",
"CD8",
"Va7.2",
"Vd2",
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
"CTLA4",
"Ki-67",
"CD127",
"IntegrinB7",
"CD56",
"CD16",
"CD161",
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

# arcsinh transform
# (big_heatmap_arsinh <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
#     geom_tile(aes(fill=Intensity), color="white")+
#     scale_fill_gradientn(colors=rev(my_palette))+
#     xlab("Cluster")+
#     ylab("Marker")+
#     ggtitle("/"Differentially Expressed\" Cluster Medians Per Volunteer in VAC69a")+
#     theme(#axis.text.x = element_text(hjust = 1),
#           panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           legend.title = element_blank(),
#           plot.title = element_text(hjust = 0.5))+
#     facet_wrap(~ Volunteer, ncol=6, labeller=as_labeller(labs))
# )
# 
# # arcsinh transform, accesible color scheme
# (big_heatmap_arsinh <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
#     geom_tile(aes(fill=Intensity), color="white")+
#     scale_fill_gradientn(colors=color_blind)+
#     xlab("Cluster")+
#     ylab("Marker")+
#     ggtitle("\"Differentially Expressed\" Cluster Medians Per Volunteer in VAC69a")+
#     theme(#axis.text.x = element_text(hjust = 1),
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.line = element_line(colour = "black"),
#       legend.title = element_blank(),
#       plot.title = element_text(hjust = 0.5))+
#     facet_wrap(~ Volunteer, ncol=5, labeller=as_labeller(labs))
# )
# 
# 

# unique(long_deg_medians_all$Volunteer)


# just makin a table







# figure only with only what's up from dod to dod6
mat <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_results/results/aggregate_cluster_CVs.csv")
mat <- select(mat, colnames(mat)[c(1, 5, 17, 25:59, 65, 67)])


aggregate_medians <- select(deg_medians_all, colnames(deg_medians_all)[c(1, 5, 17, 25:59, 65, 67, 72)])

colnames(mat) <- substr(colnames(aggregate_medians)[1:40], 8, nchar(colnames(aggregate_medians)[1:40])-10)
colnames(mat)[c(1,2)] <- c("ClusterID", "CD45")


# convert to long format

mat2 <- dplyr::select(mat, -CD45, -CD8, -CD4, -Va7.2, -Vd2, -CD3, -CD20, -TIM.3, -CXCR5, -CX3CR1, -CD103,  -IntegrinB7, -TCRgd, -CTLA4)


rownames(mat2) <- mat2$ClusterID

mat2 <- as.matrix(mat2)
tmat <- t(mat2[,2:ncol(mat2)])

corr_mat=cor(tmat, method="s")

specific_levels <- rownames(corr_mat[order(corr_mat[,2], decreasing = T),])





short <-  select(data, colnames(data[3:20]))
# 

clusters_02 <- as.integer(list_of_degs[[1]]$Cluster)
clusters_03 <- as.integer(list_of_degs[[2]]$Cluster)
clusters_05 <- as.integer(list_of_degs[[3]]$Cluster)
clusters_06 <- as.integer(list_of_degs[[4]]$Cluster)
clusters_07 <- as.integer(list_of_degs[[5]]$Cluster)
clusters_09 <- as.integer(list_of_degs[[6]]$Cluster)

abun_clusters_02 <- short[c(clusters_02),c(10,12)]
abun_clusters_02$ClusterID <- as.character(clusters_02)
abun_clusters_02$Volunteer <- substr(colnames(abun_clusters_02)[2], 2, 3) 
colnames(abun_clusters_02) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_03 <- short[c(clusters_03),c(13,15)]
abun_clusters_03$ClusterID <- as.character(clusters_03)
abun_clusters_03$Volunteer <- substr(colnames(abun_clusters_03)[2], 2, 3) 
colnames(abun_clusters_03) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_05 <- short[c(clusters_05),c(16,18)]
abun_clusters_05$ClusterID <- as.character(clusters_05)
abun_clusters_05$Volunteer <- substr(colnames(abun_clusters_05)[2], 2, 3) 
colnames(abun_clusters_05) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_06 <- short[c(clusters_06),c(1,3)]
abun_clusters_06$ClusterID <- as.character(clusters_06)
abun_clusters_06$Volunteer <- substr(colnames(abun_clusters_06)[2], 2, 3) 
colnames(abun_clusters_06) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_07 <- short[c(clusters_07),c(4,6)]
abun_clusters_07$ClusterID <- as.character(clusters_07)
abun_clusters_07$Volunteer <- substr(colnames(abun_clusters_07)[2], 2, 3) 
colnames(abun_clusters_07) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_09 <- short[c(clusters_09),c(7,9)]
abun_clusters_09$ClusterID <- as.character(clusters_09)
abun_clusters_09$Volunteer <- substr(colnames(abun_clusters_09)[2], 2, 3) 
colnames(abun_clusters_09) <- c("Baseline", "T+6", "ClusterID", "Volunteer")


abun_clusters <- rbind(abun_clusters_02,abun_clusters_03,abun_clusters_05,abun_clusters_06,abun_clusters_07,abun_clusters_09)
long_abun_clusters <- gather(abun_clusters, Timepoint, Frequency, c("Baseline", "T+6"))






##############          working figure



for(i in unique(long_deg_medians_all$Volunteer)){
  ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
  
  sub_set <- filter(long_deg_medians_all, Volunteer == i)
  sub_set <- filter(sub_set, Direction=="up")
  
  assign(paste("Volunteer_", unique(sub_set$Volunteer), sep=''),
         ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Volunteer))+
           geom_tile(aes(fill=Intensity), color="white")+
           scale_fill_gradientn(colors=rev(my_palette))+
           scale_y_discrete(position = result1)+
           xlab(NULL)+
           theme(panel.border = element_blank(),
                 axis.text.y.left = result,
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
                aes_(x=factor(sub_set2$ClusterID, levels = specific_levels), y=sub_set2$Frequency, fill=factor(sub_set2$Timepoint, levels=c("Baseline", "T+6")))
         )+
           geom_bar(stat="identity", position=position_dodge())+
           scale_fill_brewer(palette="Paired")+
           
           ylab("% of CD8+ T cells")+
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
                 axis.title.y = result,
                 #axis.title.y = element_text(size=20, color="black"),
                 axis.text.y = element_text(size=20, color="black")))
} 



# laptop
# setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
# iMac
setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("cd8_01_d6.png", grid.arrange(Volunteer_02, Volunteer_03, Volunteer_05, Volunteer_06, Volunteer_07, Volunteer_09, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                             c(4,5,6))
),  width = 40, height = 40, limitsize = F)



ggsave("v02_cd8_heatmaps_barplots.pdf", plot_grid(V02_bar, Volunteer_02, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v03_cd8_heatmaps_barplots.pdf", plot_grid(V03_bar, Volunteer_03, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v05_cd8_heatmaps_barplots.pdf", plot_grid(V05_bar, Volunteer_05, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v06_cd8_heatmaps_barplots.pdf", plot_grid(V06_bar, Volunteer_06, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v07_cd8_heatmaps_barplots.pdf", plot_grid(V07_bar, Volunteer_07, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v09_cd8_heatmaps_barplots.pdf", plot_grid(V09_bar, Volunteer_09, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)


ggsave("v02_cd8_heatmaps_barplots.png", plot_grid(V02_bar, Volunteer_02, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v03_cd8_heatmaps_barplots.png", plot_grid(V03_bar, Volunteer_03, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v05_cd8_heatmaps_barplots.png", plot_grid(V05_bar, Volunteer_05, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v06_cd8_heatmaps_barplots.png", plot_grid(V06_bar, Volunteer_06, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v07_cd8_heatmaps_barplots.png", plot_grid(V07_bar, Volunteer_07, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)
ggsave("v09_cd8_heatmaps_barplots.png", plot_grid(V09_bar, Volunteer_09, ncol=1, rel_heights = c(1,2), align="v"), height = 17, width=11.3)








#######         figures for cluster abundances

data <- read.csv("/Users/s1249052/PhD/cytof/better_gating/big_flowsoms/FlowSOM_all_cd8s_baseline_dod_t6_better_(copy)_(copy)_results/results/cluster_abundances.csv")
short <- select(data, colnames(data[3:20]))


clusters_02 <- as.integer(list_of_degs[[1]]$Cluster)
clusters_03 <- as.integer(list_of_degs[[2]]$Cluster)
clusters_05 <- as.integer(list_of_degs[[3]]$Cluster)
clusters_06 <- as.integer(list_of_degs[[4]]$Cluster)
clusters_07 <- as.integer(list_of_degs[[5]]$Cluster)
clusters_09 <- as.integer(list_of_degs[[6]]$Cluster)

abun_clusters_02 <- short[c(clusters_02),c(10,12)]
abun_clusters_02$ClusterID <- as.factor(clusters_02)
abun_clusters_02$Volunteer <- "V02"
colnames(abun_clusters_02) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_03 <- short[c(clusters_03),c(13,15)]
abun_clusters_03$ClusterID <- as.factor(clusters_03)
abun_clusters_03$Volunteer <- "V03"
colnames(abun_clusters_03) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_05 <- short[c(clusters_05),c(16,18)]
abun_clusters_05$ClusterID <- as.factor(clusters_05)
abun_clusters_05$Volunteer <- "V05"
colnames(abun_clusters_05) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_06 <- short[c(clusters_06),c(1,3)]
abun_clusters_06$ClusterID <- as.factor(clusters_06)
abun_clusters_06$Volunteer <- "V06"
colnames(abun_clusters_06) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_07 <- short[c(clusters_07),c(4,6)]
abun_clusters_07$ClusterID <- as.factor(clusters_07)
abun_clusters_07$Volunteer <- "V07"
colnames(abun_clusters_07) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_09 <- short[c(clusters_09),c(7,9)]
abun_clusters_09$ClusterID <- as.factor(clusters_09)
abun_clusters_09$Volunteer <- "V09"
colnames(abun_clusters_09) <- c("Baseline", "T+6", "ClusterID", "Volunteer")


abun_clusters <- rbind(abun_clusters_02,abun_clusters_03,abun_clusters_05,abun_clusters_06,abun_clusters_07,abun_clusters_09)
long_abun_clusters <- gather(abun_clusters, Timepoint, Frequency, c("Baseline", "T+6"))


for(i in unique(long_abun_clusters$Volunteer)){
  
(graphic <- 
  
  ggplot(data=filter(long_abun_clusters, Volunteer == "02"),
         aes(x=factor(ClusterID, levels = Volunteer_02_levels), y=Frequency, fill=Timepoint)
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
        axis.text.y = element_text(size=20, color="black")))

  Volunteer_02_bar <- graphic 
}

ggsave("heatmap_plus_abundance_02.pdf", grid.arrange(Volunteer_02, Volunteer_02_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_03.pdf", grid.arrange(Volunteer_03, Volunteer_03_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_05.pdf", grid.arrange(Volunteer_05, Volunteer_05_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)

ggsave("heatmap_plus_abundance_06.pdf", grid.arrange(Volunteer_06, Volunteer_06_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_07.pdf", grid.arrange(Volunteer_07, Volunteer_07_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_09.pdf", grid.arrange(Volunteer_09, Volunteer_09_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)

