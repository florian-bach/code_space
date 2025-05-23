library(edgeR)
library(dplyr)
library(tidyr)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)


### clean environment

remove(list = ls())

# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# map new data to same some, but have different seed sequence should produce similar, but not identical results


# read in data 

# iMac
data <- read.csv("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_02a_results/results/cluster_abundances.csv")
data2 <- read.csv("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_abundances.csv")

# laptop
# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_02a_results/results/cluster_abundances.csv")
# data2 <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_abundances.csv")


#extract number of cells in each fcs file to convert frequency to actual number
setwd("/Users/s1249052/PhD/cytof/better_gating")

# setwd("C:/Users/Florian/PhD/cytof/vac69a/T_cells_only/better_gating")


files_list <- list.files(path=".", pattern="*.fcs")

flo_set <- read.flowSet(files_list[16:20], transformation = FALSE, truncate_max_range = FALSE)


short <- data

# convert to absolute numbers: we're keeping them the same because we still need to check whether edgeR can
# adequately deal with the lib.size argument differing with the fcs file size; this is just the smallest number
# from any of the fcs files

short[,3] <- data[,3]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short[,4] <- data[,4]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short[,5] <- data[,5]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short[,6] <- data[,6]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short[,7] <- data[,7]*nrow(flo_set@frames[[files_list[19]]]@exprs)


# repeat same stuff for other dataframe 
short2 <- data2
short2[,3] <- data2[,3]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short2[,4] <- data2[,4]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short2[,5] <- data2[,5]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short2[,6] <- data2[,6]*nrow(flo_set@frames[[files_list[19]]]@exprs)
short2[,7] <- data2[,7]*nrow(flo_set@frames[[files_list[19]]]@exprs)

# combine
short <- cbind(short, short2)
rownames(short) <- short$ClusterID

short[,c(1:2,8:9)] <-NULL

# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
#groups <- rep(substr(colnames(short)[1:14],1, nchar(colnames(short))-9), times=2)

colnames(short) <- rep(c("Baseline", "C10", "C12", "DoD", "DoD6"), times=2)

colnames(short)[6:10] <- paste(colnames(short)[1:5], "_1", sep='')
groups <- rep(c("Baseline", "C10", "C12", "DoD", "DoD6"), times=2)

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)

# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short, group=groups, lib.size = colSums(short))
fit <- estimateDisp(y, design)
fit <- glmQLFit(fit, design, robust=TRUE) 


# make a contrast (matrix calculating differential expression between to or more samples) for each person
base_dod6 <- makeContrasts(groupsDoD6 - groupsBaseline, levels=design)
base_dod <- makeContrasts(groupsDoD - groupsBaseline, levels=design)
dod_dod6 <- makeContrasts(groupsDoD6 - groupsDoD, levels=design)
base_c10 <- makeContrasts(groupsC10 - groupsBaseline, levels=design)
base_c12 <- makeContrasts(groupsC12 - groupsBaseline, levels=design)


# make deg list for each person                                                                             
deg_base_dod6 <- glmQLFTest(fit, contrast=base_dod6) 
deg_base_dod <- glmQLFTest(fit, contrast=base_dod) 
deg_dod_dod6 <- glmQLFTest(fit, contrast=dod_dod6) 
deg_base_c12 <- glmQLFTest(fit, contrast=base_c12) 
deg_base_c10 <- glmQLFTest(fit, contrast=base_c10) 

# save as dataframes
deg_base_dod6 <- topTags(deg_base_dod6, n=64)$table
deg_base_dod <- topTags(deg_base_dod, n=64)$table
deg_dod_dod6 <- topTags(deg_dod_dod6, n=64)$table 
deg_base_c12 <- topTags(deg_base_c12, n=64)$table 
deg_base_c10 <- topTags(deg_base_c10, n=64)$table 


#add cluster as column 
deg_base_dod6$Cluster <- rownames(deg_base_dod6)
deg_base_dod$Cluster <- rownames(deg_base_dod)
deg_dod_dod6$Cluster <- rownames(deg_dod_dod6)
deg_base_c12$Cluster <- rownames(deg_base_c12)
deg_base_c10$Cluster <- rownames(deg_base_c10)

### this ordering is important: in the chunk later we need to match clusters from the edgeR output with the cluster abundances of FlowSOM, we decided
### to do it by indexing, so the order needs to be the same; this circumnavigates numbers v letters etc but is otherwise not very elegant, might change


deg_base_dod6 <- deg_base_dod6[order(as.numeric(deg_base_dod6$Cluster)),]
deg_base_dod <- deg_base_dod[order(as.numeric(deg_base_dod$Cluster)),]
deg_dod_dod6 <- deg_dod_dod6[order(as.numeric(deg_dod_dod6$Cluster)),]
deg_base_c12 <- deg_base_c12[order(as.numeric(deg_base_c12$Cluster)),]
deg_base_c10 <- deg_base_c10[order(as.numeric(deg_base_c10$Cluster)),]

deg_base_dod6$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), "matters_not", "matters"), "matters")
deg_base_dod$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), "matters_not", "matters"), "matters")
deg_dod_dod6$matters <- ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c12$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), ifelse(short[,2] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c10$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), ifelse(short[,3] < 0.01*nrow(flo_set@frames[[files_list[19]]]@exprs), "matters_not", "matters"), "matters")

# add volunteer, cluster column & collate deglists 
deg_base_dod6$Volunteer <- "V02"
deg_base_dod$Volunteer <- "V02"
deg_dod_dod6$Volunteer <- "V02"
deg_base_c12$Volunteer <- "V02"
deg_base_c10$Volunteer <- "V02"

deg_base_dod6$Timepoint <- "deg_base_dod6"
deg_base_dod$Timepoint <- "deg_base_dod"
deg_dod_dod6$Timepoint <- "deg_dod_dod6"
deg_base_c12$Timepoint <- "deg_base_c12"
deg_base_c10$Timepoint <- "deg_base_c10"


all_degs <- list(deg_base_dod6, deg_base_dod, deg_dod_dod6, deg_base_c12, deg_base_c10)

# make a list of lists that contain re-transformed fold change values
fold_change  <- lapply(all_degs, function(x){x$Fold_Change <- 2^x$logFC})

# add list to list of dataframes as new column; the simplify=F argument prevents conversion of every column
# to a list

all_degs <- mapply(cbind, all_degs, "Fold_Change"= fold_change, SIMPLIFY=F)

# combine them in one big file
#individual_from_all <- rbind(deg_base_dod6, deg_base_dod, deg_dod_dod6, deg_base_c8, deg_base_c10)

individual_from_all <- plyr::ldply(all_degs, rbind)
#nrow=320

individual_from_all <- dplyr::filter(individual_from_all, FDR<0.01)
#nrow=156

# subset dataframe so only fold changes over 2 and less than 0.5 are included
lower_cut_off <- dplyr::filter(individual_from_all, Fold_Change > 2)
upper_cut_off <- dplyr::filter(individual_from_all, Fold_Change < 0.5)
cut_off <- rbind(upper_cut_off, lower_cut_off)
#nrow=36

cut_off <- dplyr::filter(cut_off, matters == "matters")
nrow(cut_off)
#nrow=23

# disassemble big dataframe for making figures

list_of_degs <- split(cut_off, cut_off$Timepoint)
important_ones <- plyr::ldply(list_of_degs, rbind)
#nrow=23 ###  yaaay



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change


# iMac
setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_medians/")

# # laptop
# setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_medians/")

  
  
deg_medians_aggregate  <- read.csv("aggregate_cluster_medians.csv")

#0 1 transform of marker intensities, columnwise
deg_medians_aggregate[,5:71] <- lapply(deg_medians_aggregate[5:71], function(x){scales::rescale(x,to=c(0,1))})


### make small dataframes for each cluster comparison to then add marker expression data to it later

clusters_base_c10 <- dplyr::filter(important_ones, Timepoint=="deg_base_c10")
clusters_base_c12 <- dplyr::filter(important_ones, Timepoint=="deg_base_c12")
clusters_base_dod <- dplyr::filter(important_ones, Timepoint=="deg_base_dod")
clusters_base_dod6 <- dplyr::filter(important_ones, Timepoint=="deg_base_dod6")
clusters_dod_dod6 <- dplyr::filter(important_ones, Timepoint=="deg_dod_dod6")


### add marker expression data

medians_base_c10 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_base_c10$Cluster)  %>%
  mutate(Comparison = "base_c10") %>%
  mutate(Fold_Change = clusters_base_c10$Fold_Change)

medians_base_c12 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_base_c12$Cluster)  %>%
  mutate(Comparison = "base_c12") %>%
  mutate(Fold_Change = clusters_base_c12$Fold_Change)

medians_base_dod <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_base_dod$Cluster)  %>%
  mutate(Comparison = "base_dod") %>%
  mutate(Fold_Change = clusters_base_dod$Fold_Change)

medians_base_dod6 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_base_dod6$Cluster)  %>%
  mutate(Comparison = "base_dod6") %>%
  mutate(Fold_Change = clusters_base_dod6$Fold_Change)

medians_dod_dod6 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_dod_dod6$Cluster)  %>%
  mutate(Comparison = "dod_dod6") %>%
  mutate(Fold_Change = clusters_dod_dod6$Fold_Change)


#put it all together to make ggplots; drop UMAP channels
deg_medians_aggregate <- rbind(medians_base_c10, medians_base_c12, medians_base_dod, medians_base_dod6, medians_dod_dod6)

deg_medians_all <- select(deg_medians_aggregate, colnames(deg_medians_aggregate)[c(1, 2, 5, 16:17, 25:59, 65, 67, 72, 73)])

colnames(deg_medians_all)[4:42] <- substr(colnames(deg_medians_all)[4:42], 8, nchar(colnames(deg_medians_all)[4:42])-10)
colnames(deg_medians_all)[3] <- "CD45"


#get rid of cd45, cd3, tcrgd channels
deg_medians_all <- select(deg_medians_all, -CD45, -CD3, -CD14, -CD20, -TIM.3, -CXCR5, -CX3CR1, -CD103,  -IntegrinB7, -TCRgd, -MetaclusterID)


# order the expression datasets so that the fold change can be neatly carried over from the abundance set
deg_medians_all <- deg_medians_all[order(deg_medians_all$ClusterID),]


# make beautiful iris color palette

#my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)
#my_palette <- rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))


# make levels to reorder markers in a meaningful way

marker_levels <- c("CD4",
                   "CD8",
                   "Vd2",
                   "Va7.2",
                   "CD38",
                   "HLA.DR",
                   "ICOS",
                   "CD28",
                   "PD1",
                   "TIM.3",
                   "CD95",
                   "BCL.2",
                   "CD27",
                   "Perforin",
                   "GZB",
                   "Tbet",
                   "CTLA4",
                   "Ki.67",
                   "CD127",
                   "IntegrinB7",
                   "CD49d",
                   "CD56",
                   "CD16",
                   "CD161",
                   "CD103",
                   "CD25",
                   "FoxP3",
                   "CD39",
                   "CLA",
                   "CD14",
                   "CX3CR1",
                   "CD20",
                   "CXCR5",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")


my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")
#color_blind <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")

#######         figures for cluster abundances

#iMac
data <- read.csv("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_abundances.csv")

# laptop
# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_02b_results/results/cluster_abundances.csv")


short <- select(data, colnames(data[3:7]))

### make a dataframe for each comparison that contains the cluster abundance at the pre and post timepoint

base_dod_clusters <- short[c(clusters_base_dod$Cluster),c(4,1)]
base_dod_clusters$ClusterID <- rownames(base_dod_clusters)
base_dod_clusters$Comparison <- "base_dod"
colnames(base_dod_clusters)[1:2] <- c("post", "pre")
base_dod_clusters$Fold_Change <- base_dod_clusters$post/base_dod_clusters$pre

base_c10_clusters <- short[c(clusters_base_c10$Cluster),c(2,1)]
base_c10_clusters$ClusterID <- rownames(base_c10_clusters)
base_c10_clusters$Comparison <- "base_c10"
colnames(base_c10_clusters)[1:2] <- c("post", "pre")
base_c10_clusters$Fold_Change <- base_c10_clusters$post/base_c10_clusters$pre

base_c12_clusters <- short[c(clusters_base_c12$Cluster),c(3,1)]
base_c12_clusters$ClusterID <- rownames(base_c12_clusters)
base_c12_clusters$Comparison <- "base_c12"
colnames(base_c12_clusters)[1:2] <- c("post", "pre")
base_c12_clusters$Fold_Change <- base_c12_clusters$post/base_c12_clusters$pre

base_dod6_clusters <- short[c(as.numeric(clusters_base_dod6$Cluster)),c(5,1)]
base_dod6_clusters$ClusterID <- rownames(base_dod6_clusters)
base_dod6_clusters$Comparison <- "base_dod6"
colnames(base_dod6_clusters)[1:2] <- c("post", "pre")
base_dod6_clusters$Fold_Change <- base_dod6_clusters$post/base_dod6_clusters$pre

dod_dod6_clusters <- short[c(as.numeric(clusters_dod_dod6$Cluster)),c(5,4)]
dod_dod6_clusters$ClusterID <- rownames(dod_dod6_clusters)
dod_dod6_clusters$Comparison <- "dod_dod6"
colnames(dod_dod6_clusters)[1:2] <- c("post", "pre")
dod_dod6_clusters$Fold_Change <- dod_dod6_clusters$post/dod_dod6_clusters$pre


### put it all together

abun_clusters <- rbind(base_c10_clusters, base_c12_clusters, base_dod_clusters, base_dod6_clusters, dod_dod6_clusters)

# order the thing so the fold change carries over correctly
abun_clusters <- abun_clusters[order(as.numeric(abun_clusters$ClusterID)),]
abun_clusters$ClusterID <- as.numeric(abun_clusters$ClusterID)

long_abun_clusters <- gather(abun_clusters, Timepoint, Count, c("pre", "post"))

# add correct fold change variable
deg_medians_all$Fold_Change <- abun_clusters$Fold_Change

# long format
long_deg_medians_all <- gather(deg_medians_all, Marker, Intensity, colnames(deg_medians_all)[2:31])

# add categorical variable whether something is going up or down based on fold change
long_deg_medians_all$Direction <- ifelse(long_deg_medians_all$Fold_Change>1, "up", "down")

# reorder that variable to first the up panel is displayed in the plots
long_deg_medians_all$Directions <- factor(long_deg_medians_all$Direction, levels = c("up", "down"))

########   combine two for loops so to make a heatmap displaying median expression of each marker for each cluster; this heatmap is amde for each comparison, highlighting up and downregulated clusters. next to the heatmap will be a barplot showing cluster abundances at the pre and post timepoints of each comparison

#### barplot

for(i in unique(long_abun_clusters$Comparison)){
  
  sub_set <- dplyr::filter(long_abun_clusters, Comparison == i)
  specific_levels <- unique(sub_set[order(sub_set$Fold_Change, decreasing = TRUE),"ClusterID"])
  
  assign(paste(i,"_bar", sep=''), 
         
         ggplot(data = sub_set,
                aes_(x=factor(sub_set$ClusterID, levels = specific_levels), y=sub_set$Count, fill=factor(sub_set$Timepoint, levels=c("pre", "post")))
         )+
           geom_bar(stat="identity", position=position_dodge())+
           scale_fill_brewer(palette="Paired")+
           xlab("Cluster ID")+
           ylab("Percentage of CD3+ T cells")+
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
  
  ##### heatmap
  
  sub_set <- dplyr::filter(long_deg_medians_all, Comparison == i)
  
  assign(paste("comparison_", unique(sub_set$Comparison), sep=''),
         ggplot(data = sub_set, aes_(x=factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Comparison))+
           geom_tile(aes(fill=Intensity), color="white")+
           scale_fill_gradientn(colors=rev(my_palette))+
           scale_y_discrete(position = "left")+
           xlab("Cluster ID")+
           facet_grid(~ Directions, scales = "free")+
           ggtitle(paste("V02 ", i, sep=''))+
           theme(panel.border = element_blank(),
                 axis.text.y.left = element_text(size=35),
                 axis.line.y.left = element_blank(),
                 axis.line.y.right = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.title.x = element_text(size=33),
                 axis.text.x = element_text(size = 33),
                 axis.text.y.right = element_text(size = 35),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.title = element_blank(),
                 legend.position = "none",
                 plot.title = element_text(size = 45, hjust = 0.5),
                 plot.margin = unit(c(1,0,1,0), "cm"),
                 strip.text.x = element_text(size=28))
  )
}

setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("v02_01_heatmap_plus_abundance_base_dod.pdf", grid.arrange(comparison_base_dod, base_dod_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("v02_01_heatmap_plus_abundance_base_c10.pdf", grid.arrange(comparison_base_c10, base_c10_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("v02_01_heatmap_plus_abundance_base_c12.pdf", grid.arrange(comparison_base_c12, base_c12_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("v02_01_heatmap_plus_abundance_base_dod6.pdf", grid.arrange(comparison_base_dod6, base_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("v02_01_heatmap_plus_abundance_dod_dod6.pdf", grid.arrange(comparison_dod_dod6, dod_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)



mat <- spread(long_deg_medians_all, Marker, Intensity)

mat2 <- mat %>%
  dplyr::filter(., Comparison=="dod_dod6") %>%
  dplyr::filter(., Direction=="up")%>%
  dplyr::select(., -Comparison, -Fold_Change, -Direction, -Directions)

mat2 <- as.matrix(mat2)

tmat <- t(mat2)

colnames(tmat) <- tmat[1,]

corr_mat=cor(tmat[2:nrow(tmat),],method="s")

# figure only with only what's up from dod to dod6

##### heatmap

sub_set <- dplyr::filter(long_deg_medians_all, Comparison == "dod_dod6")
sub_set <- dplyr::filter(sub_set, Direction=="up")

sub_set$ClusterID <- as.character(sub_set$ClusterID)

specific_levels <- rownames(corr_mat[order(corr_mat[,1], decreasing = TRUE),])

up_v2 <- ggplot(data = sub_set, aes_(x=factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Comparison))+
  geom_tile(aes(fill=Intensity), color="white")+
  scale_fill_gradientn(colors=rev(my_palette))+
  scale_y_discrete(position = "left")+
  xlab("Cluster ID")+
  ggtitle("Volunteer 02")+
  theme(panel.border = element_blank(),
        axis.text.y.left = element_text(size=18),
        axis.line.y.left = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=18),
        axis.text.x = element_text(size = 18),
        axis.text.y.right = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"),
        strip.text.x = element_text(size=16))

# setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures")

ggsave("v02_up.pdf", height=7.5, width=5)
ggsave("v02_up.png", height=7.5, width=5)



#imac
# setwd("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/figures")

# laptop
setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")













up_dod_dod6 <- filter(long_abun_clusters, Comparison=="dod_dod6")
up_dod_dod6 <- filter(up_dod_dod6, Fold_Change>1)
up_dod_dod6 <- filter(up_dod_dod6, Timepoint=="post")

up_dod_dod6[,c(2,4)] <- NULL
up_dod_dod6$Count <- up_dod_dod6$Count*100
up_dod_dod6 <- up_dod_dod6[order(up_dod_dod6$Fold_Change, decreasing = T),]

up_dod_dod6$SubSet <- c("CD4", "CD4", "MAIT", "CD4", "CD8")

up_dod_dod6 <- up_dod_dod6[order(up_dod_dod6$SubSet),]
#up_dod_dod6$ClusterID <- paste0("Cluster", up_dod_dod6$ClusterID, sep=" \n")

up_dod_dod6$ymin[1]<-0
up_dod_dod6$ymin[2:nrow(up_dod_dod6)] <- cumsum(up_dod_dod6$Count)[1:nrow(up_dod_dod6)-1]
up_dod_dod6$ymax <- up_dod_dod6$ymin + up_dod_dod6$Count


ggplot(up_dod_dod6)+
  
  geom_rect(aes(fill=Fold_Change, ymin=up_dod_dod6$ymin, ymax=up_dod_dod6$ymax, xmax=2.8, xmin=2, colour=Fold_Change))+
  geom_rect(aes(fill=as.numeric(factor(ClusterID))*100, ymin=ymin, ymax=ymax, xmax=6, xmin=3))+
  geom_rect(aes(fill=as.numeric(factor(SubSet))*500, ymin=ymin, ymax=ymax, xmax=9, xmin=6.2))+
  scale_color_gradientn(colours=rev(sequential_hcl(5, "Heat")))+
  scale_fill_gradientn(guide = FALSE, colors=my_long_palette, values=c(scales::rescale(seq(1,50), to=c(0,0.8)),scales::rescale(c(500,900,1000,1500,1800,2000,2500), to=c(0.8,1))))+
  
  geom_text(aes(x=4.5, y=(ymin+ymax)/2, label = paste("Cluster", ClusterID, sep="\n")),size=2.4,fontface="bold")+
  geom_text(aes(x=7.5, y=(ymin+ymax)/2, label = paste(SubSet)), size=2.7, fontface="bold")+
  
  guides(color=guide_colorbar(barwidth = 2))+
  labs(color = "Fold Change")+
  theme(aspect.ratio=1,
        legend.title = element_text(),
        legend.title.align = 0.5,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  xlim(c(0, 9))+
  coord_polar("y")


setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
up_dod_dod6$Volunteer <- "V02"
write.csv(up_dod_dod6, "v02_up_dod_dod6.csv")





