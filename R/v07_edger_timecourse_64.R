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


# read in data 
data <- read.csv("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_abundances.csv")
data2 <- read.csv("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_07b_results/results/cluster_abundances.csv")

# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_abundances.csv")
# data2 <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_07b_results/results/cluster_abundances.csv")




  

#extract number of cells in each fcs file to convert frequency to actual number
setwd("/Users/s1249052//PhD/cytof/better_gating")
# setwd("C:/Users/Florian/PhD/cytof/vac69a/T_cells_only/better_gating")

files_list <- list.files(path=".", pattern="*.fcs")

flo_set <- read.flowSet(files_list[6:10], transformation = FALSE, truncate_max_range = FALSE)


short <- data

# convert to absolute numbers: we're keeping them the same because we still need to check whether edgeR can
# adequately deal with the lib.size argument differing with the fcs file sizel this is just the smallest number
# from any of the fcs files

short[,3] <- data[,3]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short[,4] <- data[,4]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short[,5] <- data[,5]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short[,6] <- data[,6]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short[,7] <- data[,7]*nrow(flo_set@frames[[files_list[9]]]@exprs)

# short[,3] <- data[,3]*nrow(flo_set@frames[[files_list[21]]]@exprs)
# short[,4] <- data[,4]*nrow(flo_set@frames[[files_list[22]]]@exprs)
# short[,5] <- data[,5]*nrow(flo_set@frames[[files_list[23]]]@exprs)
# short[,6] <- data[,6]*nrow(flo_set@frames[[files_list[9]]]@exprs)
# short[,7] <- data[,7]*nrow(flo_set@frames[[files_list[25]]]@exprs)


# repeat same stuff for other dataframe 
short2 <- data2
short2[,3] <- data2[,3]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short2[,4] <- data2[,4]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short2[,5] <- data2[,5]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short2[,6] <- data2[,6]*nrow(flo_set@frames[[files_list[9]]]@exprs)
short2[,7] <- data2[,7]*nrow(flo_set@frames[[files_list[9]]]@exprs)
# combine
short <- cbind(short, short2)
rownames(short) <- short$ClusterID
short[,c(1:2,8:9)] <-NULL
# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
#groups <- rep(substr(colnames(short)[1:14],1, nchar(colnames(short))-9), times=2)

colnames(short) <- rep(c("Baseline", "C8", "C10", "DoD", "DoD6"), times=2)
colnames(short)[6:10] <- paste(colnames(short)[1:5], "_1", sep='')
groups <- rep(c("Baseline", "C8", "C10", "DoD", "DoD6"), times=2)

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)

# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short, group=groups, lib.size = colSums(short))
fit <- estimateDisp(y, design)
fit <- glmQLFit(fit, design, robust=TRUE) 

# (optional: look for stuff that's differentially expressed across the board)
# pre_post_all <- makeContrasts(T.6_02 - baseline_02,
#                            T.6_03 - baseline_03,
#                            T.6_05 - baseline_05,
#                            T.6_07 - baseline_07,
#                            T.6_09 - baseline_09,
#                            levels=design)
# dge_all <- glmQLFTest(fit, contrast=pre_post_all) 
# topTags(dge_all)

# make a contrast (matrix calculating differential expression between to or more samples) for each person
base_dod6 <- makeContrasts(groupsDoD6 - groupsBaseline, levels=design)
base_dod <- makeContrasts(groupsDoD - groupsBaseline, levels=design)
dod_dod6 <- makeContrasts(groupsDoD6 - groupsDoD, levels=design)
base_c8 <- makeContrasts(groupsC8 - groupsBaseline, levels=design)
base_c10 <- makeContrasts(groupsC10 - groupsBaseline, levels=design)


# make deg list for each person                                                                             
deg_base_dod6 <- glmQLFTest(fit, contrast=base_dod6) 
deg_base_dod <- glmQLFTest(fit, contrast=base_dod) 
deg_dod_dod6 <- glmQLFTest(fit, contrast=dod_dod6) 
deg_base_c8 <- glmQLFTest(fit, contrast=base_c8) 
deg_base_c10 <- glmQLFTest(fit, contrast=base_c10) 

# save as dataframes
deg_base_dod6 <- topTags(deg_base_dod6, n=64)$table
deg_base_dod <- topTags(deg_base_dod, n=64)$table
deg_dod_dod6 <- topTags(deg_dod_dod6, n=64)$table 
deg_base_c8 <- topTags(deg_base_c8, n=64)$table 
deg_base_c10 <- topTags(deg_base_c10, n=64)$table 


#add cluser as column 
deg_base_dod6$Cluster <- rownames(deg_base_dod6)
deg_base_dod$Cluster <- rownames(deg_base_dod)
deg_dod_dod6$Cluster <- rownames(deg_dod_dod6)
deg_base_c8$Cluster <- rownames(deg_base_c8)
deg_base_c10$Cluster <- rownames(deg_base_c10)

# deg_02 <- deg_02[order(as.numeric(deg_02$Cluster)),]; deg_02$Baseline <- short[,7]; deg_02$Treatment <- short[,1]
# deg_03 <- deg_03[order(as.numeric(deg_03$Cluster)),]; deg_03$Baseline <- short[,8]; deg_03$Treatment <- short[,2]
# deg_05 <- deg_05[order(as.numeric(deg_05$Cluster)),]; deg_05$Baseline <- short[,9]; deg_05$Treatment <- short[,3]
# deg_06 <- deg_06[order(as.numeric(deg_06$Cluster)),]; deg_06$Baseline <- short[,10]; deg_06$Treatment <- short[,4]
# deg_07 <- deg_07[order(as.numeric(deg_07$Cluster)),]; deg_07$Baseline <- short[,11]; deg_07$Treatment <- short[,5]
# deg_09 <- deg_09[order(as.numeric(deg_09$Cluster)),]; deg_09$Baseline <- short[,12]; deg_09$Treatment <- short[,6]

#mark stuff that never surpasses 1% frequency
# deg_02$matters <- ifelse(deg_02$Baseline < 0.01*number_of_cells, ifelse(deg_02$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")


### this ordering is important: in the chunk later we need to match clusters from the edgeR output with the cluster abundances of FlowSOM, we decided
### to do it by indexing, so the order needs to be the same; this circumnavigates numbers v letters etc but is otherwise not very elegant, might change


deg_base_dod6 <- deg_base_dod6[order(as.numeric(deg_base_dod6$Cluster)),]
deg_base_dod <- deg_base_dod[order(as.numeric(deg_base_dod$Cluster)),]
deg_dod_dod6 <- deg_dod_dod6[order(as.numeric(deg_dod_dod6$Cluster)),]
deg_base_c8 <- deg_base_c8[order(as.numeric(deg_base_c8$Cluster)),]
deg_base_c10 <- deg_base_c10[order(as.numeric(deg_base_c10$Cluster)),]

deg_base_dod6$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), "matters_not", "matters"), "matters")
deg_base_dod$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), "matters_not", "matters"), "matters")
deg_dod_dod6$matters <- ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c8$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), ifelse(short[,2] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c10$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), ifelse(short[,3] < 0.01*nrow(flo_set@frames[[files_list[9]]]@exprs), "matters_not", "matters"), "matters")

# add volunteer, cluster column & collate deglists 
deg_base_dod6$Volunteer <- "V07"
deg_base_dod$Volunteer <- "V07"
deg_dod_dod6$Volunteer <- "V07"
deg_base_c8$Volunteer <- "V07"
deg_base_c10$Volunteer <- "V07"

deg_base_dod6$Timepoint <- "deg_base_dod6"
deg_base_dod$Timepoint <- "deg_base_dod"
deg_dod_dod6$Timepoint <- "deg_dod_dod6"
deg_base_c8$Timepoint <- "deg_base_c8"
deg_base_c10$Timepoint <- "deg_base_c10"


all_degs <- list(deg_base_dod6, deg_base_dod, deg_dod_dod6, deg_base_c8, deg_base_c10)

# make a list of lists that contain re-transformed fold change values
fold_change  <- lapply(all_degs, function(x){x$Fold_Change <- 2^x$logFC})

# add list to list of dataframes as new column; the simplify=F argument prevents conversion of every column
# to a list

all_degs <- mapply(cbind, all_degs, "Fold_Change"= fold_change, SIMPLIFY=F)

# combine them in one big file
#individual_from_all <- rbind(deg_base_dod6, deg_base_dod, deg_dod_dod6, deg_base_c8, deg_base_c10)

individual_from_all <- plyr::ldply(all_degs, rbind)
nrow(individual_from_all)
#nrow=320

individual_from_all <- dplyr::filter(individual_from_all, FDR<0.01)
nrow(individual_from_all)
#nrow=205

# subset dataframe so only fold changes over 2 and less than 0.5 are included
lower_cut_off <- dplyr::filter(individual_from_all, Fold_Change > 2)
upper_cut_off <- dplyr::filter(individual_from_all, Fold_Change < 0.5)
cut_off <- rbind(upper_cut_off, lower_cut_off)
nrow(cut_off)
#nrow=15

cut_off <- dplyr::filter(cut_off, matters == "matters")
nrow(cut_off)
#nrow=10

# disassemble big dataframe for making figures

list_of_degs <- split(cut_off, cut_off$Timepoint)
important_ones <- plyr::ldply(list_of_degs, rbind)
#nrow=42 ###  yaaay

######        the plan is to make figures showing a starplot of all their deg clusters with the fold change



setwd("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_medians/")
# setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_medians/")



deg_medians_aggregate  <- read.csv("aggregate_cluster_medians.csv")


#0 1 transform of marker intensities, columnwise
deg_medians_aggregate[,5:71] <- lapply(deg_medians_aggregate[5:71], function(x){scales::rescale(x,to=c(0,1))})


### make small dataframes for each cluster comparison
### nothing significant from base_c8 or base_c10 or base_dod (!)

clusters_base_dod6 <- dplyr::filter(important_ones, Timepoint=="deg_base_dod6")
clusters_dod_dod6 <- dplyr::filter(important_ones, Timepoint=="deg_dod_dod6")

### add marker expression data

medians_base_dod6 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_base_dod6$Cluster)  %>%
  mutate(Comparison = "base_dod6") %>%
  mutate(Fold_Change = clusters_base_dod6$Fold_Change)

medians_dod_dod6 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% clusters_dod_dod6$Cluster)  %>%
  mutate(Comparison = "dod_dod6") %>%
  mutate(Fold_Change = clusters_dod_dod6$Fold_Change)


#put it all together to make ggplots; drop UMAP channels
deg_medians_aggregate <- rbind(medians_base_dod6, medians_dod_dod6)

deg_medians_all <- select(deg_medians_aggregate, colnames(deg_medians_aggregate)[c(1, 2, 5, 16:17, 25:59, 65, 67, 72, 73)])

colnames(deg_medians_all)[4:42] <- substr(colnames(deg_medians_all)[4:42], 8, nchar(colnames(deg_medians_all)[4:42])-10)
colnames(deg_medians_all)[3] <- "CD45"



# convert to long format
# deg_medians_all$MetaclusterID <- NULL
# long_deg_medians_all <- gather(deg_medians_all, Marker, Intensity, colnames(deg_medians_all)[2:41])
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

data <- read.csv("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_abundances.csv")
# data <- read.csv("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/FlowSOM_big_timecourse_07a_results/results/cluster_abundances.csv")

short <- select(data, colnames(data[3:7]))


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

abun_clusters <- rbind(base_dod6_clusters, dod_dod6_clusters)

# order the thing so the fold change carries over correctly
abun_clusters <- abun_clusters[order(as.numeric(abun_clusters$ClusterID)),]
abun_clusters$ClusterID <- as.numeric(abun_clusters$ClusterID)

long_abun_clusters <- gather(abun_clusters, Timepoint, Count, c("pre", "post"))
# long_abun_clusters$Count <- long_abun_clusters$Count/nrow(flo_set@frames[[files_list[24]]]@exprs)

# add correct fold change variable
deg_medians_all$Fold_Change <- abun_clusters$Fold_Change

# long format
long_deg_medians_all <- gather(deg_medians_all, Marker, Intensity, colnames(deg_medians_all)[2:31])

# add categorical variable whether something is going up or down based on fold change
long_deg_medians_all$Direction <- ifelse(long_deg_medians_all$Fold_Change>1, "up", "down")

# reorder that variable to first the up panel is displayed in the plots
long_deg_medians_all$Directions <- factor(long_deg_medians_all$Direction, levels = c("up", "down"))

# 
# # makes a barplot of abundance at the pre and post timepoint for each comparison
# for(i in unique(long_abun_clusters$Comparison)){
#   print(i)
#   
#   sub_set <- dplyr::filter(long_abun_clusters, Comparison == i)
#   specific_levels <- unique(sub_set[order(sub_set$Fold_Change, decreasing = TRUE),"ClusterID"])
#   print(specific_levels)
#   
#   assign(paste(i,"_bar", sep=''), 
#          
#          ggplot(data = sub_set,
#                 aes_(x=factor(sub_set$ClusterID, levels = specific_levels), y=sub_set$Count, fill=factor(sub_set$Timepoint, levels=c("pre", "post")))
#          )+
#            geom_bar(stat="identity", position=position_dodge())+
#            scale_fill_brewer(palette="Paired")+
#            xlab("Cluster ID")+
#            scale_y_continuous()+
#            theme(legend.title = element_blank(),
#                  legend.text = element_text(size = 20),
#                  legend.position = "top", 
#                  legend.justification = "center",
#                  legend.direction = "horizontal",
#                  axis.line = element_line(colour = "black"),
#                  axis.text.x = element_text(size=20, color="black"),
#                  axis.title.x = element_text(size=24, color="black"),
#                  axis.title.y = element_text(size=24, color="black"),
#                  axis.text.y = element_text(size=20, color="black")))
# }
# 
# 
# 
# ##############          working figure
# 
# 
# for(i in unique(long_deg_medians_all$Comparison)){
#   specific_levels <- NULL
#   print(i)
#   # ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
#   # ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
#   # 
#   sub_set <- dplyr::filter(long_deg_medians_all, Comparison == i)
#   sub_set <- sub_set[order(sub_set$Fold_Change, decreasing=TRUE),]
#   specific_levels <- unique(sub_set$ClusterID)
#   
#   print(specific_levels)
#   # specific_levels <- sub_set %>% 
#   #   dplyr::filter(Marker == "CD4") %>%
#   #   arrange(desc(Intensity))
#   # 
#   # specific_levels <- c(as.character(specific_levels$ClusterID))
#   # 
#   assign(paste("comparison_", unique(sub_set$Comparison), sep=''),
#          ggplot(data = sub_set, aes_(x=factor(sub_set$ClusterID, levels = as.character(specific_levels)), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Comparison))+
#            geom_tile(aes(fill=Intensity), color="white")+
#            scale_fill_gradientn(colors=rev(my_palette))+
#            scale_y_discrete(position = "left")+
#            xlab(NULL)+
#            facet_grid(~ Direction, scales = "free")+
#            ggtitle(paste(i))+
#            theme(panel.border = element_blank(),
#                  axis.text.y.left = element_text(size=35),
#                  axis.line.y.left = element_blank(),
#                  axis.line.y.right = element_blank(),
#                  axis.ticks.y = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.text.x = element_text(size = 33),
#                  axis.text.y.right = element_text(size = 35),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  axis.line = element_line(colour = "black"),
#                  legend.title = element_blank(),
#                  legend.position = "none",
#                  plot.title = element_text(size = 45, hjust = 0.5),
#                  plot.margin = unit(c(1,0,1,0), "cm"),
#                  strip.text.x = element_text(size=28))
#   )
# } 
# 
# 
# 
# ggsave("sandbox.pdf", grid.arrange(comparison_base_dod, comparison_base_dod6, comparison_dod_dod6, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
#                                                                                                                                                              c(4,5,6))
# ),  width = 40, height = 40, limitsize = F)
# 
# 
# #ggsave("heatmap_plus_abundance_base_c8.pdf", grid.arrange(comparison_base_c8, base_c8_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
# ggsave("heatmap_plus_abundance_base_dod.pdf", grid.arrange(comparison_base_dod, base_dod_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
# ggsave("heatmap_plus_abundance_base_dod6.pdf", grid.arrange(comparison_base_dod6, base_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
# ggsave("heatmap_plus_abundance_dod_dod6.pdf", grid.arrange(comparison_dod_dod6, dod_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
# 
# 












########   combine the two for loops so the order f the bar graph matches the drapes

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
           ggtitle(paste("V07 ", i, sep=''))+
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

setwd("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/figures/")

# ggsave("v07_heatmap_plus_abundance_base_dod6.pdf", grid.arrange(comparison_base_dod6, base_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
# ggsave("v07_heatmap_plus_abundance_dod_dod6.pdf", grid.arrange(comparison_dod_dod6, dod_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)

ggsave("v07_01_heatmap_plus_abundance_base_dod6.pdf", grid.arrange(comparison_base_dod6, base_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("v07_01_heatmap_plus_abundance_dod_dod6.pdf", grid.arrange(comparison_dod_dod6, dod_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)



# figure only with only what's up from dod to dod6
mat <- spread(long_deg_medians_all, Marker, Intensity)

mat2 <- mat %>%
  dplyr::filter(., Comparison=="dod_dod6") %>%
  dplyr::filter(., Direction=="up")%>%
  dplyr::select(., -Comparison, -Fold_Change, -Direction, -Directions)

mat2 <- as.matrix(mat2)

tmat <- t(mat2)

colnames(tmat) <- tmat[1,]

corr_mat=cor(tmat[2:nrow(tmat),],method="s")


##### heatmap

sub_set <- dplyr::filter(long_deg_medians_all, Comparison == "dod_dod6")
sub_set <- dplyr::filter(sub_set, Direction=="up")

sub_set$ClusterID <- as.character(sub_set$ClusterID)

specific_levels <- rownames(corr_mat[order(corr_mat[,1], decreasing = TRUE),])

up_v7 <- ggplot(data = sub_set, aes_(x=factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Comparison))+
  geom_tile(aes(fill=Intensity), color="white")+
  scale_fill_gradientn(colors=rev(my_palette))+
  scale_y_discrete(position = "left")+
  xlab("Cluster ID")+
  ggtitle("Volunteer 07")+
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

#setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
setwd("/Users/s1249052/PhD/cytof/better_gating/double_flowsoms/figures/")

ggsave("v07_up.pdf", height=9, width=11)
ggsave("v07_up.png", height=9, width=11)






up_dod_dod6 <- filter(long_abun_clusters, Comparison=="dod_dod6")
up_dod_dod6 <- filter(up_dod_dod6, Fold_Change>1)
up_dod_dod6 <- filter(up_dod_dod6, Timepoint=="post")

up_dod_dod6[,c(2,4)] <- NULL
up_dod_dod6$Count <- up_dod_dod6$Count*100
up_dod_dod6$SubSet <- c("CD4", "CD4", "CD4")
up_dod_dod6$SubSet <- factor(up_dod_dod6$SubSet, levels= c("CD4", "CD8", "Vd2", "MAIT", "DN"))

up_dod_dod6 <- up_dod_dod6[order(up_dod_dod6$SubSet),]
#up_dod_dod6$ClusterID <- paste0("Cluster_", up_dod_dod6$ClusterID)

up_dod_dod6$ymin[1]<-0
up_dod_dod6$ymin[2:nrow(up_dod_dod6)] <- cumsum(up_dod_dod6$Count)[1:nrow(up_dod_dod6)-1]
up_dod_dod6$ymax <- up_dod_dod6$ymin + up_dod_dod6$Count

#my_long_palette <- c(rev(sequential_hcl(5, "Heat")),qualitative_hcl(100, "Cold"), qualitative_hcl(150, "Dynamic"))


my_palette <- c(qualitative_hcl(10, "Dynamic"), rev(qualitative_hcl(5, "Dark3")))
ggplot()+
  # circles made of numbers indicating fold change
  geom_text(aes(x=4.05, y=seq(0,max(up_dod_dod6$ymax), by=max(up_dod_dod6$ymax)/nrow(up_dod_dod6)/8), label=paste0("2")), size=1.2)+
  geom_text(aes(x=5.05, y=seq(0,max(up_dod_dod6$ymax), by=max(up_dod_dod6$ymax)/nrow(up_dod_dod6)/8), label=paste0("4")), size=1.2)+
  geom_text(aes(x=6.05, y=seq(0,max(up_dod_dod6$ymax), by=max(up_dod_dod6$ymax)/nrow(up_dod_dod6)/8), label=paste0("8")), size=1.2)+
  geom_text(aes(x=7.05, y=seq(0,max(up_dod_dod6$ymax), by=max(up_dod_dod6$ymax)/nrow(up_dod_dod6)/8), label=paste0("16")), size=1.2)+
  geom_text(aes(x=8.05, y=seq(0,max(up_dod_dod6$ymax), by=max(up_dod_dod6$ymax)/nrow(up_dod_dod6)/8), label=paste0("32")), size=1.2)+
  # clusters and subset pies
  geom_rect(data=up_dod_dod6, aes(fill=as.numeric(factor(ClusterID)), ymin=ymin, ymax=ymax, xmax=6+log2(up_dod_dod6$Fold_Change), xmin=3+log2(up_dod_dod6$Fold_Change)))+
  geom_rect(data=up_dod_dod6, aes(fill=as.numeric(factor(SubSet, levels=levels(SubSet)))+10, ymin=ymin, ymax=ymax, xmax=13, xmin=11))+
  #labelling those pies
  #geom_text(aes(x=4.5+log2(up_dod_dod6$Fold_Change), y=(ymin+ymax)/2, label = paste("C_", ClusterID, sep="")),size=2.4,fontface="bold")+
  geom_text(data=up_dod_dod6, aes(x=12, y=(ymin+ymax)/2, label = paste(SubSet)), size=2, fontface="bold")+
  #beautify
  #scale_fill_gradientn(colors=c(qualitative_hcl(10, "Harmonic"), qualitative_hcl(15, "Dark3")), values=c(scales::rescale(seq(1,10), to=c(0,0.66)), scales::rescale(seq(10,15), to=c(0.66,1))))+
  scale_fill_gradientn(guide="none", colors=my_palette)+
  theme(aspect.ratio=1,
        legend.title = element_text(),
        legend.title.align = 0.5,
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  xlim(c(0, 14))+
  ggtitle("Volunteer 07")+
  coord_polar("y")

setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
up_dod_dod6$Volunteer <- "V07"
write.csv(up_dod_dod6, "v07_up_dod_dod6.csv")



setwd("C:/Users/Florian/PhD/cytof/vac69a/double_flowsoms/figures")
ggsave("v07_hot_pie.pdf")
