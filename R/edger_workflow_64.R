library(edgeR)
library(dplyr)
library(tidyr)
library(flowCore)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)


# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# map new data to same some, but have different seed sequence should produce similar, but not identical results


# read in data 
data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_big_timecourse_03_(copy)_(copy)_results/results/cluster_abundances.csv")
data2 <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_big_timecourse_03_(copy)_(copy)_(copy)_results/results/cluster_abundances.csv")

#extract number of cells in each fcs file to convert frequency to actual number
setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/better_gating/")
files_list <- list.files(path=".", pattern="*.fcs")

flo_set <- read.flowSet(files_list[21:25], transformation = FALSE, truncate_max_range = FALSE)


short <- data

# convert to absolute numbers: we're keeping them the same because we still need to check whether edgeR can
# adequately deal with the lib.size argument differing with the fcs file sizel this is just the smallest number
# from any of the fcs files

short[,3] <- data[,3]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short[,4] <- data[,4]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short[,5] <- data[,5]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short[,6] <- data[,6]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short[,7] <- data[,7]*nrow(flo_set@frames[[files_list[24]]]@exprs)

# short[,3] <- data[,3]*nrow(flo_set@frames[[files_list[21]]]@exprs)
# short[,4] <- data[,4]*nrow(flo_set@frames[[files_list[22]]]@exprs)
# short[,5] <- data[,5]*nrow(flo_set@frames[[files_list[23]]]@exprs)
# short[,6] <- data[,6]*nrow(flo_set@frames[[files_list[24]]]@exprs)
# short[,7] <- data[,7]*nrow(flo_set@frames[[files_list[25]]]@exprs)


# repeat same stuff for other dataframe 
short2 <- data2
short2[,3] <- data2[,3]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short2[,4] <- data2[,4]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short2[,5] <- data2[,5]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short2[,6] <- data2[,6]*nrow(flo_set@frames[[files_list[24]]]@exprs)
short2[,7] <- data2[,7]*nrow(flo_set@frames[[files_list[24]]]@exprs)
# combine
short <- cbind(short, short2)
rownames(short) <- short$ClusterID
short[,c(1:2,8:9)] <-NULL
# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
#groups <- rep(substr(colnames(short)[1:14],1, nchar(colnames(short))-9), times=2)
groups <- substr(colnames(short), 2,6)

colnames(short) <- substr(colnames(short), 2,6)
colnames(short)[6:10] <- paste(colnames(short)[1:5], "_1", sep='')


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
base_dod6 <- makeContrasts(groups14c10 - groups14c06,levels=design)
base_dod <- makeContrasts(groups14c09 - groups14c06,levels=design)
dod_dod6 <- makeContrasts(groups14c10 - groups14c09,levels=design)
base_c8 <- makeContrasts(groups14c07 - groups14c06,levels=design)
base_c10 <- makeContrasts(groups14c08 - groups14c06,levels=design)


# make deg list for each person                                                                             
deg_base_dod6 <- glmQLFTest(fit, contrast=base_dod6) 
deg_base_dod <- glmQLFTest(fit, contrast=base_dod) 
deg_dod_dod6 <- glmQLFTest(fit, contrast=dod_dod6) 
deg_base_c8 <- glmQLFTest(fit, contrast=dod_dod6) 
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

deg_02 <- deg_02[order(as.numeric(deg_02$Cluster)),]; deg_02$Baseline <- short[,7]; deg_02$Treatment <- short[,1]
deg_03 <- deg_03[order(as.numeric(deg_03$Cluster)),]; deg_03$Baseline <- short[,8]; deg_03$Treatment <- short[,2]
deg_05 <- deg_05[order(as.numeric(deg_05$Cluster)),]; deg_05$Baseline <- short[,9]; deg_05$Treatment <- short[,3]
deg_06 <- deg_06[order(as.numeric(deg_06$Cluster)),]; deg_06$Baseline <- short[,10]; deg_06$Treatment <- short[,4]
deg_07 <- deg_07[order(as.numeric(deg_07$Cluster)),]; deg_07$Baseline <- short[,11]; deg_07$Treatment <- short[,5]
deg_09 <- deg_09[order(as.numeric(deg_09$Cluster)),]; deg_09$Baseline <- short[,12]; deg_09$Treatment <- short[,6]

#mark stuff that never surpasses 1% frequency
deg_02$matters <- ifelse(deg_02$Baseline < 0.01*number_of_cells, ifelse(deg_02$Treatment < 0.01*number_of_cells, "matters_not", "matters"), "matters")


deg_base_dod6 <- deg_base_dod6[order(as.numeric(deg_base_dod6$Cluster)),]
deg_base_dod <- deg_base_dod[order(as.numeric(deg_base_dod$Cluster)),]
deg_dod_dod6 <- deg_dod_dod6[order(as.numeric(deg_dod_dod6$Cluster)),]
deg_base_c8 <- deg_base_c8[order(as.numeric(deg_base_c8$Cluster)),]
deg_base_c10 <- deg_base_c10[order(as.numeric(deg_base_c10$Cluster)),]

deg_base_dod6$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), "matters_not", "matters"), "matters")
deg_base_dod$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), "matters_not", "matters"), "matters")
deg_dod_dod6$matters <- ifelse(short[,4] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), ifelse(short[,5] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c8$matters <- ifelse(short[,2] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), ifelse(short[,2] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), "matters_not", "matters"), "matters")
deg_base_c10$matters <- ifelse(short[,1] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), ifelse(short[,3] < 0.01*nrow(flo_set@frames[[files_list[24]]]@exprs), "matters_not", "matters"), "matters")

# add volunteer, cluster column & collate deglists 
deg_base_dod6$Volunteer <- "V03"
deg_base_dod$Volunteer <- "V03"
deg_dod_dod6$Volunteer <- "V03"
deg_base_c8$Volunteer <- "V03"
deg_base_c10$Volunteer <- "V03"

deg_base_dod6$Timepoint <- "deg_base_dod6"
deg_base_dod$Timepoint <- "deg_base_dod"
deg_dod_dod6$Timepoint <- "deg_dod_dod6"
deg_base_c8$Timepoint <- "deg_base_c8"
deg_base_c10$Timepoint <- "deg_base_c10"


# deg_02$Cluster <- rownames(deg_02)
# deg_03$Cluster <- rownames(deg_03)
# deg_05$Cluster <- rownames(deg_05)
# deg_06$Cluster <- rownames(deg_06)
# deg_07$Cluster <- rownames(deg_07)
# deg_09$Cluster <- rownames(deg_09)

# combine them in one big file
individual_from_all <- rbind(deg_base_dod6, deg_base_dod, deg_dod_dod6, deg_base_c8, deg_base_c10)

# subset dataframe so only fold changes over 2 and less than 0.5 are included
upper_cut_off <- dplyr::filter(individual_from_all, logFC > 1)
lower_cut_off <- dplyr::filter(individual_from_all, individual_from_all$logFC < -1)
cut_off <- rbind(upper_cut_off, lower_cut_off)
cut_off <- dplyr::filter(cut_off, matters == "matters")
nrow(cut_off)

# disassemble big dataframe for making figures

list_of_degs <- split(cut_off, cut_off$Timepoint)



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change



setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_big_timecourse_03_(copy)_(copy)_results/results/cluster_medians/")

deg_medians_aggregate  <- read.csv("aggregate_cluster_medians.csv")

# super convoluted way of doing it but it works: restrict cluster medians to clusters that can be found in deg list with cutoff of log2 1 and -1
deg_medians_02 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[1]]$Cluster)) %>%
  mutate(Comparison = list_of_degs[[1]]$Timepoint) %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[1]]$logFC)),0,4), ")", sep=''))

deg_medians_03 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[2]]$Cluster))  %>%
  mutate(Comparison = list_of_degs[[2]]$Timepoint) %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[2]]$logFC)),0,4), ")", sep=''))

deg_medians_05 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[3]]$Cluster))  %>%
  mutate(Comparison = list_of_degs[[3]]$Timepoint) %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[3]]$logFC)),0,4), ")", sep=''))

deg_medians_06 <- deg_medians_aggregate %>%
  dplyr::filter(ClusterID %in% as.numeric(list_of_degs[[4]]$Cluster))  %>%
  mutate(Comparison = list_of_degs[[4]]$Timepoint) %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[4]]$logFC)),0,4), ")", sep=''))


#put it all together to make ggplots; drop UMAP channels
deg_medians_aggregate <- rbind(deg_medians_02, deg_medians_03, deg_medians_05, deg_medians_06)

deg_medians_all <- select(deg_medians_aggregate, colnames(deg_medians_aggregate)[c(1, 2, 5, 16:17, 25:59, 65, 67, 72)])


colnames(deg_medians_all)[4:42] <- substr(colnames(deg_medians_all)[4:42], 8, nchar(colnames(deg_medians_all)[4:42])-10)
colnames(deg_medians_all)[3] <- "CD45"
# convert to long format
deg_medians_all$MetaclusterID <- NULL
long_deg_medians_all <- gather(deg_medians_all, Marker, Intensity, colnames(deg_medians_all)[2:41])

# make beautiful iris color palette

#my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)
#my_palette <- rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))


# make levels to reorder markers in a meaningful way

marker_levels <- c("CD4",
"CD8",
"Vd2",
"Va7.2",
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
"Ki.67",
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


my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")
#color_blind <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")



##############          working figure


for(i in unique(long_deg_medians_all$Comparison)){
  specific_levels <- NULL
  # ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  # ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
  # 
  sub_set <- dplyr::filter(long_deg_medians_all, Comparison == i)
  
  specific_levels <- sub_set %>% 
    dplyr::filter(Marker == "CD4") %>%
    arrange(desc(Intensity))
  
  specific_levels <- c(as.character(specific_levels$ClusterID))
  
  assign(paste("comparison_", unique(sub_set$Comparison), sep=''),
         ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Comparison))+
           geom_tile(aes(fill=Intensity), color="white")+
           scale_fill_gradientn(colors=rev(my_palette))+
           scale_y_discrete(position = "left")+
           xlab(NULL)+
           ggtitle(paste("Volunteer ", i,sep=''))+
           theme(panel.border = element_blank(),
                 axis.text.y.left = element_text(size=35),
                 axis.line.y.left = element_blank(),
                 axis.line.y.right = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size = 33),
                 axis.text.y.right = element_text(size = 35),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.title = element_blank(),
                 legend.position = "none",
                 plot.title = element_text(size = 45, hjust = 0.5),
                 plot.margin = unit(c(1,0,1,0), "cm"))
  )
} 



ggsave("sandbox.pdf", grid.arrange(comparison_deg_base_c8, comparison_deg_base_dod, comparison_deg_base_dod6, comparison_deg_dod_dod6, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                             c(4,5,6))
),  width = 40, height = 40, limitsize = F)

#######         figures for cluster abundances

data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_big_timecourse_03_(copy)_(copy)_results/results/cluster_abundances.csv")
short <- select(data, colnames(data[3:14]))


clusters_base_c8 <- as.integer(list_of_degs[[1]]$Cluster)
clusters_base_dod <- as.integer(list_of_degs[[2]]$Cluster)
clusters_base_dod6 <- as.integer(list_of_degs[[3]]$Cluster)
clusters_dod_dod6 <- as.integer(list_of_degs[[4]]$Cluster)


abun_clusters_02 <- short[c(clusters_base_c8),c(2,1)]
abun_clusters_02$ClusterID <- as.factor(clusters_base_c8)
abun_clusters_02$Comparison <- "base_c8"
colnames(abun_clusters_02) <- c("post", "pre", "ClusterID", "Comparison")

abun_clusters_03 <- short[c(clusters_base_dod),c(4,1)]
abun_clusters_03$ClusterID <- as.factor(clusters_base_dod)
abun_clusters_03$Comparison <- "base_dod"
colnames(abun_clusters_03) <- c("post", "pre", "ClusterID", "Comparison")

abun_clusters_05 <- short[c(clusters_base_dod6),c(5,1)]
abun_clusters_05$ClusterID <- as.factor(clusters_base_dod6)
abun_clusters_05$Comparison <- "base_dod6"
colnames(abun_clusters_05) <- c("post", "pre", "ClusterID", "Comparison")

abun_clusters_06 <- short[c(clusters_dod_dod6),c(5,4)]
abun_clusters_06$ClusterID <- as.factor(clusters_dod_dod6)
abun_clusters_06$Comparison <- "dod_dod6"
colnames(abun_clusters_06) <- c("post", "pre", "ClusterID", "Comparison")


abun_clusters <- rbind(abun_clusters_02, abun_clusters_03, abun_clusters_05, abun_clusters_06)

long_abun_clusters <- gather(abun_clusters, Timepoint, Count, c("pre", "post"))
long_abun_clusters$Count <- long_abun_clusters$Count/nrow(flo_set@frames[[files_list[24]]]@exprs)


for(i in unique(long_abun_clusters$Comparison)){
  assign(paste(i,"_bar", sep=''), 
    
  ggplot(data = dplyr::filter(long_abun_clusters, Comparison == i),
         aes(x=factor(ClusterID), y=Count, fill=factor(Timepoint, levels=c("pre", "post")))
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
}

ggsave("heatmap_plus_abundance_base_c8.pdf", grid.arrange(comparison_deg_base_c8, base_c8_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_base_dod.pdf", grid.arrange(comparison_deg_base_dod, base_dod_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_base_dod6.pdf", grid.arrange(comparison_deg_base_dod6, base_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)
ggsave("heatmap_plus_abundance_dod_dod6.pdf", grid.arrange(comparison_deg_dod_dod6, dod_dod6_bar, layout_matrix = rbind(c(1,1,NA),c(1,1,2),c(1,1,NA))), height = 20, width=28)



########   combine the two for loops so the order f the bar graph matches the drapes
