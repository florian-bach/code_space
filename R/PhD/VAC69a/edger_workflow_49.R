library(edgeR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)


# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# read in data 
data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_vac69a_t_cells_49cluster_prepost_1_results/results/cluster_abundances.csv")
data2 <-read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_vac69a_t_cells_49cluster_prepost_2_results/results/cluster_abundances.csv")

# get rid of uneccessary columns, convert fractions to counts
data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_vac69a_t_cells_49cluster_prepost_1_results/results/cluster_abundances.csv")
short <- select(data, colnames(data[3:14]))

number_of_cells <- 5534

short <- short*number_of_cells # this is the number of cells from each fcs file

# repeat same stuff for other dataframe 
short2 <- select(data2, colnames(data[3:14]))
short2 <- short2*number_of_cells # this is the number of cells from each fcs file

# combine
short <- cbind(short, short2)

# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
groups <- rep(substr(colnames(short)[1:12],1, nchar(colnames(short))-9), times=2)

#groups <- gsub("T.6", "DoD.6", groups, fixed=T) 

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)
colnames(design) <- substr(colnames(short)[1:12],1, nchar(colnames(short))-9)

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
pre_post_02 <- makeContrasts(T.6_02 - baseline_02,levels=design)
pre_post_03 <- makeContrasts(T.6_03 - baseline_03,levels=design)
pre_post_05 <- makeContrasts(T.6_05 - baseline_05,levels=design)
pre_post_06 <- makeContrasts(T.6_06 - baseline_06,levels=design)
pre_post_07 <- makeContrasts(T.6_07 - baseline_07,levels=design)
pre_post_09 <- makeContrasts(T.6_09 - baseline_09,levels=design)

# make deg list for each person                                                                             
deg_02 <- glmQLFTest(fit, contrast=pre_post_02) 
deg_03 <- glmQLFTest(fit, contrast=pre_post_03) 
deg_05 <- glmQLFTest(fit, contrast=pre_post_05) 
deg_06 <- glmQLFTest(fit, contrast=pre_post_06) 
deg_07 <- glmQLFTest(fit, contrast=pre_post_07) 
deg_09 <- glmQLFTest(fit, contrast=pre_post_09) 

# save as dataframes
deg_02 <- topTags(deg_02, n=49)$table
deg_03 <- topTags(deg_03, n=49)$table
deg_05 <- topTags(deg_05, n=49)$table 
deg_06 <- topTags(deg_06, n=49)$table 
deg_07 <- topTags(deg_07, n=49)$table 
deg_09 <- topTags(deg_09, n=49)$table

#add cluser as column 
deg_02$Cluster <- rownames(deg_02)
deg_03$Cluster <- rownames(deg_03)
deg_05$Cluster <- rownames(deg_05)
deg_06$Cluster <- rownames(deg_06)
deg_07$Cluster <- rownames(deg_07)
deg_09$Cluster <- rownames(deg_09)

deg_02 <- deg_02[order(as.numeric(deg_02$Cluster)),]; deg_02$Baseline <- short[,7]; deg_02$Treatment <- short[,1]
deg_03 <- deg_03[order(as.numeric(deg_03$Cluster)),]; deg_03$Baseline <- short[,8]; deg_03$Treatment <- short[,2]
deg_05 <- deg_05[order(as.numeric(deg_05$Cluster)),]; deg_05$Baseline <- short[,9]; deg_05$Treatment <- short[,3]
deg_06 <- deg_06[order(as.numeric(deg_06$Cluster)),]; deg_06$Baseline <- short[,10]; deg_06$Treatment <- short[,4]
deg_07 <- deg_07[order(as.numeric(deg_07$Cluster)),]; deg_07$Baseline <- short[,11]; deg_07$Treatment <- short[,5]
deg_09 <- deg_09[order(as.numeric(deg_09$Cluster)),]; deg_09$Baseline <- short[,12]; deg_09$Treatment <- short[,6]

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

# disassemble big dataframe for making figures

list_of_degs <- split(cut_off, cut_off$Volunteer)



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change

setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_vac69a_t_cells_49cluster_prepost_1_results/results/cluster_medians")

median_02 <- read.csv("T+6_02_UMAP_cluster_medians.csv")
median_03 <- read.csv("T+6_03_UMAP_cluster_medians.csv")
median_05 <- read.csv("T+6_05_UMAP_cluster_medians.csv")
median_06 <- read.csv("T+6_06_UMAP_cluster_medians.csv")
median_07 <- read.csv("T+6_07_UMAP_cluster_medians.csv")
median_09 <- read.csv("T+6_09_UMAP_cluster_medians.csv")

# super convoluted way of doing it but it works: restrict cluster medians to clusters that can be found in deg list with cutoff of log2 1 and -1
deg_medians_02 <- median_02 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[1]]$Cluster)) %>%
  mutate(Volunteer = "02") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[1]]$logFC)),0,4), ")", sep=''))

deg_medians_03 <- median_03 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[2]]$Cluster)) %>%
  mutate(Volunteer = "03") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[2]]$logFC)),0,4), ")", sep=''))

deg_medians_05 <- median_05 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[3]]$Cluster)) %>%
  mutate(Volunteer = "05") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[3]]$logFC)),0,4), ")", sep=''))

deg_medians_06 <- median_06 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[4]]$Cluster)) %>%
  mutate(Volunteer = "06") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[4]]$logFC)),0,4), ")", sep=''))


deg_medians_07 <- median_07 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[5]]$Cluster)) %>%
  mutate(Volunteer = "07") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[5]]$logFC)),0,4), ")", sep=''))

deg_medians_09 <- median_09 %>%
  filter(ClusterID %in% as.numeric(list_of_degs[[6]]$Cluster)) %>%
  mutate(Volunteer = "09") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[6]]$logFC)),0,4), ")", sep=''))

#put it all together to make ggplots; drop UMAP channels
deg_medians_all <- rbind(deg_medians_02, deg_medians_03, deg_medians_05, deg_medians_06, deg_medians_07, deg_medians_09)
data_medians_all <- select(deg_medians_all, colnames(deg_medians_all)[c(1:36, 39:40)])

#make column names marker names
colnames(data_medians_all) <- c('ClusterID', 'MetaclusterID', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
                          '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
                          '156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
                          '165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
                          '174Yb_CLA', '175Lu_Perforin', '198Pt_CD8', '209Bi_CD16', 'Volunteer', 'Cluster (Fold Change)')

colnames(data_medians_all)[3:36] <- substr(colnames(data_medians_all)[3:36], 7, 30)

# convert to long format
data_medians_all$MetaclusterID <- NULL
long_deg_medians_all <- gather(data_medians_all, Marker, Intensity, colnames(data_medians_all)[2:35])

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

#0-1 transform is pointless


my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")
#color_blind <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")

# arcsinh transform
# (big_heatmap_arsinh <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
#     geom_tile(aes(fill=Intensity), color="white")+
#     scale_fill_gradientn(colors=rev(my_palette))+
#     xlab("Cluster")+
#     ylab("Marker")+
#     ggtitle("\"Differentially Expressed\" Cluster Medians Per Volunteer in VAC69a")+
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


##############          working figure


for(i in unique(long_deg_medians_all$Volunteer)){
  specific_levels <- NULL
  ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
  ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
  
  sub_set <- filter(long_deg_medians_all, Volunteer == i)
  
  specific_levels <- sub_set %>% 
    filter(Marker == "CD4") %>%
    arrange(desc(Intensity))
  
  specific_levels <- c(as.character(specific_levels$ClusterID))
  
  assign(paste("Volunteer_", unique(sub_set$Volunteer), sep=''),
         ggplot(data = sub_set, aes_(x = factor(sub_set$ClusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$Volunteer))+
           geom_tile(aes(fill=Intensity), color="white")+
           scale_fill_gradientn(colors=rev(my_palette))+
           scale_y_discrete(position = result1)+
           xlab(NULL)+
           ggtitle(paste("Volunteer ", i,sep=''))+
           theme(panel.border = element_blank(),
                 axis.text.y.left = result,
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



ggsave("sandbox.pdf", grid.arrange(Volunteer_02, Volunteer_03, Volunteer_05, Volunteer_06, Volunteer_07, Volunteer_09, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
                                                                                                                                                             c(4,5,6))
),  width = 40, height = 40, limitsize = F)

#######         figures for cluster abundances

data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_vac69a_t_cells_49cluster_prepost_1_results/results/cluster_abundances.csv")
short <- select(data, colnames(data[3:14]))


clusters_02 <- as.integer(list_of_degs[[1]]$Cluster)
clusters_03 <- as.integer(list_of_degs[[2]]$Cluster)
clusters_05 <- as.integer(list_of_degs[[3]]$Cluster)
clusters_06 <- as.integer(list_of_degs[[4]]$Cluster)
clusters_07 <- as.integer(list_of_degs[[5]]$Cluster)
clusters_09 <- as.integer(list_of_degs[[6]]$Cluster)

abun_clusters_02 <- short[c(clusters_02),c(7,1)]
abun_clusters_02$ClusterID <- as.factor(clusters_02)
abun_clusters_02$Volunteer <- substr(colnames(abun_clusters_02)[2], 5, 6) 
colnames(abun_clusters_02) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_03 <- short[c(clusters_03),c(8,2)]
abun_clusters_03$ClusterID <- as.factor(clusters_03)
abun_clusters_03$Volunteer <- substr(colnames(abun_clusters_03)[2], 5, 6) 
colnames(abun_clusters_03) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_05 <- short[c(clusters_05),c(9,3)]
abun_clusters_05$ClusterID <- as.factor(clusters_05)
abun_clusters_05$Volunteer <- substr(colnames(abun_clusters_05)[2], 5, 6) 
colnames(abun_clusters_05) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_06 <- short[c(clusters_06),c(10,4)]
abun_clusters_06$ClusterID <- as.factor(clusters_06)
abun_clusters_06$Volunteer <- substr(colnames(abun_clusters_06)[2], 5, 6) 
colnames(abun_clusters_06) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_07 <- short[c(clusters_07),c(11,5)]
abun_clusters_07$ClusterID <- as.factor(clusters_07)
abun_clusters_07$Volunteer <- substr(colnames(abun_clusters_07)[2], 5, 6) 
colnames(abun_clusters_07) <- c("Baseline", "T+6", "ClusterID", "Volunteer")

abun_clusters_09 <- short[c(clusters_09),c(12,6)]
abun_clusters_09$ClusterID <- as.factor(clusters_09)
abun_clusters_09$Volunteer <- substr(colnames(abun_clusters_09)[2], 5, 6) 
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

