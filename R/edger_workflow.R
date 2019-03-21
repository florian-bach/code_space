library(edgeR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# read in data 
data <- read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_first_flowsom_(copy)_(copy)_results/results/cluster_abundances.csv")
data2 <-read.csv("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/FlowSOM_first_flowsom_(copy)_(copy)_(copy)_results/results/cluster_abundances.csv")

# get rid of uneccessary columns, convert fractions to counts
short <- select(data, colnames(data[3:12]))
short <- short*26027 # this is the number of cells from each fcs file

# add rownames (why does that matter?)
rownames(short)<-paste("Cluster_", data$ClusterID, sep='')

# repeat same stuff for other dataframe 
short2 <- select(data2, colnames(data[3:12]))
short2 <- short2*26027 # this is the number of cells from each fcs file
rownames(short2)<-paste("Cluster_", data2$ClusterID, sep='')

# combine
short <- cbind(short, short2)



# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
groups <- rep(substr(colnames(short)[1:10],1, nchar(colnames(short))-9), times=2)
groups <- gsub("T.6", "DoD.6", groups, fixed=T) 

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)
colnames(design) <- substr(colnames(short)[1:10],1, nchar(colnames(short))-9)

# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short, group=groups, lib.size= colSums(short))
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
pre_post_07 <- makeContrasts(T.6_07 - baseline_07,levels=design)
pre_post_09 <- makeContrasts(T.6_09 - baseline_09,levels=design)

# make deg list for each person                                                                             
dge_02 <- glmQLFTest(fit, contrast=pre_post_02) 
dge_03 <- glmQLFTest(fit, contrast=pre_post_03) 
dge_05 <- glmQLFTest(fit, contrast=pre_post_05) 
dge_07 <- glmQLFTest(fit, contrast=pre_post_07) 
dge_09 <- glmQLFTest(fit, contrast=pre_post_09) 

# save as dataframes
dge_02 <- topTags(dge_02, n=36)$table
dge_03 <- topTags(dge_03, n=36)$table
dge_05 <- topTags(dge_05, n=36)$table 
dge_07 <- topTags(dge_07, n=36)$table 
dge_09 <- topTags(dge_09, n=36)$table

# add volunteer, cluster column & collate dgelists 
dge_02$Volunteer <- "V_02"
dge_03$Volunteer <- "V_03"
dge_05$Volunteer <- "V_05"
dge_07$Volunteer <- "V_07"
dge_09$Volunteer <- "V_09"

dge_02$Cluster <- rownames(dge_02)
dge_03$Cluster <- rownames(dge_03)
dge_05$Cluster <- rownames(dge_05)
dge_07$Cluster <- rownames(dge_07)
dge_09$Cluster <- rownames(dge_09)

# combine them in one big file
individual_from_all <- rbind(dge_02, dge_03, dge_05, dge_07, dge_09)

# subset dataframe so only fold changes over 2 and less than 0.5 are included
upper_cut_off <- filter(individual_from_all, logFC > 1)
lower_cut_off <- filter(individual_from_all, logFC < -1)
cut_off <- rbind(upper_cut_off, lower_cut_off)

# disassemble big dataframe for making figures

list_of_degs <- split(cut_off, cut_off$Volunteer)



######        the plan is to make figures showing a starplot of all their deg clusters with the fold change

setwd("/Users/s1249052/PhD/flow data/vac69a/t cells only/FlowSOM_first_flowsom_(copy)_(copy)_results/results/cluster_medians")

median_02 <- read.csv("T+6_02_UMAP_cluster_medians.csv")
median_03 <- read.csv("T+6_03_UMAP_cluster_medians.csv")
median_05 <- read.csv("T+6_05_UMAP_cluster_medians.csv")
median_07 <- read.csv("T+6_07_UMAP_cluster_medians.csv")
median_09 <- read.csv("T+6_09_UMAP_cluster_medians.csv")

# super convoluted way of doing it but it works: restrict cluster medians to clusters that can be found in deg list with cutoff of log2 1 and -1
deg_medians_02 <- median_02 %>%
  filter(ClusterID %in% as.numeric(substr(list_of_degs[[1]]$Cluster, 9,11))) %>%
  mutate(Volunteer = "02") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[1]]$logFC)),0,4), ")", sep=''))

deg_medians_03 <- median_03 %>%
  filter(ClusterID %in% as.numeric(substr(list_of_degs[[2]]$Cluster, 9,11))) %>%
  mutate(Volunteer = "03") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[2]]$logFC)),0,4), ")", sep=''))

deg_medians_05 <- median_05 %>%
  filter(ClusterID %in% as.numeric(substr(list_of_degs[[3]]$Cluster, 9,11))) %>%
  mutate(Volunteer = "05") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[3]]$logFC)),0,4), ")", sep=''))

deg_medians_07 <- median_07 %>%
  filter(ClusterID %in% as.numeric(substr(list_of_degs[[4]]$Cluster, 9,11))) %>%
  mutate(Volunteer = "07") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[4]]$logFC)),0,4), ")", sep=''))

deg_medians_09 <- median_09 %>%
  filter(ClusterID %in% as.numeric(substr(list_of_degs[[5]]$Cluster, 9,11))) %>%
  mutate(Volunteer = "09") %>%
  mutate(FC = paste("Cluster_", ClusterID, " (", substr(as.character(2^as.numeric(list_of_degs[[5]]$logFC)),0,4), ")", sep=''))

#put it all together to make ggplots; drop UMAP channels
deg_medians_all <- rbind(deg_medians_02, deg_medians_03, deg_medians_05, deg_medians_07, deg_medians_09)
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

my_palette <- colorRampPalette(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))(n=20)
my_palette <- rev(c("#ffd4c4", "#ffb5e2", "#ce70e8", "#a347e5", "#6a39ff"))


# make levels to reorder markers in a meaningful way

levels <- c("CD4",
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

labs <- c(
  "02" = "Volunteer 02",
  "03" = "Volunteer 03",
  "05" = "Volunteer 05",
  "07" = "Volunteer 07",
  "09" = "Volunteer 09"
)

#0-1 transform
(big_heatmap_01 <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
  geom_tile(aes(fill=scales:::rescale(long_deg_medians_all$Intensity, to=c(0,1))), color="white")+
  scale_x_discrete()+
  scale_fill_gradientn(colors=rev(my_palette))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")+
   facet_wrap(~ Volunteer, ncol=5, labeller=as_labeller(labs))
)

my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")
color_blind <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")

# arcsinh transform
(big_heatmap_arsinh <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
    geom_tile(aes(fill=Intensity), color="white")+
    scale_fill_gradientn(colors=rev(my_palette))+
    xlab("Cluster")+
    ylab("Marker")+
    ggtitle("\"Differentially Expressed\" Cluster Medians Per Volunteer in VAC69a")+
    theme(#axis.text.x = element_text(hjust = 1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Volunteer, ncol=5, labeller=as_labeller(labs))
)

# arcsinh transform, accesible color scheme
(big_heatmap_arsinh <- ggplot(data=long_deg_medians_all, aes(x=factor(ClusterID), y=factor(Marker, levels=rev(levels)), group=Volunteer))+
    geom_tile(aes(fill=Intensity), color="white")+
    scale_fill_gradientn(colors=color_blind)+
    xlab("Cluster")+
    ylab("Marker")+
    ggtitle("\"Differentially Expressed\" Cluster Medians Per Volunteer in VAC69a")+
    theme(#axis.text.x = element_text(hjust = 1),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Volunteer, ncol=5, labeller=as_labeller(labs))
)



# unique(long_deg_medians_all$Volunteer)


# just makin a table

whatevs <- select(data, colnames(data[3:12]))
significant <- whatevs[c(6,11,12,18,24,28,33,34),]

baseline_abundance <- significant[,1:5]
infected_abundance <- significant[,6:10]

base_table <- as.data.frame(
  cbind(
    apply(baseline_abundance*1000, 1, mean),
    apply(baseline_abundance*1000, 1, sd),
    rownames(baseline_abundance)
    ))

colnames(base_table) <- c("Mean", "STDV", "Cluster")

inf_table <- as.data.frame(cbind(
  apply(infected_abundance*1000, 1, mean),
  apply(infected_abundance*1000, 1, sd),
  rownames(infected_abundance)
))

colnames(inf_table) <- c("Mean", "STDV", "Cluster")

base_table$Timepoint <- "Baseline"
inf_table$Timepoint <- "T+6"

# combine table, convert level to float (stupid r...)
one_table <- rbind(base_table, inf_table)
one_table$Mean <- round(as.numeric(as.matrix(one_table)[,1]), digits=2)
one_table$STDV <- round(as.numeric(as.matrix(one_table)[,2]), digits=2)

#tables <- c(grid.table(base_table),grid.table(inf_table))

grid.table(one_table)

ft <- flextable(one_table)


(figura <- grid.arrange(tableGrob(base_table[,1:2]),
                        tableGrob(inf_table[,1:2]),
             big_heatmap_arsinh,
             widths=c(5,5),
             heights=c(5, 10),
             layout_matrix = rbind(c(1, 2),
                                   c(3, 3)),
             ncol=2,
             nrow=2)
)






