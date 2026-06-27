library(WGCNA)
library(flashClust)
library(tidyr)
library(dplyr)
library(ggplot2)

`%notin%`=Negate(`%in%`)
time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")

infectiontype_cols <- list("A"="darkslateblue",
                           "S" = "darkred")


# read in RNA data and select variable genes ####
transcriptomic_data <- readRDS("~/Downloads/dataObject_MI_part1_part2.rds")

rna_expression_data <- data.frame(transcriptomic_data$expr)%>%
  mutate(targetName=rownames(transcriptomic_data$expr))%>%
  pivot_longer(cols = colnames(transcriptomic_data$expr), names_to = "rna_sample", values_to = "expression")

#select most variable genes
keep <- rowMeans(transcriptomic_data$expr) > 0.1   # or 25th percentile
datExpr_f <- transcriptomic_data$expr[keep, ]

mad_values <- apply(datExpr_f, 1, mad)

cutoff <- quantile(mad_values, 0.7)
#0.8~4000 genes; 0.9~2000; 0.95~1000
most_variable_genes <- names(mad_values)[mad_values >= cutoff]

# put data together ####
#transforming the data.frame so columns now represent genes and rows represent samples
slim_rna_expression_data <- rna_expression_data%>%
  filter(targetName %in% most_variable_genes)%>%
  mutate(target_name=paste("rna_", targetName, sep=""))
  

# read in protein data ####
protein_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

# put data together ####
#transforming the data.frame so columns now represent genes and rows represent samples
slim_protein_data<- protein_data%>%
  #these samples are being removed because the look like outliers in the network;
  # we are removing them here because the code below strops the expression.data object
  # of its rownames, which we need, in order to link the network data back to metadata
  filter(targetName !="CTSS",
         # timepoint!="day7",
         sample_id%notin% c("496 day7 NM", "363 day14 A", "496 day0 S"))%>%
  mutate(sample_id2=paste(id, timepoint_imm, infectiontype))%>%
  mutate(expression=z_conc)%>%
  mutate(targetName=paste("protein_", targetName, sep=""))%>%
  select(sample_id2, targetName, expression)#%>%
  #pivot_wider(names_from = targetName, values_from = z_conc)

# create map that links RNA sample names with proteomics sample names
rna_dictionary <- transcriptomic_data$pheno%>%
  mutate(rna_sample=`novogene code`)%>%
  mutate(sample_id2=paste(cohortid, timepoint, infection))%>%
  select(rna_sample, sample_id2)

slim_rna_expression_data2 <- slim_rna_expression_data%>%
  left_join(rna_dictionary, by="rna_sample")

# combine ####

long_expression.data <- slim_protein_data%>%
  filter(sample_id2 %in% rna_dictionary$sample_id2)%>%
  bind_rows(., slim_rna_expression_data2)%>%
  select(sample_id2, targetName, expression)%>%
  filter(sample_id2 %in% slim_protein_data$sample_id2)
  

expression.data <- as.data.frame(long_expression.data%>%
  pivot_wider(names_from = targetName, values_from = expression))

# expression.data <- expression.data[,!grepl("^protein", colnames(expression.data))]



# get rid of non-expression data
rownames(expression.data) <- expression.data$sample_id2
expression.data <- expression.data[,-1]


# get rid of non-expression data
# look at dendrogram of smaples to look for outliers
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#167 looks weird so let's cut tree at height 22
# 
# #Setting the graphical parameters
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# #draw on line to show cutoff height
# abline(h = 22, col = "red");
# 
# #cut the tree, remove outlier from dataset
# cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 22) #returns numeric vector
# expression.data <- expression.data[cut.sampleTree==1, ]


# set parameters for network construction #### 
spt <- pickSoftThreshold(expression.data) 

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")


softPower <- 4

adjacency <- adjacency(expression.data, power = softPower)

TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM

geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

# define modules ####
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(Modules)

#assigns each module number a color
ModuleColors <- labels2colors(Modules) 

table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)


#A ME (Module Eigengene) is the standardized gene expression profile for a given module.
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
MEs$sample_id2 <- rownames(MEs)

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity




merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs


# merge to metadata ####
musical_metadata <- transcriptomic_data$pheno%>%
  mutate(sample_id2=paste(cohortid, timepoint, infection))%>%
  mutate(log_qpcr = log10(qpcr+0.001))

traitRows <- musical_metadata%>%
  right_join(., MEs, by="sample_id2")

datTraits <- traitRows%>%
  select(temperature, log_qpcr, ageyrs)

nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation


# figures ####

##  heatmap of correlations ####

module_cor_df <- data.frame(module.trait.correlation,
                            "value"="R",
                            "network"=rownames(module.trait.correlation))%>%
  pivot_longer(cols = c(temperature, log_qpcr, ageyrs), names_to = "covariate", values_to = "metric")

module_p_df <- data.frame(module.trait.Pvalue,
                          "value"="p",
                          "network"=rownames(module.trait.Pvalue))%>%
  pivot_longer(cols=c(temperature, log_qpcr, ageyrs), names_to = "covariate", values_to = "metric")

heatmap_df <- bind_rows(module_cor_df, module_p_df)%>%
  pivot_wider(names_from = value, values_from = metric)%>%
  mutate(covariate=case_match(covariate, "ageyrs"~"age", "log_qpcr"~"parasitaemia", .default = covariate))

(module_cor_heatmap <- heatmap_df%>%
  # mutate(network=case_match(network, "MEturquoise"~"module 2", "MEred"~"module 1", "MEthistle3"~"module 3", "MEwhite"~"module 4"))%>%
  ggplot(., aes(x=covariate, y=network, fill=R))+
  # geom_point(aes(size=(-log10(p)), color=R))+
  geom_tile()+
  geom_text(color=ifelse(heatmap_df$R>=0.2, "black", "white"),
            aes(label=paste(round(R, digits = 2), "\n(p=",
                            signif(p, 1), ")", sep = "")))+
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#000000", "#FFA500", "#FFFFB2"),  # blue → black → yellow
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = c(-1, 1)
  )+
  scale_color_manual(values=c("black","white"))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.title = element_text(hjust = 0.16)))


# module eigengenes through time ####

parasite_networks <- c("MEblue", "MEgrey", "MEsalmon", "MEred")
fever_networks <- c("MEgreen", "MEyellow", "MEturqoise")


musical_with_networks <- traitRows%>%
  pivot_longer(cols = colnames(traitRows)[188:ncol(traitRows)], names_to = "network", values_to = "network_value")


## timepoint vs network value ####
(musical_with_networks%>%
   filter(timepoint_category%in%c("Baseline", "Day 0","Day 14"))%>%
   filter(infection%in%c("A", "S"), network%in%c(parasite_networks, fever_networks))%>%
   mutate(timepoint=factor(timepoint, levels=c("Baseline", "Day 0","Day 14")))%>%
   # mutate(network=case_match(network, "MEturquoise"~"fever",
   #                           "MEred"~"parasite density",
   #                           "MEmagenta"~"parasite clearance", 
   #                           "MEbrown"~"combined fever and parasites"))%>%
   ggplot(., aes(x=timepoint_category, y=network_value, fill=timepoint_category))+
   geom_violin(draw_quantiles = 0.5)+
   geom_smooth(method="lm")+
   geom_point()+
   ggpubr::stat_compare_means(aes(group=timepoint_category), comparisons = list(c("Day 0", "Baseline")),
                              label = "p.signif", hide.ns = F)+
   facet_wrap(~infection+network, scales="free_x", nrow=2)+
   ylab("module eigengene")+
   scale_fill_manual(values=unname(unlist(time_cols)))+
   theme_minimal()
)

## infectiontype vs network values####
musical_with_networks%>%
  filter(timepoint_category%in%c("Baseline"))%>%
  filter(infection%in%c("A", "S"), network%in%c(parasite_networks, fever_networks))%>%
  mutate(timepoint_category=factor(timepoint_category, levels=c("Baseline")))%>%
  ggplot(., aes(x=infection, y=network_value, fill=infection))+
  geom_violin(draw_quantiles = 0.5)+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_compare_means(label = "p.signif", hide.ns = F)+
  facet_wrap(~timepoint_category+network, scales="free_x")+
  ylab("Module Eigengene")+
  scale_fill_manual(values=infectiontype_cols)+
  theme_minimal()

## timepoint_category vs network value parasite controllers and non controllers####
(asymp_ctrl <- musical_with_networks%>%
   filter(infection%in%c("A", "S"))%>%
   mutate(day14_para=if_else(timepoint_category=="Day 14"&infection=="A" & qpcr > 10, "parasitemic_Day 14", "no_parasites_Day 14"))%>%
   group_by(cohortid)%>%
   mutate(class2= if_else(any(day14_para=="parasitemic_Day 14"), "parasitaemic at Day 14", "parasites cleared"))%>%
   filter(timepoint_category%in%c("Baseline", "Day 0","Day 14"))%>%
   mutate(timepoint_category=factor(timepoint_category, levels=c("Baseline", "Day 0","Day 14")))%>%
   filter(!is.na(class2))%>%
   # mutate(network=case_match(network, "MEturquoise"~"fever",
   #                           "MEred"~"parasite density",
   #                           "MEmagenta"~"parasite clearance", 
   #                           "MEbrown"~"combined fever and parasites"))%>%
   ggplot(., aes(x=timepoint_category, y=network_value, color=timepoint_category, fill=class2))+
   geom_violin(draw_quantiles = 0.5)+
   geom_point(aes(color=timepoint_category), position = position_dodge(width=0.9))+
   ggpubr::stat_compare_means(label = "p.signif", hide.ns = T, size=13, vjust=1)+
   facet_wrap(~infection+network, scales="free_x")+
   scale_color_manual(values=unname(unlist(time_cols)))+
   scale_fill_manual(values=c("#636363", "darkgrey"))+
   ylab("module eigengene")+
   theme_minimal()+
   guides(fill=guide_legend(override.aes = list(color=NA)), color=guide_none())+
   theme(axis.title.x = element_blank(),
         legend.title = element_blank(),
         legend.position = "bottom"))
# ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/asymp_ctrl.png", asymp_ctrl, width = 12, height=8, dpi=444)

(symp_time <- musical_with_networks%>%
    
    # filter(infection%in%c("S"), network%in%c(fever_networks))%>%
    filter(infection%in%c("S", "A"), network%in%c(fever_networks))%>%
    
    mutate(day14_para=if_else(timepoint_category=="Day 14"&infection=="A" & qpcr > 10, "parasitemic_Day 14", "no_parasites_Day 14"))%>%
    group_by(cohortid)%>%
    mutate(class2=if_else(any(day14_para=="parasitemic_Day 14"), "parasitaemic at Day 14", "parasites cleared"))%>%
    filter(timepoint_category%in%c("Baseline", "Day 0","Day 14"))%>%
    mutate(timepoint_category=factor(timepoint_category, levels=c("Baseline", "Day 0","Day 14")))%>%
    filter(!is.na(class2))%>%
    ggplot(., aes(x=timepoint_category, y=network_value, fill=timepoint_category))+
    geom_violin(draw_quantiles = 0.5)+
    geom_point(position = position_dodge(width=0.9))+
    ggpubr::stat_compare_means(vjust=2,comparisons = list(
      c("Day 0", "Baseline")),
      label = "p.signif", hide.ns = F)+
    facet_wrap(~infection+network, scales="free_x")+
    scale_color_manual(values=unname(unlist(time_cols)))+
    scale_fill_manual(values=unname(unlist(time_cols)))+
    ylab("Module Eigengene")+
    theme_minimal(base_size = 20)+
    guides(fill=guide_legend(override.aes = list(color=NA)), color=guide_none())+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          legend.title = element_blank()))




## visualise "loadings" ####
# Compute module membership (kME)
kME <- cor(expression.data, MEs, use = "p")
# Convert to data frame for inspection
kME_df <- as.data.frame(kME)
long_kME_df <- kME_df%>%
  mutate(targetName=rownames(.))%>%
  pivot_longer(cols = grep("^ME", colnames(.), value=T), names_to = "network", values_to = "network_value")


# now calculate correlations between networks and individual proteins 


# Compute module membership (kME)
kME <- cor(expression.data, MEs, use = "p")

# Convert to data frame for inspection
kME_df <- as.data.frame(kME)

corPvalueStudent <- function(cor, nSamples) {
  t = cor * sqrt((nSamples - 2) / (1 - cor^2))
  2 * pt(-abs(t), df = nSamples - 2)
}
kME_pvalues <- corPvalueStudent(kME, nSamples)

kME_pvalues <- data.frame(kME_pvalues)%>%
  mutate(targetName=rownames(.))%>%
  select(-sample_id)%>%
  pivot_longer(cols = colnames(traitRows)[188:ncol(traitRows)], names_to = "network", values_to = "p")%>%
  group_by(network)%>%
  mutate(padj=p.adjust(p, method="BH"))


long_kme_p <- kME_df%>%
  select(-sample_id)%>%
  mutate(targetName=rownames(.))%>%
  pivot_longer(cols = colnames(traitRows)[188:ncol(traitRows)], names_to = "network", values_to = "Pearson R")%>%
  left_join(., kME_pvalues, by=c("targetName", "network"))

(
  curve_plot1 <- long_kme_p%>%
    # filter(network%in%fever_networks)%>%
    filter(network %in% c("MEblack"))%>%
    # mutate(network=case_match(network, "MEturquoise"~"fever",
    #                           "MEred"~"parasite density",
    #                           "MEmagenta"~"parasite clearance", 
    #                           "MEbrown"~"combined fever and parasites"))%>%
    ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network))+
    geom_point(aes(color=network))+
    facet_wrap(~network)+
    ggrepel::geom_text_repel(
      data = . %>% 
        # filter(network%in%fever_networks)%>%
        group_by(network)%>%
        slice_min(n = 50, order_by = padj),
      aes(label = targetName),   # or whatever your gene column is called
      force = 10,
      nudge_x = -0.4, box.padding = 0.2
    ) +
    # scale_color_manual(values=c("tomato4", "turquoise4"))+
    # theme_minimal(base_size = 18)+
    theme(legend.position = "none"))


(
  curve_plot2 <- long_kme_p%>%
    # filter(network%in%fever_networks)%>%
    filter(network %in% fever_networks)%>%
    # mutate(network=case_match(network, "MEturquoise"~"fever",
    #                           "MEred"~"parasite density",
    #                           "MEmagenta"~"parasite clearance", 
    #                           "MEbrown"~"combined fever and parasites"))%>%
    ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network))+
    geom_point(aes(color=network))+
    facet_wrap(~network)+
    ggrepel::geom_text_repel(
      data = . %>% 
        # filter(network%in%fever_networks)%>%
        group_by(network)%>%
        slice_min(n = 40, order_by = padj),
      aes(label = targetName),   # or whatever your gene column is called
      force = 10,
      nudge_x = -0.4, box.padding = 0.2
    ) +
    # scale_color_manual(values=c("tomato4", "turquoise4"))+
    # theme_minimal(base_size = 18)+
    theme(legend.position = "none"))


