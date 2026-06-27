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
genes_to_keep <- rowMeans(transcriptomic_data$expr) > 0.1   # or 25th percentile
# datExpr_f <- transcriptomic_data$expr[genes_to_keep, ]

samples_to_keep <- transcriptomic_data$pheno%>%
  filter(infection=="A" | infection=="asymptomatic" | timepoint_category=="Baseline")%>%
  pull(`novogene code`)

datExpr_f <- transcriptomic_data$expr[genes_to_keep, samples_to_keep]

mad_values <- apply(datExpr_f, 1, mad)

cutoff <- quantile(mad_values, 0.85)
#0.8~4000 genes; 0.9~2000; 0.95~1000
most_variable_genes <- names(mad_values)[mad_values >= cutoff]

# put data together ####
#transforming the data.frame so columns now represent genes and rows represent samples
wide_rna_expression_data <- data.frame(rna_expression_data%>%
                                         filter(targetName %in% most_variable_genes, rna_sample%in% samples_to_keep, rna_sample!="B209")%>%
                                         pivot_wider(names_from = "targetName", values_from = "expression"))

rownames(wide_rna_expression_data) <- wide_rna_expression_data$rna_sample
expression.data <- wide_rna_expression_data[,-match("rna_sample", colnames(wide_rna_expression_data))]




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

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")


softPower <- 6
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

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

png("~/postdoc/stanford/abstracts/immunology_retreat_2025/wcgna_dendro.png", height = 3, width=6, units = "in", res = 444)
par(cex = 5) 
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#A ME (Module Eigengene) is the standardized gene expression profile for a given module.
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
MEs$sample_id <- rownames(MEs)

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

# construct a cluster tree. Any branches below 0.25 are more than 75% related,
# thus will be merged; none were merged here
# METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
# par(mar = c(0,4,2,0)) #seting margin sizes
# par(cex = 0.6);#scaling the graphic
# plot(METree)
# abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
# merge to metadata ####

musical_metadata <- transcriptomic_data$pheno%>%
  mutate(sample_id=`novogene code`)%>%
  mutate(log_qpcr = log10(qpcr+0.001))

traitRows <- musical_metadata%>%
  inner_join(., MEs, by="sample_id")

datTraits <- traitRows%>%
  select(temperature, log_qpcr, ageyrs)

nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation


#Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
# par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 1,
               plotLegend = FALSE,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
temperature = as.data.frame(datTraits$temperature)
names(temperature) = "temperature"

log_qpcr = as.data.frame(datTraits$log_qpcr)
names(log_qpcr) = "log_qpcr"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p", method="spearman"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, temperature, use = "p", method="spearman"))
temp_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(temperature), sep="")
names(temp_GSPvalue) = paste("p.GS.", names(temperature), sep="")
temp_GSPvalue$targetName <- rownames(temp_GSPvalue)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, log_qpcr, use = "p"))
qpcr_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(log_qpcr), sep="")
names(qpcr_GSPvalue) = paste("p.GS.", names(log_qpcr), sep="")
qpcr_GSPvalue$targetName <- rownames(qpcr_GSPvalue)
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


module_cor_heatmap <- heatmap_df%>%
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
        legend.title = element_text(hjust = 0.16))

# ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/module_cor_heatmap.png", module_cor_heatmap, width=5.5, height=4.7, dpi=444)


## integration with rest of data ####

parasite_networks <- c("MEblack", "MEgrey", "MEturquoise")
# fever_networks <- c("MEblue", "MEpink", "MEred")


musical_with_networks <- traitRows%>%
  pivot_longer(cols = colnames(traitRows)[188:ncol(traitRows)], names_to = "network", values_to = "network_value")


## timepoint vs network value ####
a_timepoints <- musical_with_networks%>%
  filter(timepoint_category%in%c("Baseline", "Day 0","Day 14", "Day 28"))%>%
  filter(infection%in%c("A"), network%in%parasite_networks)%>%
  mutate(timepoint=factor(timepoint, levels=c("Baseline", "Day 0","Day 14")))%>%
  # mutate(network=case_match(network, "MEturquoise"~"fever",
  #                           "MEred"~"parasite density",
  #                           "MEmagenta"~"parasite clearance", 
  #                           "MEbrown"~"combined fever and parasites"))%>%
  ggplot(., aes(x=timepoint_category, y=network_value, fill=timepoint_category))+
  geom_violin(draw_quantiles = 0.5)+
  geom_smooth(method="lm")+
  geom_point()+
  facet_wrap(~infection+network, scales="free_x", nrow=1)+
  ggpubr::stat_compare_means(aes(group=timepoint_category), comparisons = list(c("Baseline", "Day 0"),
                                                                               c("Baseline", "Day 14"),
                                                                               c("Baseline", "Day 28")),
                             label = "p.signif", hide.ns = F)+
  
  ylab("module eigengene")+
  scale_fill_manual(values=unname(unlist(time_cols)))+
  theme_minimal(base_size = 20)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

s_timepoints <- musical_with_networks%>%
  filter(timepoint_category%in%c("Baseline", "Day 0"))%>%
  filter(infection%in%c("S"), network%in%parasite_networks)%>%
  mutate(timepoint=factor(timepoint, levels=c("Baseline", "Day 0","Day 14")))%>%
  # mutate(network=case_match(network, "MEturquoise"~"fever",
  #                           "MEred"~"parasite density",
  #                           "MEmagenta"~"parasite clearance", 
  #                           "MEbrown"~"combined fever and parasites"))%>%
  ggplot(., aes(x=timepoint_category, y=network_value, fill=timepoint_category))+
  geom_violin(draw_quantiles = 0.5)+
  geom_smooth(method="lm")+
  geom_point()+
  facet_wrap(~infection+network, scales="free_x", nrow=1)+
  ggpubr::stat_compare_means(aes(group=timepoint_category), comparisons = list(c("Baseline", "Day 0")),
                             label = "p.signif", hide.ns = F)+
  
  ylab("module eigengene")+
  scale_fill_manual(values=unname(unlist(time_cols)))+
  theme_minimal(base_size = 20)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

combo_timepoints <- a_timepoints + s_timepoints + plot_layout(widths = c(1.5, 1))
ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/parasite_networks.png", combo_timepoints, width = 17.5, height=7, dpi=444)

ggsave()
## infectiontype vs network values####
musical_with_networks%>%
  filter(timepoint_category%in%c("Baseline"))%>%
  filter(infection%in%c("A", "S"))%>%
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
asymp_ctrl <- musical_with_networks%>%
  filter(infection%in%c("A"), network %in% c(parasite_networks, "MEred"))%>%
  mutate(day14_para=if_else(timepoint_category=="Day 14"&infection=="A" & qpcr > 10, "parasitemic_Day 14", "no_parasites_Day 14"))%>%
  group_by(cohortid)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_Day 14"), "parasitaemic at Day 14", "parasites cleared"))%>%
  filter(timepoint_category%in%c("Baseline", "Day 0","Day 14"))%>%
  mutate(timepoint_category=factor(timepoint_category, levels=c("Baseline", "Day 0","Day 14", "Day 28")))%>%
  filter(!is.na(class2))%>%
  # mutate(network=case_match(network, "MEturquoise"~"fever",
  #                           "MEred"~"parasite density",
  #                           "MEmagenta"~"parasite clearance", 
  #                           "MEbrown"~"combined fever and parasites"))%>%
  ggplot(., aes(x=timepoint_category, y=network_value, color=timepoint_category, fill=class2))+
  geom_violin(draw_quantiles = 0.5)+
  geom_point(aes(color=timepoint_category), position = position_dodge(width=0.9))+
  ggpubr::stat_compare_means(label = "p.signif", hide.ns = T, size=13, vjust=1)+
  facet_wrap(~infection+network, scales="free_x", nrow=1)+
  scale_color_manual(values=unname(unlist(time_cols)))+
  scale_fill_manual(values=c("#636363", "darkgrey"))+
  ylab("module eigengene")+
  theme_minimal()+
  guides(fill=guide_legend(override.aes = list(color=NA)), color=guide_none())+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/parasite_networks_parsite_control.png", asymp_ctrl, width = 10, height=5, dpi=444)

(symp_time <- musical_with_networks%>%
    
    # filter(infection%in%c("S"), network%in%c(fever_networks))%>%
    filter(infection%in%c("S"), network%in%c("MEblue", "MEpink"))%>%
    
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
ggsave("~/postdoc/stanford/abstracts/immunology_retreat_2025/asymp_ctrl.png", asymp_ctrl, width = 10, height=5, dpi=444)



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


curve_plot1 <- long_kme_p%>%
  # filter(network%in%fever_networks)%>%
  filter(network %in% c("MEgreen", parasite_networks))%>%
  mutate(network_color=gsub("ME", "", network))%>%
  # mutate(network=case_match(network, "MEturquoise"~"fever",
  #                           "MEred"~"parasite density",
  #                           "MEmagenta"~"parasite clearance", 
  #                           "MEbrown"~"combined fever and parasites"))%>%
  ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network_color))+
  geom_point()+
  facet_wrap(~network, nrow=1)+
  ggrepel::geom_text_repel(
    data = . %>% 
      # filter(network%in%fever_networks)%>%
      group_by(network)%>%
      slice_min(n = 15, order_by = padj),
    aes(label = targetName),   # or whatever your gene column is called
    force = 10,
    nudge_x = -0.4, box.padding = 0.2
  ) +
  theme_minimal()+
  scale_color_manual(values=c("black", "grey", "red", "turquoise"))+
  # theme_minimal(base_size = 18)+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/parasite_networks_curve_plot.png", curve_plot1, width = 10, height=5, dpi=444)


curve_plot2 <- long_kme_p%>%
  # filter(network%in%fever_networks)%>%
  filter(network %in% c("MEblue", "MEpink"))%>%
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
      slice_min(n = 30, order_by = padj),
    aes(label = targetName),   # or whatever your gene column is called
    force = 10,
    nudge_x = -0.4, box.padding = 0.2
  ) +
  scale_color_manual(values=c("blue", "orchid1"))+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/fever_networks_curve_plot.png", curve_plot2, width = 5, height=5, dpi=444)


module_cor_heatmap1 <- heatmap_df%>%
  filter(network %in% fever_networks)%>%
  mutate(text_color=ifelse(.$R>=0.2, "black", "white"))%>%
  ggplot(., aes(x=covariate, y=network, fill=R))+
  # geom_point(aes(size=(-log10(p)), color=R))+
  geom_tile()+
  geom_text(aes(color=text_color,
                label=paste(round(R, digits = 2), "\n(p=",
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
        legend.title = element_text(hjust = 0.16))

fever_combo <- module_cor_heatmap1+curve_plot1+patchwork::plot_layout(widths = c(1,2))




curve_plot2 <- long_kme_p%>%
  filter(network%in%parasite_networks)%>%
  ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network))+
  geom_point(aes(color=network))+
  facet_wrap(~network)+
  ggrepel::geom_text_repel(
    data = . %>% 
      filter(network%in%parasite_networks)%>%
      group_by(network)%>%
      slice_min(n = 20, order_by = padj),
    aes(label = targetName),   # or whatever your gene column is called
    force = 10,
    nudge_x = -0.4, box.padding = 0.2
  ) +
  # scale_color_manual(values=c("tomato4", "turquoise4"))+
  # theme_minimal(base_size = 18)+
  theme(legend.position = "none")

module_cor_heatmap2 <- heatmap_df%>%
  filter(network %in% parasite_networks)%>%
  mutate(text_color=ifelse(.$R>=0.2, "black", "white"))%>%
  ggplot(., aes(x=covariate, y=network, fill=R))+
  # geom_point(aes(size=(-log10(p)), color=R))+
  geom_tile()+
  geom_text(aes(color=text_color,
                label=paste(round(R, digits = 2), "\n(p=",
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
        legend.title = element_text(hjust = 0.16))

parasite_combo <- module_cor_heatmap2+curve_plot2+patchwork::plot_layout(widths = c(1,2))




module_cor_heatmap1 <- heatmap_df%>%
  filter(network %in% fever_networks)%>%
  mutate(text_color=ifelse(.$R>=0.2, "black", "white"))%>%
  ggplot(., aes(x=covariate, y=network, fill=R))+
  # geom_point(aes(size=(-log10(p)), color=R))+
  geom_tile()+
  geom_text(aes(color=text_color,
                label=paste(round(R, digits = 2), "\n(p=",
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
        legend.title = element_text(hjust = 0.16))

fever_combo <- module_cor_heatmap1+curve_plot1+patchwork::plot_layout(widths = c(1,2))




curve_plot2 <- long_kme_p%>%
  filter(network%in%parasite_networks)%>%
  ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network))+
  geom_point(aes(color=network))+
  facet_wrap(~network)+
  ggrepel::geom_text_repel(
    data = . %>% 
      filter(network%in%parasite_networks)%>%
      group_by(network)%>%
      slice_min(n = 20, order_by = padj),
    aes(label = targetName),   # or whatever your gene column is called
    force = 10,
    nudge_x = -0.4, box.padding = 0.2
  ) +
  # scale_color_manual(values=c("tomato4", "turquoise4"))+
  # theme_minimal(base_size = 18)+
  theme(legend.position = "none")

module_cor_heatmap2 <- heatmap_df%>%
  filter(network %in% parasite_networks)%>%
  mutate(text_color=ifelse(.$R>=0.2, "black", "white"))%>%
  ggplot(., aes(x=covariate, y=network, fill=R))+
  # geom_point(aes(size=(-log10(p)), color=R))+
  geom_tile()+
  geom_text(aes(color=text_color,
                label=paste(round(R, digits = 2), "\n(p=",
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
        legend.title = element_text(hjust = 0.16))

parasite_combo <- module_cor_heatmap2+curve_plot2+patchwork::plot_layout(widths = c(1,2))












blue_top_hits <- long_kme_p %>%
  filter(network=="MEturquoise")%>%
  slice_min(n = 25, order_by = padj)

brown_top_hits <- long_kme_p %>%
  filter(network=="MEred")%>%
  slice_min(n = 12, order_by = padj)



clean_data%>%
  filter(targetName %in% brown_top_hits$targetName)%>%
  filter(infection%in%c("A", "S"))%>%
  filter(timepoint%in%c("Baseline", "Day 0", "Day 14"))%>%
  ggplot(., aes(x=qpcr, y=concentration, color=infection))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  scale_x_log10()+
  scale_color_manual(values=infection_cols)+
  facet_wrap(~factor(targetName, levels=brown_top_hits$targetName), nrow=2)+
  theme_minimal()

clean_data%>%
  filter(targetName %in% blue_top_hits$targetName)%>%
  mutate(day14_para=if_else(timepoint=="Day 14"&infection=="A" & qpcr > 10, "parasitemic_Day 14", "no_parasites_Day 14"))%>%
  group_by(cohortid)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_Day 14"), "parasitaemic at Day 14", "parasites cleared"))%>%
  filter(infection%in%c("A"))%>%
  filter(timepoint%in%c("Baseline"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=class2))+
  geom_violin()+
  ggpubr::stat_compare_means(label = "p.signif")+
  # scale_color_manual(values=infection_cols)+
  facet_wrap(~factor(targetName, levels=blue_top_hits$targetName), nrow=5, scales="free")+
  theme_minimal()

clean_data%>%
  filter(targetName %in% brown_top_hits$targetName)%>%
  mutate(day14_para=if_else(timepoint=="Day 14"&infection=="A" & qpcr > 10, "parasitemic_Day 14", "no_parasites_Day 14"))%>%
  group_by(cohortid)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_Day 14"), "parasitaemic at Day 14", "parasites cleared"))%>%
  filter(infection%in%c("A"))%>%
  filter(timepoint%in%c("Baseline"))%>%
  ggplot(., aes(x=timepoint, y=concentration, fill=class2))+
  geom_violin()+
  ggpubr::stat_compare_means(label = "p.signif")+
  # scale_color_manual(values=infection_cols)+
  facet_wrap(~factor(targetName, levels=brown_top_hits$targetName), nrow=2, scales="free")+
  theme_minimal()



musical_with_networks%>%
  filter(network %in% c("MEbrown", "MEblue"))%>%
  ggplot(., aes(x=temperature, y=network_value))+
  geom_point()+
  ggpubr::stat_cor()+
  geom_smooth(method="lm")+
  facet_wrap(~network)+
  theme_minimal()

musical_with_networks%>%
  filter(infection%in%"A")%>%
  # filter(network %in% c("MEbrown", "MEturquoise"))%>%
  ggplot(., aes(x=qpcr+0.001, y=network_value))+
  geom_point()+
  scale_x_log10()+
  ggpubr::stat_cor()+
  geom_smooth(method="lm")+
  facet_wrap(~network)+
  theme_minimal()


musical_with_networks%>%
  filter(network %in% c("MEbrown", "MEblue"))%>%
  ggplot(., aes(x=ageyrs, y=network_value))+
  geom_point()+
  geom_smooth(method="lm")+
  ggpubr::stat_cor()+
  facet_wrap(~network)+
  theme_minimal()
# simple corr ####

View(clean_data%>%
       group_by(targetName)%>%
       nest()%>%
       mutate(pearson_cor = map(data, ~cor.test(.$log_qpcr, .$concentration, method="p")))%>%
       mutate(pearson_cor_p=map_dbl(pearson_cor, ~.$p.value),
              pearson_cor_r=map_dbl(pearson_cor, ~.$estimate))%>%
       arrange(pearson_cor_p))


# QC ####
# sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance
# 
# #Setting the graphical parameters
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# 
# #Plotting the cluster dendrogram
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# 
# #167 looks weird so let's cut tree at height 35
# 
# #Setting the graphical parameters
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# #draw on line to show cutoff height
# abline(h = 35, col = "red");
# sandbox ####

# 
# combined_gsp <- inner_join(qpcr_GSPvalue, temp_GSPvalue, by="targetName")
# colnames(combined_gsp) <- c("qpcr_p", "targetName", "temperature_p")
# 
# sig_targets <- combined_gsp%>%
#   filter(qpcr_p<0.05 | temperature_p<0.05)
# 
# qpcr_only_sigs <- sig_targets%>%
#   filter(qpcr_p<0.05 & temperature_p>0.05)
# 
# fever_only_sigs <- sig_targets%>%
#   filter(qpcr_p>0.05 & temperature_p<0.05)
# 
# 
# qpcr_only_highlights <- qpcr_only_sigs%>%
#   top_n(n = 9, wt = desc(qpcr_p))
# 
# fever_only_highlights <- fever_only_sigs%>%
#   top_n(n = 9, wt = desc(temperature_p))
# 
# combo_highlights <- sig_targets%>%
#   pivot_longer(cols = c(qpcr_p, temperature_p), names_to = "correlate", values_to = "p")%>%
#   group_by(targetName)%>%
#   top_n(n = 1, wt = desc(p))%>%
#   ungroup()%>%
#   top_n(n = 9, wt = desc(p))
# 
#   
# # 
# clean_data%>%
#   filter(targetName %in% qpcr_only_highlights$targetName)%>%
#   ggplot(., aes(x=qpcr+0.001, y=concentration, color=temperature))+
#   geom_point()+
#   scale_x_log10()+
#   geom_smooth(method="lm")+
#   ggtitle("qpcr highlights")+
#   ggpubr::stat_cor(method="spearman")+
#   facet_wrap(~targetName, scales="free")+
#   scale_color_viridis_c()+
#   theme_minimal()
# 
# clean_data%>%
#   filter(targetName %in% fever_only_highlights$targetName)%>%
#   ggplot(., aes(x=temperature, y=concentration, color=temperature))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   ggtitle("fever highlights")+
#   ggpubr::stat_cor(method="spearman")+
#   scale_color_viridis_c()+
#   facet_wrap(~targetName, scales="free")+
#   theme_minimal()
# 
# 
