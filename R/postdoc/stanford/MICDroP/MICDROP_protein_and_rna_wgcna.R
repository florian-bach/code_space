library(WGCNA)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# prepping the data for WGCNA ####
caro_data <- readRDS("/Users/fbach/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Immunology/RNA-seq/clean_data/MICDROP_dataObject.pf.geno.pheno.rds")
dim(caro_data$pheno)
dim(caro_data$expr)


#select most variable genes
keep <- rowMeans(caro_data$expr) > 0.1   # or 25th percentile
datExpr_f <- caro_data$expr[keep, ]

mad_values <- apply(datExpr_f, 1, mad)

cutoff <- quantile(mad_values, 0.8)
#0.8~4000 genes; 0.9~2000; 0.95~1000
most_variable_genes <- names(mad_values)[mad_values >= cutoff]


long_expr <- data.frame(caro_data$expr)%>%
  mutate(gene=rownames(.))%>%
  filter(gene %in% most_variable_genes)%>%
  pivot_longer(cols = names(data.frame(caro_data$expr)), names_to = "novogene code", values_to = "reads")%>%
  pivot_wider(names_from = "gene", values_from = reads)

expr_frame <- caro_data$pheno%>%
  select(id, mstatus, ageinwks, timepoint_num, `novogene code`)%>%
  left_join(., long_expr, by="novogene code")%>%
  pivot_longer(cols = names(long_expr)[2:length(names(long_expr))], names_to = "targetName", values_to = "RNA")%>%
  mutate(target_name=paste("rna_", targetName, sep=""),
         sample_id=paste(id, timepoint_num, sep="_"),
         expression=RNA)%>%
  select(sample_id, target_name, expression)


nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

slim_nulisa_data <- nulisa_data %>%
  mutate(target_name=paste("protein_", targetName, sep=""),
         expression=conc)%>%
  filter(timepoint_num!=24)%>%
  mutate(sample_id=paste(id, timepoint_num, sep="_"))%>%
  select(sample_id, target_name, expression)

combo_data <- bind_rows(slim_nulisa_data, expr_frame)

expression.data <- combo_data%>%
  pivot_wider(names_from = target_name, values_from = expression)%>%
  tibble::column_to_rownames("sample_id")





# running WGCNA ####

# get rid of non-expression data
# look at dendrogram of smaples to look for outliers
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# set parameters for network construction #### 
spt <- pickSoftThreshold(expression.data) 

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")


softPower <- 5

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
MEs$sample_id <- rownames(MEs)

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity


merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs


# merge to metadata ####
musical_metadata <- caro_data$pheno%>%
  mutate(sample_id=paste(id, timepoint_num, sep="_"))

traitRows <- musical_metadata%>%
  right_join(., MEs, by="sample_id")

datTraits <- traitRows%>%
  select(temp, log_qpcr, ageinwks, treatmentarm)

nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation


# figures ####

##  heatmap of correlations ####

module_cor_df <- data.frame(module.trait.correlation,
                            "value"="R",
                            "network"=rownames(module.trait.correlation))%>%
  pivot_longer(cols = c(temp, log_qpcr, ageinwks), names_to = "covariate", values_to = "metric")

module_p_df <- data.frame(module.trait.Pvalue,
                          "value"="p",
                          "network"=rownames(module.trait.Pvalue))%>%
  pivot_longer(cols=c(temp, log_qpcr, ageinwks, treatmentarm), names_to = "covariate", values_to = "metric")

heatmap_df <- bind_rows(module_cor_df, module_p_df)%>%
  pivot_wider(names_from = value, values_from = metric)%>%
  mutate(covariate=case_match(covariate, "ageinwks"~"age", "log_qpcr"~"parasitaemia", "temp"~"body temperature", .default = covariate))

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





## timepoint vs treatmentarm ####
micdrop_with_networks <- traitRows%>%
  pivot_longer(cols = starts_with("ME"), names_to = "network", values_to = "network_value")

wilcox_treatment <- micdrop_with_networks%>%
  filter(!is.na(timepoint_num))%>%
  pivot_wider(names_from = treatmentarm, values_from = network_value)%>%
  group_by(network, timepoint_num)%>%
  nest()%>%
  mutate(wilcox=map(data, ~wilcox.test(x = .$Placebo, y = .$`DP 1 year`)))%>%
  mutate(p=map_dbl(wilcox, ~.$p.value))%>%
  group_by(timepoint_num)%>%
  mutate(padj=p.adjust(p, method="fdr"))


(micdrop_with_networks%>%
    filter(!is.na(treatmentarm), network %in% c("MEblack", "MElightgreen"))%>%
    ggplot(., aes(x=factor(timepoint_num), y=network_value, fill=treatmentarm))+
    geom_boxplot(outliers = F)+
    ggpubr::stat_compare_means(label = "p.signif", hide.ns = T)+
    facet_wrap(~network, scales="free", nrow=1)+
    ylab("module eigengene")+
    scale_fill_manual(values=unname(unlist(time_cols)))+
    theme_minimal()
)


kl## visualise "loadings" ####
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
  pivot_longer(cols =starts_with("ME"), names_to = "network", values_to = "p")%>%
  group_by(network)%>%
  mutate(padj=p.adjust(p, method="BH"))


long_kme_p <- kME_df%>%
  select(-sample_id)%>%
  mutate(targetName=rownames(.))%>%
  pivot_longer(cols = starts_with("ME"), names_to = "network", values_to = "Pearson R")%>%
  left_join(., kME_pvalues, by=c("targetName", "network"))

(
  curve_plot1 <- long_kme_p%>%
    # filter(network%in%fever_networks)%>%
    filter(network %in% c("MEblack", "MElightgreen"))%>%
    # mutate(network=case_match(network, "MEturquoise"~"fever",
    #                           "MEred"~"parasite density",
    #                           "MEmagenta"~"parasite clearance", 
    #                           "MEbrown"~"combined fever and parasites"))%>%
    ggplot(., aes(x=`Pearson R`, y=-log10(padj), color=network))+
    geom_point(aes(color=network))+
    facet_wrap(~network, scales="free")+
    ggrepel::geom_text_repel(
      data = . %>% 
        # filter(network%in%fever_networks)%>%
        group_by(network)%>%
        slice_min(n = 15, order_by = padj),
      aes(label = targetName),   # or whatever your gene column is called
      force = 10,
      nudge_x = -0.4, box.padding = 0.2
    ) +
    scale_color_manual(values=c("black", "darkgreen"))+
    theme_minimal(base_size = 18)+
    theme(legend.position = "none"))

# simple correlations ####

coexpression_purf <- combo_data%>%
  group_by(targetName, timepoint_num)%>%
  nest()%>%
  mutate(spearman=map(data, ~cor.test(x = .$RNA, y = .$protein, method="spearman")))%>%
  mutate(rho=map_dbl(spearman, ~.$estimate))%>%
  mutate(p=map_dbl(spearman, ~.$p.value))%>%
  group_by(timepoint_num)%>%
  mutate(padj=p.adjust(p, method="fdr"))

top9_genes <- coexpression_purf%>%
  ungroup()%>%
  slice_min(n = 15, order_by = padj, with_ties = F)%>%
  pull(targetName)


combo_data%>%
  # filter(targetName%in%c("IL10", "IL15", "TLR3"))%>%
  filter(targetName%in%top9_genes)%>%
  
  ggplot(., aes(x=RNA, y=protein, color=factor(timepoint_num)))+
  geom_smooth(method="lm", alpha = 0.15)+
  ggpubr::stat_cor(method = "spearman")+
  geom_point()+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()
