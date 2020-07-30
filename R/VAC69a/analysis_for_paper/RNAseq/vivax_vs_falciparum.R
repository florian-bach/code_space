# data <- read.csv("~/PhD/RNAseq/vac63c/FirstDoD_genelist_padj_0.05.csv", header = T, stringsAsFactors = F, row.names = 1)
# 
# data_up <- subset(data, data$log2FoldChange>0)
# 
# data_down <- subset(data, data$log2FoldChange<0)
# 
# write.table(data$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_all.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# write.table(data_up$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_up.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# write.table(data_down$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_down.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# 
# 

library(magrittr)

# DoD sheet
data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_dod_all.csv", header = T, stringsAsFactors = F)
# T6 sheet
data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_t6_all.csv", header = T, stringsAsFactors = F)

# positive means enriched in vivax, negative means enriched in falciparum
data$Cluster_Difference <- data$`%Genes Cluster #2`-data$`%Genes Cluster #1`


falci_rich <- subset(data, data$Cluster_Difference < -20)
vivax_diff <- subset(data, data$Cluster_Difference > 20)
shared <- subset(data, data$Cluster_Difference < 20 & data$Cluster_Difference > -20)

list_of_rich <- list("falci_rich"=falci_rich, "vivax_diff"=vivax_diff, "shared"=shared)

list_of_rich_genes <- lapply(list_of_rich, function(x)x$`Associated Genes Found`)


sapply(list_of_rich, nrow)

# DoD
# falci_rich vivax_diff     shared 
# 4            10             275 

# T6
# falci_rich vivax_diff     shared 
# 130          0        105 


#split gene list so that each gene becomes a list entry
list_unique_genes <- lapply(list_of_rich_genes, function(x)
  unlist(
    lapply(x, function(y) strsplit(y, ','))
    )
)

# get rid of special characters, order list and get rid of duplicats  
list_unique_genes <- lapply(list_unique_genes, function(x) x %<>%
                              gsub(" ", "", .) %>%
                              gsub("]", "", ., fixed = T) %>%
                              gsub("[", "", ., fixed = T) %>%
                              .[order(.)] %>%
                              unique(.))

sapply(list_unique_genes, length)

# DoD
# falci_rich vivax_diff     shared 
# 141          282         3019 

# T6
# falci_rich vivax_diff     shared 
# 861              0           395

#DoD genes of note:

falci_unique <- list_unique_genes$falci_rich[!list_unique_genes$falci_rich  %in% list_unique_genes$shared]

vivax_unique <- list_unique_genes$vivax_diff[list_unique_genes$falci_rich  %in% list_unique_genes$shared]

vivax_unique <- list(`Unique to Vivax DoD`=vivax_unique)

# 473 out of 861

falci_dod_faves <- list("BCL6", "BCR", "CCR2", "CGAS", "HIF1A", "HLA-F", "GATA3", "HLA-G", "IL7R", "JAK2", "JAK3",
                        "STAT1", "STAT3", "STAT5B", "TBX21")

vivax_dod_faves <- list("CCL5", "OAS2", "TNF", "TNFSF13")

# ALL falci_rich genes AND all vivax_rich genes are also SHARED
shared_dod_faves <- unlist(subset(falci_dod_faves, falci_dod_faves %in% list_unique_genes$shared))

# T6 genes of note:
# falciparum:
falci_t6_faves <- list("CCR5", "CD19", "CD28", "CD38", "CD79A", "CD79B", "CTLA4", "CX3CR1", "CXCL1", "CXCL8", "CXCR3",
"CXCR4", "CXCR6", "ICOS", "IFNG", "IL12RB1", "IL12RB2", "IL21", "IL32", "IRF8", "LAG3", "TBX21", "TNFRSF13B",
"TNFRSF13C")  

# shared
shared_faves <- unlist(subset(falci_faves, falci_faves %in% list_unique_genes$shared))
# "CD28"  "CD38"  "CXCL8" "IFNG"  "TBX21"

vivax_specific_dod <- DoD_Baseline[!scan("~/PhD/RNAseq/vac69a/all/xls/gene_lists/DoD_Baseline_ALL_significant_symbol_only.txt", what = "") %in% scan("~/PhD/RNAseq/vac63c/vac63a_b_sig_dod_baseline_all.txt", what=""),]


# up down plots ####
falciparum_up_down <- data.frame("Timepoint"=c("DoD", "DoD", "T6", "T6"), "Direction"=c(rep(c("up", "down"), times=2)))
falciparum_up_down$Genes <- c(1861, -1121, 721, -221)


falciparum_up_down_plot <- ggplot(falciparum_up_down, aes(x=Timepoint, y=Genes))+
  geom_bar(stat="identity", aes(fill=Timepoint))+
  theme_minimal()+
  ylim(-1200, 2000)+
  geom_text(aes(label=abs(Genes), vjust= -0.2), data = subset(falciparum_up_down, falciparum_up_down$Direction=="up"))+
  geom_text(aes(label=abs(Genes), vjust= 1.2), data = subset(falciparum_up_down, falciparum_up_down$Direction=="down"))+
  scale_fill_manual(values=colorspace::sequential_hcl(6, palette = "Purple Yellow")[c(3,4)])+
  theme(legend.position = "none",
      axis.title=element_blank(),
      axis.text.x = element_text(angle=45, hjust=1))

ggsave("./figures/falciparum_up_down_plot.png", falciparum_up_down_plot)  

combo_plot <- cowplot::plot_grid(sig_gene_count_plots, falciparum_up_down_plot, rel_widths = c(5, 2))
ggsave("./figures/combo_plot.png", combo_plot)  
