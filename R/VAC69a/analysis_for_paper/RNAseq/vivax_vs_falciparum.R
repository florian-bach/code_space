data <- read.csv("~/PhD/RNAseq/vac63c/FirstDoD_genelist_padj_0.05.csv", header = T, stringsAsFactors = F, row.names = 1)

data_up <- subset(data, data$log2FoldChange>0)

data_down <- subset(data, data$log2FoldChange<0)

write.table(data$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_all.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(data_up$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_up.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(data_down$Symbol, "~/PhD/RNAseq/vac63c/vac63c_sig_t6_baseline_down.txt", sep = "\t", quote = F, row.names = F, col.names = F)




data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_t6_all.csv", header = T, stringsAsFactors = F)

# positive means enriched in vivax, negative means enriched in falciparum
data$Cluster_Difference <- data$`%Genes Cluster #2`-data$`%Genes Cluster #1`

#data$Cluster_Difference <- data$X.Genes.Cluster..2-data$`%Genes Cluster #1`


diff_data <- subset(data, abs(data$Cluster_Difference)>20)


genes <- diff_data$`Associated Genes Found`
unique_genes <- unique(unlist(lapply(genes, function(x)as.list(unlist(strsplit(x, ','))))))

unique_genes <- gsub(" ", "", unique_genes)
unique_genes <- gsub("[", "", unique_genes, fixed = T)
unique_genes <- gsub("]", "", unique_genes, fixed = T)

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
