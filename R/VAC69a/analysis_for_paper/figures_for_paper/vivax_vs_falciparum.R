# dod hist plot ####


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
library(ggplot2)
library(ggrepel)


# DoD sheet
dod_data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_dod_all.csv", header = T, stringsAsFactors = F)
dod_data <- subset(dod_data, grepl("5", dod_data$GOLevels))




# positive means enriched in vivax, negative means enriched in falciparum
dod_data$Cluster_Difference <- dod_data$`%Genes Cluster #2`-dod_data$`%Genes Cluster #1`



threshold <- 30

#DoD
dod_falci_rich_data <- subset(dod_data, dod_data$Cluster_Difference< -threshold)
dod_vivax_rich_data <- subset(dod_data, dod_data$Cluster_Difference> threshold)

vivax_x_limits <- c(80, NA)
vivax_y_limits <- c(20, NA)

dod_dot_plot <- ggplot(dod_data, aes(x=`%Genes Cluster #2`, y=`%Genes Cluster #1`))+
  theme_minimal()+
  scale_color_gradient2(high="#fec200", low="#db0085", midpoint = 0)+
  scale_y_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  scale_x_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
  ylab(expression('GO term enrichement'~italic("P. falciparum")~'at Diagnosis'))+
  xlab(expression('GO term enrichement'~italic("P. vivax")~'at Diagnosis'))+
  ggitle("Diagnosis")+
  ggforce::geom_circle(aes(x0=50, y0=50, r=sqrt((0.5*threshold)^2+(0.5*threshold)^2)), fill="grey", color=NA, alpha=0.2, inherit.aes = F)+
  # geom_rect(aes(ymin=50-(threshold*0.5), ymax=50+(threshold*0.5),
  #               xmin=50-(threshold*0.5), xmax=50+(threshold*0.5)),fill="grey", color=NA, alpha=0.2, inherit.aes = F )+
  geom_point(aes(color=Cluster_Difference))+
  # geom_label_repel(data=dod_vivax_rich_data, aes(label = stringr::str_wrap(GOTerm, 25)), size=1.8,
  #                  box.padding   = 0.35, nudge_x=-20, nudge_y=-30,
  #                  point.padding = 0.5, segment.alpha = 0.2, ylim  = vivax_y_limits, xlim  = vivax_x_limits)+
  theme(legend.position = "none",
        axis.text = element_text(size=6, ),
        axis.title = element_text(size=8),
        plot.margin=unit(c(0.5,0,0.5,0.5),"cm"))

dod_hist_plot <- ggplot(dod_data, aes(y=`%Genes Cluster #1`))+
  theme_minimal()+
  xlab("# GO Terms")+
  scale_y_continuous(breaks = seq(0,100, by=20), limits=c(-5, 105))+
  scale_x_continuous(breaks = seq(10,60, by=20))+
  geom_histogram(aes(fill = ..y..), orientation = "y", binwidth = 1, color="darkgrey", size=0.15)+
  scale_fill_gradient2(low="#fec200", high="#db0085", midpoint = 50)+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=6),
        legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0.3),"cm"))


vivax_falciparum_dod <- cowplot::plot_grid(dod_dot_plot, dod_hist_plot, ncol=2, rel_widths = c(3,1), align="h", axis="b")
ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_dod_GO.png", vivax_falciparum_dod, height = 3, width=3.9)


#ggsave("~/PhD/RNAseq/vac69a/all/xls/figures/vivax_falciparum_dod_GO_var.png", vivax_falciparum_dod, height = 7, width=9)









  # T6 sheet
t6_data <- data.table::fread("~/PhD/RNAseq/vac69a/cytoscape/vivax_falciparum_t6_all.csv", header = T, stringsAsFactors = F)

#t6_data <- subset(t6_data, grepl("5", t6_data$GOLevels))

# positive means enriched in vivax, negative means enriched in falciparum
t6_data$Cluster_Difference <- t6_data$`%Genes Cluster #2`-t6_data$`%Genes Cluster #1`

threshold <- 30

  
  #T6 for few labels on dot plot
  t6_falci_rich_data <- subset(t6_data, t6_data$Cluster_Difference< -50)
  t6_vivax_rich_data <- subset(t6_data, t6_data$Cluster_Difference> 50)
  
  
  t6_vivax_x_limits <- c(NA, NA)
  t6_vivax_y_limits <- c(NA, NA)
  
  t6_dot_plot <- ggplot(t6_data, aes(x=`%Genes Cluster #2`, y=`%Genes Cluster #1`))+
    theme_minimal()+
    #coord_cartesian(ylim = c(0,100), expand=F, xlim = c(0, 100))+
    scale_color_gradient2(high="#fec200", low="#db0085", midpoint = 0)+
    scale_y_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
    scale_x_continuous(breaks=seq(0,100,by=10), labels = scales::label_number(suffix="%"), limits = c(-5, 105), expand=c(0,0))+
    ylab(expression('GO term enrichement'~italic("P. falciparum")~'at T6'))+
    xlab(expression('GO term enrichement'~italic("P. vivax")~'at T6'))+
    ggitle("T6")+
    ggforce::geom_circle(aes(x0=50, y0=50, r=sqrt((0.5*threshold)^2+(0.5*threshold)^2)), fill="grey", color=NA, alpha=0.2, inherit.aes = F)+
    geom_point(aes(color=Cluster_Difference))+
    theme(legend.position = "none",
          axis.text = element_text(size=6, ),
          axis.title = element_text(size=8),
          plot.margin=unit(c(0.5,0,0.5,0.5),"cm"))
    # geom_label_repel(data=t6_vivax_rich_data, aes(label = stringr::str_wrap(GOTerm, 25)), size=2.8,
    #                  box.padding   = 0.35,
    #                  point.padding = 0.5, nudge_x = 40, nudge_y = 40, segment.alpha = 0.2)+
    # geom_label_repel(data=t6_falci_rich_data, aes(label = stringr::str_wrap(GOTerm, 20)), ylim  = t6_vivax_y_limits, xlim  = t6_vivax_x_limits, size=2.8,
    #                  box.padding   = 0.35,
    #                  point.padding = 0.5, nudge_x = -10, segment.alpha = 0.2)+

  
  t6_hist_plot <- ggplot(t6_data, aes(y=`%Genes Cluster #1`))+
    theme_minimal()+
    xlab("# GO Terms")+
    scale_y_continuous(breaks = seq(0,100, by=20), limits=c(-5, 105))+
    scale_x_continuous(breaks = seq(10,60, by=20))+
    geom_histogram(aes(fill = ..y..), orientation = "y", binwidth = 1, color="darkgrey", size=0.15)+
    scale_fill_gradient2(low="#fec200", high="#db0085", midpoint = 50)+
    theme(axis.title.y =  element_blank(),
          axis.text.y =  element_blank(),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size=6),
          legend.position = "none",
          plot.margin=unit(c(0.5,0.5,0.5,0.3),"cm"))
  
  
vivax_falciparum_t6 <- cowplot::plot_grid(t6_dot_plot, t6_hist_plot, ncol=2, rel_widths = c(3,1), align = "h", axis="b")
# ggsave("~/PhD/RNAseq/vac69a/all/xls/figures/vivax_falciparum_dod_GO.png", vivax_falciparum_dod, height = 7, width=9)


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vivax_falciparum_t6_GO_var.png", vivax_falciparum_t6, height = 3, width=3.9)








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
falci_t6_faves <- list("T Cell Genes of Interest"=c("CCR5", "CD19", "CD28", "CD38", "CD79A", "CD79B", "CTLA4", "CX3CR1", "CXCL1", "CXCL8", "CXCR3",
"CXCR4", "CXCR6", "ICOS", "IFNG", "IL12RB1", "IL12RB2", "IL21", "IL32", "IRF8", "LAG3", "TBX21", "TNFRSF13B",
"TNFRSF13C")  )

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


# vivax vs falciparum gene heatmaps ####
library(ComplexHeatmap)
library(tidyr)
library(dplyr)


falci_dod_data <- read.csv("~/PhD/RNAseq/vac63c/FirstDoD_genelist_padj_0.05.csv", header = T, stringsAsFactors = F, row.names = 1)
falci_t6_data <- read.csv("~/PhD/RNAseq/vac63c/T6vBas_First_DEGs_log2fc_cutoff.csv", header = T, stringsAsFactors = F, row.names = 1)

falci_dod_data$file_name <- "DoD_Baseline"
falci_t6_data$file_name <- "T6_Baseline"

falci_data <- falci_t6_data

falci_t6_faves <- list("T Cell Genes of Interest"=c("CCR5", "CD19", "CD20", "CD28", "IgD", "CD27", "CD38", "CD79A", "CXCR5", "CD79B", "CTLA4", "CX3CR1", "CXCR3",
                                                    "CXCR4", "CXCR6", "ICOS", "IFNG", "IL12RB1", "IL12RB2", "IL21", "IL32","MKI67", "LAG3", "TBX21"))

slimmed_falci_data <- subset(falci_data, falci_data$Symbol %in% falci_t6_faves$`T Cell Genes of Interest`)

slimmed_falci_data <- slimmed_falci_data[order(slimmed_falci_data$log2FoldChange,decreasing = TRUE),]



plot_data <- slimmed_falci_data

# plot_levels <- plot_data %>%
#   dplyr::filter(file_name=="T6_Baseline") %>%
#   dplyr::arrange(log2FoldChange) %>%
#   dplyr::select(Symbol)

definitive_falci_faves <- list("definitive_falci_faves"=plot_data$Symbol)
# 
# plot_plot <- ggplot(plot_data, aes(x=factor(file_name, levels=c("DoD_Baseline", "T6_Baseline")),
#                                    #y=factor(Symbol, levels = c(vivax_plot_levels$Symbol))))+
#                                    y=factor(Symbol, levels = rev(Symbol))))+
#   geom_tile(aes(fill=log2FoldChange, width=0.92, height=0.92), size=0.4,
#             #color=ifelse(plot_data$padj<0.05, "darkgreen", "black")
#             )+
#   scale_fill_gradientn(name="log2FC",
#                        values = scales::rescale(c(min(plot_data$log2FoldChange, na.rm = TRUE), 0, max(plot_data$log2FoldChange, na.rm = TRUE)), to=c(0,1)),
#                        colors = c("#0859C6","black","#FFA500"))+
#   theme_void()+
#   ggtitle("T cell genes at falciparum T6\n")+
#   guides(fill=guide_colorbar(nbin=30,
#                              breaks=seq(min(plot_data$log2FoldChange),max(plot_data$log2FoldChange), by=1)))+
#   coord_fixed(ratio = 1)+
#   theme(
#     axis.text.x = element_text(angle=45, hjust = 1, vjust = 1, size = 9),
#     axis.text.y = element_text(),
#     plot.title = element_text(hjust=0.5),
#     legend.title = element_text(),
#     legend.margin=margin(0,0,0,0),
#     legend.position = "right",
#     legend.box.margin=margin(0,0,0,0))
# 
# ggsave(paste("~/PhD/RNAseq/vac63c/T cell Genes of Interest falciparum.png"), height=6,width=6, plot_plot)



# combo plot




falci_faves_in_vivax <- falci_faves_in_vivax %>%
  select(colnames(plot_data)[-7]) %>%
  filter(Symbol %in% plot_data$Symbol)

combo_data <- rbind(cbind(plot_data[,-7],"Infection"="P. falciparum"), cbind(falci_faves_in_vivax, "Infection"="P. vivax"))

plot_plot2 <- ggplot(combo_data, aes(x=factor(Infection, levels=c("P. vivax", "P. falciparum")),
                                    #y=factor(Symbol, levels = c(vivax_plot_levels$Symbol))))+
                                    y=factor(Symbol, levels=rev(plot_data$Symbol))))+
  geom_tile(aes(fill=log2FoldChange, width=0.92, height=0.92), size=0.4,
            #color=ifelse(combo_data$padj<0.05, "red", "black")
  )+
  scale_fill_gradientn(name="log2FC",
                       values = scales::rescale(c(min(combo_data$log2FoldChange, na.rm = TRUE), 0, max(combo_data$log2FoldChange, na.rm = TRUE)), to=c(0,1)),
                       colors = c("#0859C6","black","#FFA500"))+
  theme_void()+
  ggtitle("T cell genes at T6\nrelative to Baseline\n")+
  guides(fill=guide_colorbar(nbin=30,
                             breaks=seq(min(combo_data$log2FoldChange),max(combo_data$log2FoldChange), by=1)))+
  coord_fixed(ratio = 1)+
  #facet_grid(~Infection)+
  theme(
    axis.text.x = element_text(angle=45, hjust = 1, vjust = 1, size = 7, face = "italic"),
    axis.text.y = element_text(size=9),
    plot.title = element_text(hjust=0.5),
    legend.title = element_text(hjust=0, size=8),
    legend.margin=margin(0,0,0,0),
    legend.position = "right",
    legend.box.margin=margin(0,0,0,0))



ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/falci_vivax_rna_t6.pdf", plot_plot2, height = 4.5, width=1.7)



gene_matrix <- select(combo_data, Symbol, log2FoldChange, Infection)
gene_matrix <- data.frame(tidyr::pivot_wider(gene_matrix, names_from = Infection, values_from = log2FoldChange))

rownames(gene_matrix) <- gene_matrix$Symbol

gene_matrix$Symbol <- NULL

gene_matrix <- as.matrix(gene_matrix[,c(2,1)])
colnames(gene_matrix) <-c("P. vivax", "P. falciparum")


col_fun_rna <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), c("#0859C6", "black", "#FFA500"))











gene_heatmap <- Heatmap(matrix = gene_matrix,
        cluster_rows = FALSE,
        show_heatmap_legend = TRUE,
        name = "log2FC",
        cluster_columns = FALSE,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        row_names_side = "left",
        col = col_fun_rna,
        column_names_rot = 45)

pdf("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/falci_vivax_rna_t6.pdf", height = 4.5, width=1.6)
draw(gene_heatmap
)
dev.off()


