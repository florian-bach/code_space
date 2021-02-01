library(ComplexHeatmap)
library(tidyr)
library(dplyr)


vivax_t6_data <- read.csv("~/PhD/RNAseq/vac69a/all/xls/all_unique_genes_cleaned.csv")

vivax_t6_data$Infection <- "P. vivax"


faves_in_vivax <- vivax_t6_data %>%
  filter(file_name %in% c("T6_Baseline", "DoD_Baseline")) %>%
  #filter(padj<=0.05)%>%
  #filter(abs(log2FoldChange)>=log2(1.5))%>%
  select(Symbol, log2FoldChange, padj, file_name) %>%
  arrange(log2FoldChange)
  


diana_list <- scan("~/PhD/RNAseq/vac69a/all/xls/gene_lists/diana_cell_cycle_genes.txt", what="NULL", sep = ",")


inflammatory_list <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
                       "TICAM1", "TICAM2", "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
                       "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
                       "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "SECTM1", "ICAM1",
                       "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
                       "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
                       "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7", 
                       "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
                       "CXCL8")





combo_data <- subset(faves_in_vivax, faves_in_vivax$Symbol %in% c(diana_list, inflammatory_list))

combo_data$log2FoldChange <- ifelse(combo_data$padj>0.05, 0, combo_data$log2FoldChange)


gene_matrix <- data.frame(combo_data %>%
  select(Symbol, file_name, log2FoldChange) %>%
  pivot_wider(names_from = file_name, values_from = log2FoldChange)) %>%
  arrange(DoD_Baseline)

gene_matrix <- gene_matrix[-c(50:100),]

gene_matrix <- gene_matrix[ifelse(gene_matrix$T6_Baseline==0&gene_matrix$DoD_Baseline==0, FALSE, TRUE), ]

rownames(gene_matrix) <- gene_matrix$Symbol

gene_matrix$Symbol <- NULL

gene_matrix <- t(gene_matrix)


inferno <- c(colorspace::sequential_hcl("inferno", n=8))

maxfc <- max(gene_matrix)

col_fun_rna <- circlize::colorRamp2(seq(from = min(gene_matrix),
                                        to = maxfc,
                                        by=maxfc/(length(inferno)-1)
                                        ),
                                       c("#0859C6", inferno))
                                    
col_fun_rna2 <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), c("#0859C6", "black", "#FFA500"))


#col_fun_rna <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), c("#0859C6", "black", "#FFA500"))



gene_heatmap <- Heatmap(matrix = gene_matrix,
                        cluster_rows = FALSE,
                        show_heatmap_legend = TRUE,
                        column_title ="Selected Cell Cycle Genes at T6",
                        heatmap_legend_param = list(title = "log2FC"),
                        cluster_columns = FALSE,
                        column_names_gp = gpar(fontsize = 8),
                        row_names_gp = gpar(fontsize = 8),
                        row_names_side = "left",
                        col = col_fun_rna,
                        column_names_rot = 45)

# pdf("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/falci_vivax_rna_t6.pdf", height = 1.6, width=4.5)
# draw(gene_heatmap
# )
# dev.off()
# 


png("/home/flobuntu/PhD/figures_for_thesis/vivax_rna_dod_t6_heat.png", height = 1.6, width=8, res = 900, units = "in");draw(gene_heatmap);dev.off()

