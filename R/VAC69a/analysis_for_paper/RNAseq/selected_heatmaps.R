library(ComplexHeatmap)
library(tidyr)
library(dplyr)


vivax_t6_data <- read.csv("~/PhD/RNAseq/vac69a/all/xls/all_unique_genes_cleaned.csv")

vivax_t6_data$Infection <- "P. vivax"


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



tp_keep <- "T6_Baseline";chosen_list <- diana_list
tp_keep <- "DoD_Baseline";chosen_list <- inflammatory_list



faves_in_vivax <- vivax_t6_data %>%
  filter(file_name %in% tp_keep) %>%
  #filter(file_name %in% "DoD_Baseline") %>%
  #filter(padj<=0.05)%>%
  #filter(abs(log2FoldChange)>=log2(1.5))%>%
  select(Symbol, log2FoldChange, padj, file_name) %>%
  arrange(log2FoldChange)
  








combo_data <- subset(faves_in_vivax, faves_in_vivax$Symbol %in% c(chosen_list))

combo_data$log2FoldChange <- ifelse(combo_data$padj>0.05, 0, combo_data$log2FoldChange)


gene_matrix <- data.frame(combo_data %>%
  select(Symbol, file_name, log2FoldChange) %>%
  pivot_wider(names_from = file_name, values_from = log2FoldChange)) %>%
  arrange(tp_keep)

gene_matrix <- gene_matrix[gene_matrix[,2]!=0,]

#gene_matrix <- gene_matrix[ifelse(gene_matrix$T6_Baseline==0&gene_matrix$DoD_Baseline==0, FALSE, TRUE), ]

rownames(gene_matrix) <- gene_matrix$Symbol

gene_matrix$Symbol <- NULL

gene_matrix <- t(gene_matrix)


inferno <- c(colorspace::sequential_hcl("inferno", n=8))

maxfc <- max(gene_matrix)

col_fun_rna <- circlize::colorRamp2(breaks=seq(from = min(gene_matrix),
                                                 to = maxfc,
                                                 by = maxfc/(length(inferno)-2)),
                                    colors=c("#0859C6", inferno[1:7]))
                                    
col_fun_rna2 <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), colors = c("#0859C6", "black", "#FFA500"))


#col_fun_rna <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), c("#0859C6", "black", "#FFA500"))

rownames(gene_matrix) <- gsub("_", " relative to\n ", rownames(gene_matrix))
rownames(gene_matrix) <- gsub("DoD", "Diagnosis", rownames(gene_matrix))


gene_heatmap <- Heatmap(matrix = gene_matrix,
                        cluster_rows = FALSE,
                        show_heatmap_legend = TRUE,
                        # column_title ="Selected Differentially Expressed Cell Cycle Genes at T6",
                        column_title ="Selected Differentially Expressed Inflammatory Genes at Diagnosis",
                        
                        heatmap_legend_param = list(title = "log2FC"),
                        cluster_columns = FALSE,
                        column_names_gp = gpar(fontsize = 8),
                        
                        #row_names_gp = gpar(fontsize = 8, just="center"),
                        show_row_names = FALSE,
                        col = col_fun_rna,
                        column_names_rot = 45)

pdf("/home/flobuntu/PhD/figures_for_thesis/vivax_rna_dod_heat.pdf", height = 1.6, width=8)
draw(gene_heatmap, padding=unit(c(2,6,2,2), "mm")
)
dev.off()



png("/home/flobuntu/PhD/figures_for_thesis/vivax_rna_dod_heat.png", height = 1.6, width=8, res = 900, units = "in");draw(gene_heatmap);dev.off()

