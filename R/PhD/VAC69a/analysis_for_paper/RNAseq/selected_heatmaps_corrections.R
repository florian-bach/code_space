library(ComplexHeatmap)
library(tidyr)
library(dplyr)


vivax_t6_data <- read.csv("~/PhD/RNAseq/vac69a/all/xls/all_unique_genes_cleaned.csv")

vivax_t6_data$Infection <- "P. vivax"

# 
diana_list <- scan("~/PhD/RNAseq/vac69a/all/xls/gene_lists/diana_cell_cycle_genes.txt", what="NULL", sep = ",")
# 
# 
# inflammatory_list <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
#                        "TICAM1", "TICAM2", "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
#                        "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
#                        "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "SECTM1", "ICAM1",
#                        "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
#                        "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
#                        "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7", 
#                        "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
#                        "CXCL8", "AIM2", "OASL", "MX1", "MX2", "DDX58", "IL18", "CXCL9", "IL6", "IL21")
# 
# #shorter
# inflammatory_list <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
#                        "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
#                        "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
#                        "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "ICAM1",
#                        "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
#                        "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
#                        "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7", 
#                        "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
#                        "CXCL8", "AIM2", "OASL", "MX1", "MX2", "DDX58", "IL18", "CXCL9", "IL6", "IL21")

inflammatory_list <- c("IL1RN","IL1B",  "TNF", "IL15", "IL27", "STAT1", "STAT2", "GBP1", "GBP2", "OASL", "MX1", "MX2",
                      "IRF1", "IRF9", "MYD88", "DDX58", "SOD2", "HIF1A", "CD274", "PDCD1LG2", "CASP1", "CASP3","CASP8",
                      "BAK1", "BCL2L1", "BCL2L13", "ADAM17", "CCR1", "CCR2", "CCL2", "CCL3", #"CCL4", "CCL5",
                      "CSF1", "CXCL10", "IDO1", "TNFa", "TLR", "NLRP", "MDA5", "MAVS", "AIM2")

diana_list <- unlist(strsplit(c("POLA1, MKI67, POLA2, POLD1, POLD3, POLE2, PCNA, FANCE, FANCG, FANCI, FANCL, RAD51, RAD51AP1, RAD54L, CNE2, CDC25A, CDC6, CCNA2, CCNB1, CCNB2, CDC25B, CDC25C, CENPA, CENPE, CENPF, CENPI"),
                              split = ", ",))
                     

tp_keep <- "T6_Baseline";chosen_list <- diana_list
tp_keep <- "DoD_Baseline";chosen_list <- inflammatory_list



faves_in_vivax <- vivax_t6_data %>%
  #filter(file_name %in% tp_keep) %>%
  filter(file_name == "T6_Baseline") %>%
  filter(padj<=0.05)%>%
  filter(abs(log2FoldChange)>=log2(1.5))%>%
  select(Symbol, log2FoldChange, padj, file_name) %>%
  arrange(log2FoldChange)
  


chosen_list <- c(diana_list,
                 #grep("*cdc*", faves_in_vivax$Symbol, value = TRUE, ignore.case = TRUE),
                 diana_list2)





combo_data <- subset(faves_in_vivax, faves_in_vivax$Symbol %in% c(chosen_list))

combo_data$log2FoldChange <- ifelse(combo_data$padj>0.05, 0, combo_data$log2FoldChange)


gene_matrix <- data.frame(combo_data %>%
  arrange(desc(log2FoldChange)) %>%                          
  select(Symbol, file_name, log2FoldChange) %>%
  filter(log2FoldChange>log2(1.5))%>%
  pivot_wider(names_from = file_name, values_from = log2FoldChange))

gene_matrix <- gene_matrix[gene_matrix[,2]!=0,]

#gene_matrix <- gene_matrix[ifelse(gene_matrix$T6_Baseline==0&gene_matrix$DoD_Baseline==0, FALSE, TRUE), ]

rownames(gene_matrix) <- gene_matrix$Symbol

gene_matrix$Symbol <- NULL

gene_matrix <- as.matrix(gene_matrix)


inferno <- c(colorspace::sequential_hcl("inferno", n=8))

maxfc <- max(gene_matrix)

col_fun_rna_t6 <- circlize::colorRamp2(breaks=seq(from = min(gene_matrix),
                                                 to = maxfc,
                                                 by = maxfc/length(inferno)),
                                                  colors=inferno[1:7])
                                    

col_fun_rna_dod <- circlize::colorRamp2(breaks=seq(from = min(gene_matrix),
                                   to = maxfc,
                                   by = maxfc/(length(inferno))),
                                  colors=inferno[1:7])



#col_fun_rna2 <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), colors = c("#0859C6", "black", "#FFA500"))


col_fun_rna <- circlize::colorRamp2(c(min(gene_matrix), 0, max(gene_matrix)), c("#0859C6", "black", "#FFA500"))

# rownames(gene_matrix) <- gsub("_", " relative to\n ", rownames(gene_matrix))
# rownames(gene_matrix) <- gsub("DoD", "Diagnosis", rownames(gene_matrix))

colnames(gene_matrix) <- gsub("_", " relative to\n ", colnames(gene_matrix))

gene_heatmap <- Heatmap(matrix = gene_matrix,
                        cluster_rows = FALSE,
                        show_heatmap_legend = TRUE,
                        column_title ="Selected Cell\nCycle Genes",
                        column_title_gp = gpar(just="center", fontsize=10),
                        cluster_columns = FALSE,
                        row_names_gp = gpar(fontsize = 10, fontface="italic"),
                        heatmap_legend_param = list(legend_position = "bottom",
                                                    col=col_fun_rna_dod,
                                                    title = "log2FC",
                                                    legend_direction = "horizontal",
                                                    title_position = "topcenter",
                                                    #legend_width = unit(6.2, "cm"),
                                                    border = FALSE),
                        column_names_gp = gpar(fontsize = 8, just="center"),
                        show_row_names = TRUE,
                        row_names_side = "left",
                        #col = col_fun_rna_t6,
                        col=col_fun_rna_dod,
                        column_names_rot = 45)



pdf("/home/flobuntu/PhD/manuscripts/vac69a/jci_corrections/vivax_t6_heat.pdf", width = 2, height=8.2)
draw(gene_heatmap,  heatmap_legend_side = "bottom",
     padding=unit(c(2,8,2,8), "mm")
)
dev.off()

#b, l, t, 

# 
# 
# 
# pdf("/home/flobuntu/PhD/figures_for_thesis/chapter_1/vivax_rna_dod_heat.pdf", height = 1.8, width=8.2)
# draw(gene_heatmap,  heatmap_legend_side = "bottom",
#      padding=unit(c(2,8,2,2), "mm")
# )
# dev.off()



