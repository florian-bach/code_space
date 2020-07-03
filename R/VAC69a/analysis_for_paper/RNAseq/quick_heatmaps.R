library(ggplot2)
library(data.table)
library(dplyr)
library(vac69a.cytof)

setwd("/home/flobuntu/PhD/RNAseq/vac69a/all/xls")


all_unique_genes <- fread("all_unique_genes_cleaned.csv", header = T, stringsAsFactors = F)

list_of_all_unique <- split(all_unique_genes, all_unique_genes$file_name)
list_of_all_sig_unique <- lapply(list_of_all_unique, function(x)filter(x, padj<0.05))

C14_Baseline_sig <- list_of_all_sig_unique[[1]]
C56_Baseline_sig <- list_of_all_sig_unique[[2]]
DoD_Baseline_sig <- list_of_all_sig_unique[[3]]
T6_Baseline_sig <- list_of_all_sig_unique[[4]]
T6_DoD_sig <- list_of_all_sig_unique[[5]]

C14_Baseline_all <- list_of_all_unique[[1]]
C56_Baseline_all <- list_of_all_unique[[2]]
DoD_Baseline_all <- list_of_all_unique[[3]]
T6_Baseline_all <- list_of_all_unique[[4]]
T6_DoD_all <- list_of_all_unique[[5]]


# simple heatmaps ####

# capitalisation or lack thereof here will be reflected in the heatmap title
search_terms <- c("Cytokine", "Interferon", "Interleukin", "Tumor Necrosis Factor", "TNF", "Apoptosis", "Infection", "C-C")
  
search_results <- search_genes(search_terms)
search_results$`Tumor Necrosis Factor` <- c(search_results$`Tumor Necrosis Factor`, search_results$TNF)
search_results$TNF <- NULL

search_results

cleaned_search_results <- subset(search_results, unname(lapply(search_results, length)) != 0 )
  
Inflammation_markers <- list(`Inflammation Markers` = c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
                                                       "TICAM1", "TICAM2", "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
                                                       "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
                                                       "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "SECTM1", "ICAM1",
                                                       "CD40", "PDCD1", "CD274", "PDCD1LG2"))  

More_Inflammation_Markers <- list(`More Inflammation Markers` = c("CXCL11", "CXCL10", "CCL2", "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN",
                                                                  "TNF", "IL15", "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18",
                                                                  "IL7", "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
                                                                  "CXCL8")
)


Even_More_Inflammation_Markers <- list(`Even More Inflammation Markers` = c("IL31RA", "IL15RA", "IFNGR2", "CCR2", "TNFRSF1A", "IFNGR1", "CSF1R",
                                                                            "TNFRSF1B", "CXCR1", "CXCR3")
)


quick_gene_heatmaps(cleaned_search_results)
quick_gene_heatmaps(Inflammation_markers)
quick_gene_heatmaps(More_Inflammation_Markers)
quick_gene_heatmaps(Even_More_Inflammation_Markers)


  
# MOST SIGNIFICANT GENES HEATMAPS ####
  
big_table <- fread("all_unique_genes_cleaned.csv", header = TRUE, stringsAsFactors = TRUE)


most_sig <-  big_table %>%
  group_by(file_name) %>%
  top_n(n = -50, wt = padj) %>%
  select(Symbol, file_name)


list_of_top50 <- split(most_sig, most_sig$file_name)
list_of_top50_genes <- lapply(list_of_top50, function(x) as.character(x$Symbol))

quick_gene_heatmaps(list_of_top50_genes)  



absolute_most_sig <-  big_table %>%
  top_n(n = -50, wt = padj) %>%
  select(Symbol, file_name)

abs_sig <- list(`Most Significant Genes in Dataset`=as.character(absolute_most_sig$Symbol))

quick_gene_heatmaps(abs_sig)

  