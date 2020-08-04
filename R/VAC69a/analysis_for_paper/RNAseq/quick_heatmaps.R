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



vivax_t6_faves <- list(`T6 Specific T cell Genes` = c("CD38", "CTLA4", "CXCR6", "GZMA", "ICOS", "IL21", "MKI67"))
quick_gene_heatmaps(vivax_t6_faves, sort_by = "T6_Baseline")

falci_t6_faves <- list(`T cell genes at vivax T6`=c("CCR5", "CD19", "CD28",
                                  "CD38", "CD79A", "CD79B", "CTLA4", "CX3CR1", "CXCL1", "CXCL8", "CXCR3","CXCR4", "CXCR6", "ICOS", "IFNG",
                                  "IL12RB1", "IL12RB2", "IL21", "IL32", "IRF8", "LAG3", "TBX21", "TNFRSF13B","TNFRSF13C"))


ifnas <- list(`Interferon Alpha Genes`=paste0("IFNA", c(1,2,4,5,6,7,8,10,13,14,16,21), sep=""))

ifnbs <- list(`Interferon Beta Genes`=paste0("IFNB", seq(1,100), sep=""))


quick_gene_heatmaps(cleaned_search_results, sort_by = "DoD_Baseline")
quick_gene_heatmaps(Inflammation_markers, sort_by = "DoD_Baseline")
quick_gene_heatmaps(More_Inflammation_Markers, sort_by = "DoD_Baseline")
quick_gene_heatmaps(Even_More_Inflammation_Markers, sort_by = "DoD_Baseline")
quick_gene_heatmaps(ifnas, sort_by = "DoD_Baseline")
quick_gene_heatmaps(ifnbs, sort_by = "DoD_Baseline")
quick_gene_heatmaps(falci_t6_faves, sort_by = "T6_Baseline")
  
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

  