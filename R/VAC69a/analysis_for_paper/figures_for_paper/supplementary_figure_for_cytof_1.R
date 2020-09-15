# Panel A: boxplots of all cluster frequencies ####

merged_daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))

da_t6 <- diffcyt(merged_daf,
                 design = design,
                 contrast = pairwise_contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "flo_merge",
                 verbose = T)

all_cluster_freqs <- diffcyt_boxplot(da_t6, merged_daf, counts=F, FDR=1, logFC = 0)

ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_all_clusters_freqs.png", all_cluster_freqs, height = 12, width=18)# works

# Panel B: ds_limma results of differential marker expression on coarse clusters ####
library(CATALYST)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(diffcyt)

daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/most_coarse_merge.csv", header=T, stringsAsFactors = F)

merged_daf<- mergeClusters(daf, k = "mod_meta45", table = coarse_table, id = "coarse_merge")

ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "batch"))

levels(ei$timepoint)

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,1)))

states <- names(marker_classes(merged_daf))
states <- states[-c(1:21,
                    match(c("CCR7", "CD127", "CD20", "CXCR5", "CD103", "TCRgd", "TIM3", "IntegrinB7", "CD56", "CD3", "CD49d", "CD45RA", "CD4", "Vd2", "Va72", "CD161", "FoxP3", "CD45RO"), states),
                    53, 54, 56:length(states)
)
]

metadata(merged_daf)$id_state_markers <- states


ds_t6 <- diffcyt(merged_daf,                                            
                 design = design, contrast = pairwise_contrast_t6,                    
                 analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",         
                 clustering_to_use = "coarse_merge", verbose = TRUE)


res_DA_t6 <- topTable(ds_t6, all = TRUE, show_logFC = T)
table(res_DA_t6$p_adj <= 0.05)

sigs_t6 <- subset(res_DA_t6, p_adj<=0.05)

spider_sigs_t6 <- data.frame(sigs_t6$cluster_id, sigs_t6$marker_id, sigs_t6$logFC)
spider_sigs_t6$direction <- ifelse(sigs_t6$logFC>=0, "up", "down")

trimmed_spider_sigs_t6 <- subset(spider_sigs_t6, 2^spider_sigs_t6$sigs_t6.logFC>= 1.1 |  2^spider_sigs_t6$sigs_t6.logFC<=0.9)

(ds_limma_t6_heatmap <- ggplot(trimmed_spider_sigs_t6, aes(x= sigs_t6.cluster_id, y=sigs_t6.marker_id, label=round(2^sigs_t6.logFC, digits = 2)))+
    #geom_tile(aes(fill=rescale(sigs_t6.logFC, to=c(-5, 5))))+
    geom_tile(aes(fill=round(2^sigs_t6.logFC, digits = 2)))+
    geom_text(parse = TRUE, aes(color=2^sigs_t6.logFC>0.8))+
    scale_color_manual(values = c("TRUE"="black", "FALSE"="white"), guide=FALSE)+
    scale_fill_gradientn("Fold Change", values = scales::rescale(c(2^min(trimmed_spider_sigs_t6$sigs_t6.logFC), 1, 2^max(trimmed_spider_sigs_t6$sigs_t6.logFC)), to=c(0,1)), colors = c("#0859C6","white","#FFA500"))+
    scale_x_discrete(labels=c("gamma delta"=expression(paste(gamma, delta))))+
    theme_minimal()+
    ggtitle("Fold Change in Intensity of Marker Expression\nper Cluster at T6 relative to Baseline")+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=14, color = "black"),
          axis.text.x = element_text(hjust=1, angle = 60),
          legend.position = "right",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
)


ggsave("/home/flobuntu/PhD/cytof/vac69a/final_figures_for_paper/supp_ds_limma.png", ds_limma_t6_heatmap, width = 7, height = 7)
