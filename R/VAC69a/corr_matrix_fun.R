# SPECIFIC ORDER OR NOT???
# 


cluster_correlation_matrix <- function(daf, da_analysis, tp, correlation, distance){
  
          all_frequencies <- data.frame(topTable(da_analysis, all=T, show_counts = T))
          
          all_frequencies <- gather(all_frequencies, sample_id, count, colnames(all_frequencies)[4:ncol(all_frequencies)])
          all_frequencies$volunteer <- stringr::str_match(all_frequencies$sample_id, "V[0-9]*")[, 1]
          all_frequencies$timepoint <- substr(all_frequencies$sample_id, 12,nchar(all_frequencies$sample_id))
          
          all_frequencies$sample_id <- gsub("counts_", "", all_frequencies$sample_id)
          counts <- n_cells(daf)
          all_frequencies$frequency <- all_frequencies$count / counts[all_frequencies$sample_id] *100
          
          baseline_freq_matrix <- all_frequencies %>%
            dplyr::filter(timepoint==paste(tp)) %>%
            select(cluster_id, volunteer, timepoint, frequency)
          
          baseline_freq_matrix <- spread(baseline_freq_matrix, cluster_id, frequency)
          baseline_spearman <- cor(baseline_freq_matrix[,3:ncol(baseline_freq_matrix)], method = correlation)
          
          baseline_dist <- dist(baseline_spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
          baseline_hclust <- hclust(baseline_dist)
          
          #check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
          baseline_spearman_df  <- data.frame(baseline_spearman, check.names = FALSE)
          baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)
          
          long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])
          
          corr_matrix_theme <-
            theme(axis.title = element_blank(),
                  axis.text.x = element_text(angle = 60, hjust = 1),
                  plot.title = element_text(hjust=0.5),
                  axis.text = element_text(size=7))
          
          colnames(baseline_spearman)[baseline_hclust$order]
          
          # ggplot(long_baseline_spearman, aes(x=factor(long_baseline_spearman$cluster_id_x, levels = specific_order), y=factor(long_baseline_spearman$cluster_id_y, levels=specific_order)))+
          ggplot(long_baseline_spearman, aes(x=factor(long_baseline_spearman$cluster_id_x, levels = colnames(baseline_spearman)[baseline_hclust$order]), y=factor(long_baseline_spearman$cluster_id_y, levels=colnames(baseline_spearman)[baseline_hclust$order])))+
            geom_tile(aes(fill=long_baseline_spearman$ro))+
            scale_fill_viridis(option="A")+
            labs(fill = expression(rho))+
            ggtitle(paste(tp, " (", correlation, ", ", distance,")", sep=''))+
            corr_matrix_theme
          
          }

base_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "Baseline", "pearson", "euclidean")
dod_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "DoD", "pearson", "euclidean")
t6_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "T6", "pearson", "euclidean")

lgd <- get_legend(base_corr)

pearson <- plot_grid(base_corr+theme(legend.position = "none"), dod_corr+theme(legend.position = "none"), t6_corr+theme(legend.position = "none"),lgd, ncol=4, rel_widths = c(5,5,5,1))

ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/tcell_cluster_correlation_matrix/pearson_all.png", pearson, height = 7, width=20)

ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/tcell_cluster_correlation_matrix/base_corr.png", base_corr)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/tcell_cluster_correlation_matrix/dod_corr.png", dod_corr)
ggsave("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/tcell_cluster_correlation_matrix/t6_corr.png", t6_corr)


#### loop for trying different parameters

methd <- c("euclidean", "maximum", "manhattan", "canberra", "minkowski")

mthds <- list()

for (i in methd){
  
  unit <- dist(baseline_spearman, method = i, diag = FALSE, upper = FALSE, p = 2)
  units <- hclust(unit)
  mthds[[paste(i)]] <- list(units$order)
  # colnames(mthds) <- c(colnames(mthds), paste(i))
}

orders <- c(mthds$euclidean, mthds$maximum, mthds$manhattan, mthds$canberra, mthds$minkowski)
#baselin_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
baseline_spearman_df  <- data.frame(baseline_spearman, check.names = FALSE)
baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)

long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])

# hclust_levels <- colnames(baseline_spearman)[baselin_hclust$order]

corr_matrix_theme <-
  theme(axis.title = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none")



for (i in 1:5){
  assign(paste("pearson", names(mthds)[i], sep='_' ), ggplot(long_baseline_spearman, aes_(x=factor(long_baseline_spearman$cluster_id_x, levels = colnames(baseline_spearman)[orders[[i]]]), y=factor(long_baseline_spearman$cluster_id_y, levels=colnames(baseline_spearman)[orders[[i]]])))+
           geom_tile(aes(fill=long_baseline_spearman$ro))+
           scale_fill_viridis(option="A")+
           ggtitle("Baseline")+
           labs(fill = expression(paste("Pearson ", rho)))+
           corr_matrix_theme)
}

spearman_methods <- plot_grid(spearman_canberra, spearman_euclidean, spearman_manhattan, spearman_maximum, spearman_minkowski, ncol=3)
pearson_methods <- plot_grid(pearson_canberra, pearson_euclidean, pearson_manhattan, pearson_maximum, pearson_minkowski, ncol=3)


jah_plots <- c("spearman_canberra", "spearman_euclidean", "spearman_manhattan", "spearman_maximum", "spearman_minkowski", "pearson_canberra", "pearson_euclidean", "pearson_manhattan", "pearson_maximum", "pearson_minkowski")

for (i in jah_plots){ggsave(paste("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/", i, ".png", sep=''), get(i))}

