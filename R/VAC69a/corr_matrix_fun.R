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
          
          # could make a function that does all this..
          
          print(unique(baseline_freq_matrix$timepoint))
          
          baseline_freq_matrix <- spread(baseline_freq_matrix, cluster_id, frequency)
          baseline_spearman <- cor(baseline_freq_matrix[,3:ncol(baseline_freq_matrix)], method = correlation)
          
          baseline_dist <- dist(baseline_spearman, method = distance, diag = FALSE, upper = FALSE, p = 2)
          baseline_hclust <- hclust(baseline_dist)
          
          #check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
          baseline_spearman_df  <- data.frame(baseline_spearman, check.names = FALSE)
          baseline_spearman_df$cluster_id_x <- rownames(baseline_spearman_df)
          
          long_baseline_spearman <- gather(baseline_spearman_df, cluster_id_y, ro, colnames(baseline_spearman_df)[1:ncol(baseline_spearman_df)-1])
          
          hclust_levels <- colnames(baseline_spearman)[baseline_hclust$order]
          
          corr_matrix_theme <-
            theme(axis.title = element_blank(),
                  axis.text.x = element_text(angle = 60, hjust = 1),
                  plot.title = element_text(hjust=0.5))
          
          ggplot(long_baseline_spearman, aes(x=factor(long_baseline_spearman$cluster_id_x, levels = hclust_levels), y=factor(long_baseline_spearman$cluster_id_y, levels=hclust_levels)))+
            geom_tile(aes(fill=long_baseline_spearman$ro))+
            scale_fill_viridis(option="A")+
            labs(fill = expression(rho))+
            corr_matrix_theme
          
          }

base_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "Baseline", "pearson", "euclidean")
dod_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "DoD", "pearson", "euclidean")
t6_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "T6", "pearson", "euclidean")

pearson <- plot_grid(base_corr, dod_corr, t6_corr, ncol=3)

base_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "Baseline", "spearman", "euclidean")
dod_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "DoD", "spearman", "euclidean")
t6_corr <- cluster_correlation_matrix(merged_daf, da_baseline, "T6", "spearman", "euclidean")


spearman <- plot_grid(base_corr, dod_corr, t6_corr, ncol=3)

plot_grid(pearson, spearman, ncol=1)
