plotClusterHeatmap(daf, hm2 = "abundances",
                   #m = "meta35",
                   k = "meta35",
                   cluster_anno = TRUE,
                   draw_freqs = TRUE,
                   scale=T)



ei <- metadata(daf)$experiment_info

da_formula1 <- createFormula(ei, cols_fixed = "timepoint", cols_random = "volunteer")
# 
# levels(ei$timepoint)

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

# contrast_baseline <- createContrast(c(1, rep(0, 3)))
 contrast_dod<- createContrast(c(0,1,0,0,0,0))
 contrast_t6 <- createContrast(c(0,0,1,0,0,0))
 contrast_c45 <- createContrast(c(rep(0, 3), 1, 0, 0))

# contrast_baseline <- createContrast(c(1, rep(0, 6)))
# contrast_dod<- createContrast(c(0,1,rep(0, 5)))
# contrast_t6 <- createContrast(c(0, 1, 0, rep(0, 4)))
# contrast_c45 <- createContrast(c(rep(0, 3), 1, rep(0, 3)))
# 
# contrast_t6 <- createContrast(c(0, 1, rep(0,3)))

data.frame(parameters = colnames(design), contrast_t6)


da_baseline <- diffcyt(daf,
                       #formula = da_formula1,
                       experiment_info = ei,
                       marker_info = panel,
                       contrast = contrast_baseline,
                       analysis_type = "DA",
                       # method_DA = "diffcyt-DA-GLMM",
                       method_DA = "diffcyt-DA-edgeR",
                       clustering_to_use = "meta35",
                       verbose = T)


da_dod <- diffcyt(daf,
                  design = design,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "som100",
                  verbose = F)

da_t6 <- diffcyt(daf,
                 #formula = da_formula1,
                 design = design,
                 contrast = contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "som100",
                 verbose = F)


da_inf <- diffcyt(daf,
                 formula = da_formula1,
                 design = design,
                 contrast = contrast_inf,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "meta35",
                 verbose = F)



da_c45 <- diffcyt(daf,
                  design = design,
                  contrast = contrast_c45,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta35",
                  verbose = F)

table(rowData(da_dod$res)$p_adj < FDR_cutoff)
table(rowData(da_t6$res)$p_adj < FDR_cutoff)
#table(rowData(da_inf$res)$p_adj < FDR_cutoff)


plotDiffHeatmap(daf, da_dod, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T)
plotDiffHeatmap(daf, da_c45, th = FDR_cutoff, normalize = TRUE, hm1 = T)




density_plot <- plotDR(daf, "UMAP", color_by="CCR7")+ facet_wrap(c("volunteer", "timepoint"))+
  stat_density2d(bins=100, size=0.11, colour="maroon")

density_plot$layers[[1]] <- NULL
density_plot <- density_plot+xlim(c(-15, 8))+ylim(c(-10,15))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank())






###########    glmm     #######



da_dod <- diffcyt(daf,
                  #design = design,
                  formula = da_formula1,
                  contrast = contrast_dod,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-GLMM",
                  clustering_to_use = "meta35",
                  verbose = F)


da_t6 <- diffcyt(daf,
                 formula = da_formula1,
                 #design = design,
                 contrast = contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-GLMM",
                 clustering_to_use = "meta35",
                 verbose = F)

da_inf <- diffcyt(daf,
                 formula = da_formula1,
                 #design = design,
                 contrast = contrast_t6,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-GLMM",
                 clustering_to_use = "meta35",
                 verbose = F)

table(rowData(da_inf$res)$p_adj < FDR_cutoff)
table(rowData(da_t6$res)$p_adj < FDR_cutoff)


plotDiffHeatmap(daf, da_inf, th = FDR_cutoff, normalize = TRUE, hm1 = T, all = T)
plotDiffHeatmap(daf, da_t6, th = FDR_cutoff, normalize = TRUE, hm1 = T, all = T)


topTable(da_t6, show_props = TRUE, format_vals = TRUE, digits = 2)
2