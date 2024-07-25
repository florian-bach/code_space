day14_only <- filterSCE(kids_only, timepoint=="day 14")

ei <- metadata(day14_only)$experiment_info

colnames(design)
# [1] "(Intercept)"      "timepointday 0"   "timepointday 14"  "timepointday 7"   "classsymptomatic" "subject_id161"   
# [7] "subject_id176"    "subject_id268"    "subject_id324"    "subject_id353"    "subject_id363"    "subject_id571"  
#timepoint=baseline is dummy; 
design <- createDesignMatrix(ei, c("class", "subject_id"))

FDR_cutoff <- 0.05

class_contrast <- createContrast(c(c(0, 1, 0, 0), rep(0,5)))

da_class_contrast <- diffcyt(day14_only,
                             design = design,
                             contrast = class_contrast,
                             analysis_type = "DA",
                             method_DA = "diffcyt-DA-edgeR",
                             clustering_to_use = "meta50",
                             verbose = T)

table(rowData(da_class_contrast$res)$p_adj < FDR_cutoff)

da <- rowData(da_class_contrast$res)
plotDiffHeatmap(day14_only, da, all=TRUE, )

vac69a.cytof::diffcyt_boxplot(da_class_contrast, day14_only, logFC = 0, FDR=1)


results_table <- data.frame(diffcyt::topTable(da_class_contrast, all=T, show_counts = TRUE, show_props = TRUE))

long_results <- results_table %>%
  pivot_longer(cols = starts_with("props"), names_to = "sample_prop", values_to = "freq")%>%
  pivot_longer(cols = starts_with("counts"), names_to = "sample_counts", values_to = "count")%>%
  mutate("subject_id"=substr(sample_prop, 20, 22),
         # "timepoint"=substr(sample_prop, 20, 22),
         "class"=substr(sample_prop, 24, 30))
  

ggplot(long_results, aes(x=class, y=freq))+
  geom_point(aes(color=subject_id))+
  geom_line(aes(color=subject_id, group=subject_id))+
  # geom_boxplot(aes(fill=class), outlier.shape = NA)+
  facet_wrap(~cluster_id, scales="free")
