diagnosis_plasma <- log10(plasma_data[, c(2,5,8,11,14)])
diagnosis_plasma <- t(data.frame(subset(diagnosis_plasma, rownames(diagnosis_plasma) %in% sig_analytes[1:12])))
#diagnosis_plasma$analyte <- rownames(diagnosis_plasma)

corr_cytof <- cytof_data %>%
  filter(timepoint=="T6") %>%
  select(cluster_id, volunteer, frequency) 



for (analyte in colnames(diagnosis_plasma)){
  
  cbind(corr_cytof, paste(analyte)=)
  
  
}


group_by(cluster_id) %>%
  do(broom::tidy(cor.test(.$frequency, .$distance, method="pearson")))


correlation_dfs <- split(subset(cytof_data, , cytof_data$cluster_id)
names(correlation_dfs) <- lapply(correlation_dfs, function(x)unique(x$cluster_id))
correlation_results <- lapply(correlation_dfs, function(x)cor.test(x$frequency, x$alt, method="pearson"))



