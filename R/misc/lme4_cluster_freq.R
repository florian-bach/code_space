library(ggplot2)
library(nlme)
library(dplyr)
library(tidyr)
library(lme4)

cluster_freqs <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header=T, stringsAsFactors = F)
cluster_freqs$frequency <- cluster_freqs$frequency/100

cluster_freqs$timepoint <- factor(cluster_freqs$timepoint)
cluster_freqs$volunteer <- factor(cluster_freqs$volunteer)
cluster_freqs$sample_id <- factor(cluster_freqs$sample_id)

contrasts(cluster_freqs$timepoint) <- contr.treatment(4, base = 1)

list_of_clusters <- split(cluster_freqs, cluster_freqs$cluster_id)


# test <- data.frame(list_of_clusters[[1]])
# test_model1 <- lme(frequency~timepoint, random = c(~1|volunteer, ~1|sample_id), data=test)
# test_model2 <- lme(frequency~timepoint, random = c(~1|sample_id), data=test)

freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint+(1|timepoint)+(1|volunteer), data=x, family="binomial"))
freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint+volunteer+(1|sample_id), data=x, family="binomial"))
freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint+(1|volunteer)+(1|sample_id), data=x, family = "binomial"))
freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint+volunteer+(1|volunteer)+(1|sample_id), data=x, family = "binomial"))

freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint+(1+timepoint|volunteer), data=x, family = "binomial"))
freq_time_models <- lapply(list_of_clusters, function(x) glmer(frequency~timepoint*volunteer+(1|volunteer), data=x, family = "binomial"))

  # freq_time_models <- lapply(list_of_clusters, function(x) lme(frequency~timepoint, random = c(~1|volunteer, ~1|sample_id), data=x))
# freq_time_models <- lapply(list_of_clusters, function(x) lme(frequency~volunteer+timepoint, random = ~1|sample_id, data=x))

freq_time_results <- lapply(freq_time_models, FUN=function(x)(summary(x)))

pvals <- lapply(freq_time_results, function(x) x$tTable)
pvals <- lapply(pvals, function(x){data.frame(x)})
pvals <- mapply(cbind, pvals, "model"= names(pvals), SIMPLIFY = F)
pvals <- lapply(pvals, function(x) data.frame(x, "timepoint"=rownames(x)))

results_df <- data.table::rbindlist(pvals)
results_df <- subset(results_df, results_df$timepoint!="(Intercept)")
results_df$p_adjust <-  p.adjust(results_df$p.value)

stable_through_time <- subset(results_df, p_adjust>0.05)
varies_through_time <- subset(results_df, p_adjust<0.05)

paste(varies_through_time$model, varies_through_time$timepoint)

# model: frequency~timepoint+(1|sample_id)
# [1] "activated  CD4 CM timepoint4"           "activated  MAIT  timepoint4"           
# [3] "activated PD1+ CD4 EM timepoint4"       "activated PD1+HLADR+ CD4 EM timepoint4"

# model: frequency~timepoint+volunteer+(1|sample_id)
# [1] "activated  CD4 CM timepoint"            "activated  CD8 EM timepoint"           
# [3] "activated  MAIT  timepoint"             "activated  Vd2+  timepoint"            
# [5] "activated PD1+ CD4 EM timepoint"        "activated PD1+HLADR+ CD4 EM timepoint" 
# [7] "activated PD1+HLADR+ CD8 EM timepoint"  "CLA+ CD4 EM timepoint"                 
# [9] "resting  MAIT  timepoint"  

# model: frequency~timepoint+ (1|volunteer)+(1|sample_id)
# [1] "activated  CD4 CM timepointT6"           "activated  CD8 EM timepointT6"          
# [3] "activated  MAIT  timepointT6"            "activated  Vd2+  timepointT6"           
# [5] "activated PD1+ CD4 EM timepointT6"       "activated PD1+HLADR+ CD4 EM timepointT6"
# [7] "activated PD1+HLADR+ CD8 EM timepointT6" "resting  MAIT  timepointDoD"            
# [9] "resting  MAIT  timepointT6"             





freq_time_anova <- lapply(freq_time_models, FUN=function(x)anova(x))
anova_results <- lapply(freq_time_anova, function(x){data.frame(x)})
anova_pvals <- mapply(cbind, anova_results, "model"= names(pvals), SIMPLIFY = F)
anova_pvals <- lapply(anova_pvals, function(x) data.frame(x, "timepoint"=rownames(x)))

anova_df <- data.table::rbindlist(anova_pvals)
anova_df <- subset(anova_df, anova_df$timepoint!="(Intercept)")
anova_df$p_adjust <-  p.adjust(anova_df$p.value)

stable_through_time <- subset(anova_df, p_adjust>0.05)
varies_through_time <- subset(anova_df, p_adjust<0.05)

paste(varies_through_time[order(varies_through
                                _time$timepoint),]$model, varies_through_time[order(varies_through_time$timepoint),]$timepoint)




ggplot(cluster_freqs, aes(x=timepoint, y=frequency))+
  geom_boxplot(aes(fill=timepoint))+
  #geom_point(aes(shape=volunteer))+
  facet_wrap(~cluster_id, scales="free")

ggplot(cluster_freqs, aes(x=timepoint, y=count))+
  geom_boxplot(aes(fill=timepoint))+
  #geom_point(aes(shape=volunteer))+
  facet_wrap(~cluster_id, scales="free")