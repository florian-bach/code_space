library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)

cluster_freqs <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_counts_and_freqs.csv", header=T, stringsAsFactors = F)

modelz <- lmer(frequency~timepoint+volunteer+(1|sample_id), data = cluster_freqs)

modelz <- lm(frequency~timepoint, data = cluster_freqs)

list_of_models <- data.frame()

list_of_models <- cluster_freqs %>%
  group_by(cluster_id) %>%
  do(model = lmer(frequency~timepoint+(1|volunteer), data=.)) %>%
  do(rbind(list_of_models, data.frame("model_id"=unique(.$cluster_id), "coefficients"=rownames(summary(model)$coefficients),  summary(model)$coefficients)))

list_of_models$P_adjust <- p.adjust(list_of_models$Pr...t.., "BH")

vary_through_time <- subset(list_of_models, list_of_models$P_adjust<0.05)

stable_through_time <- subset(list_of_models, list_of_models$P_adjust>0.05)



