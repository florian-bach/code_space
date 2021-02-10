# Plasma GLMs ####

library(tidyr)
library(dplyr)
library(ggplot2)

setwd("~/PhD/plasma/vac69a/")


data3 <- read.csv("big_plasma_table.csv")


data3[,3:ncol(data3)] <- log10(data3[,3:ncol(data3)])


data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T))

#data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "DoD", "T6", "C45"))


#censor volutneer 03
#data3 <- subset(data3, data3$Volunteer != "v03")


#list_of_dfs_for_glm <- lapply(colnames(data3)[3:ncol(data3)], function(x) data.frame(select(data3, Volunteer, timepoint, x)))
long_data3 <- tidyr::gather(data3, analyte, concentration, colnames(data3)[3:ncol(data3)])
long_data3[,1:3] <- lapply(long_data3[,1:3], as.character)

list_of_dfs_for_glm <- split(long_data3, long_data3$analyte)

list_of_models <- lapply(list_of_dfs_for_glm, function(x) lm(concentration~timepoint+Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))
#list_of_models <- lapply(list_of_dfs_for_glm, function(x) lme4::lmer(concentration~timepoint+(1|Volunteer), data=x))

# 
dod_contrast <- t(diffcyt::createContrast(c(0,0,1,0, rep(0,4))))
t6_contrast <- t(diffcyt::createContrast(c(0,0,0,1, rep(0,4))))
# c45_contrast <- t(diffcyt::createContrast(c(0,1,0,0,0,0,0,0)))
#
 list_of_tests <- lapply(list_of_models, function(x) multcomp::glht(x, t6_contrast))
 list_of_pvalues <- sapply(list_of_tests, function(x) summary(x)$test$pvalues)
#
#
 list_of_adj_pvalues <- p.adjust(list_of_pvalues, method = "fdr")
#
#
list_of_adj_pvalues[order(list_of_adj_pvalues,decreasing = F)]
# 
# write.table(names(list_of_adj_pvalues[order(list_of_adj_pvalues,decreasing = F)]), "~/PhD/plasma/vac69a/analytes_sorted_by_padj_better.txt", sep = "\t", col.names = FALSE, row.names = FALSE)

#list_of_summaries <- lapply(list_of_models, function(x) cbind(summary(x)$coefficients, names(x$data)[3]))
list_of_summaries <- lapply(list_of_models, function(x)cbind(summary(x)$tTable, "analyte"=unique(x$data$analyte)))


df_of_model_results <- data.frame(do.call(rbind, list_of_summaries))
#colnames(df_of_model_results) <- c("Estimate", "SE", "t_value", "raw_p", "Analyte")
df_of_model_results$Coefficient <- rownames(df_of_model_results)


df_of_model_results <- df_of_model_results[!grepl("Intercept", df_of_model_results$Coefficient, fixed=TRUE),]
#df_of_model_results <- df_of_model_results[!grepl("v0", df_of_model_results$Coefficient),]

df_of_model_results$p_adj <- p.adjust(as.numeric(as.character(df_of_model_results$p.value)), method = "fdr")

sig_hits <- subset(df_of_model_results, df_of_model_results$p_adj<0.1)


sig_levels<- as.character(sig_hits[order(sig_hits$p_adj),]$Analyte)



siggy_hits <- sig_hits %>%
  group_by(Analyte) %>%
  top_n(n = -1, wt = p_adj)

# sig_levels <- as.character(siggy_hits[order(siggy_hits$p_adj),]$Analyte)
# # this next step is necessary because TGFbeta is undetectable/unchanged and top_n returns it three times because reasons
# sig_levels <- unique(sig_levels) 
# write.table(sig_levels, "analytes_sorted_by_padj.txt", row.names = FALSE, col.names = "Analyte")

sig_glm_data <- subset(long_data, long_data$Analyte %in% sig_hits$Analyte)
sig_glm_data$AnalyteF <- factor(sig_glm_data$Analyte, levels=sig_levels)


sig_glm_data$p_adj <- sig_hits$p_adj[match(sig_glm_data$Analyte, sig_hits$Analyte)]

sig_glm_plot <- ggplot(sig_glm_data, aes(x=factor(timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer))+
  geom_point()+
  geom_line(aes(group=Volunteer))+
  facet_wrap(~ AnalyteF, scales = "free", ncol=4)+
  scale_y_log10()+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

ggsave(filename = "./figures/glm_sig_analytes_fdr_10-e-1.png", sig_glm_plot, width = 16, height=9)
