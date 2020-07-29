list_of_models1 <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]*x[,3]+x[,3]+x[,1], data=x))
# this model showed that IL1RA, MPO were higher in second, TF higher in third
list_of_models2 <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]+x[,3]+x[,1], data=x))


plot(list_of_models1[[7]])

plot(list_of_models2[[7]])


# then we make a list where we extract the glm results as a matrix and add a column with the name of the analyte
list_of_summaries1 <- lapply(list_of_models1, function(x) cbind(summary(x)$coefficients, names(x$data)[4]))
list_of_summaries2 <- lapply(list_of_models2, function(x) cbind(summary(x)$coefficients, names(x$data)[4]))

# now we pull the list of matrices together as one big data frame
df_of_model_results <- data.frame(do.call(rbind, list_of_summaries1))
#rename the df and add a column for the coefficient name (atm these are just rownames)
colnames(df_of_model_results) <- c("Estimate", "SE", "t_value", "raw_p", "Analyte")
df_of_model_results$Coefficient <- rownames(df_of_model_results)

#here we get rid of the coefficients for the intercept value (we don't care about those)
df_of_model_results <- df_of_model_results[!grepl("Intercept", df_of_model_results$Coefficient),]
# here we get rid of the coefficients that are specific to the volunteer fixed effect; a low p value here indicates that
# a volunteer has significantly more or less of analyte X at Baseline- we kind of care about this, but not for this part of the analysis
# I think it's relevant to exlude those things because each coefficient gets a p value and for adjusting for multiple testing
# having irrelevant tests in there adds noise to the results
df_of_model_results <- df_of_model_results[!grepl("V", df_of_model_results$Coefficient),]


interactive_terms <- df_of_model_results[grepl("...", df_of_model_results$Coefficient),]
df_of_second<- interactive_terms[grepl("Second", interactive_terms$Coefficient, fixed=T),]
df_of_third <- interactive_terms[grepl("Third", interactive_terms$Coefficient, fixed=T),]

combo_df <- rbind(df_of_second, df_of_third)


#interactive_terms <- interactive_terms[grepl("DoD", interactive_terms$Coefficient),]


combo_df$p_adj <- p.adjust(as.numeric(as.character(combo_df$raw_p)), method = "fdr")


sig_hits <- subset(combo_df, combo_df$p_adj<=0.1)













####################################



list_of_models1 <- lapply(list_of_dfs_for_glm, function(x) lme4::glmer(x[,4]~timepoint_with_DoD*N_infection+(1|Volunteer_code), data=x))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes

list_of_models2 <- lapply(list_of_dfs_for_glm, function(x) lme4::glmer(x[,4]~timepoint_with_DoD+N_infection+(1|Volunteer_code), data=x, family = ))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes

list_of_models1 <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(x[,4]~x[,2]*x[,3]+x[,1], random=1~Volunteer_code, data=x))
list_of_models2 <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(x[,4]~x[,2]*x[,3]+x[,1], random=1~Volunteer_code, data=x))

# list_of_modelsX <- vector(mode = "list", length = 39)
# this approach shows no change between infections
list_of_models_times <- list()

list_of_models_plus <- list()
list_of_lm_models_plus <- list()

for(i in 1:length(list_of_dfs_for_glm)){
  dat <- list_of_dfs_for_glm[[i]]
  print(unique(colnames(dat)[4]))
  biomarker <- colnames(dat)[4]
  colnames(dat)[4] <- "Analyte"
  lm1 <- nlme::lme(Analyte~timepoint_with_DoD+N_infection, random= ~1 | Volunteer_code, data=dat)
  #lm1 <- lme4::lmer(Analyte~timepoint_with_DoD*N_infection+(1|Volunteer_code), data=dat)
  
  # lm2 <- lm(Analyte~timepoint_with_DoD+N_infection+Volunteer_code, data=dat)
  list_of_models_plus[[paste(biomarker)]] <- lm1
  # list_of_lm_models_plus[[paste(biomarker)]] <- lm2
}

AIC_diff <- cbind(sapply(list_of_models_plus, AIC), sapply(list_of_lm_models_plus, AIC))
AIC_diff <- subset(AIC_diff, AIC_diff> 3)


#the purely additive model is the best for most, except
# CXCL10    CXCL9     IFNγ     IL10    IL1RA     IL21      MPO     PAI1 
#having the additional additive term in the interaction model has 0 change on AIC
  
resx <- lapply(list_of_models_plus, function(x) anova(x))

ress <- do.call(rbind, resx)
#ress <- ress[!grepl("(Intercept)", rownames(ress)),]


#ress <- ress[grepl(":", rownames(ress), fixed = T),]

ress <- ress[!grepl("(Intercept)", rownames(ress)) & !grepl("N_infection", rownames(ress), fixed=TRUE),]

#ress <- ress[!grepl("N_infection", rownames(ress), fixed = T),]

ress$p_adj <- p.adjust(as.numeric(as.character(ress$`p-value`)), method = "fdr")
ress <- ress[ress$p_adj<=0.05,]

# the interaction model returns 23 analytes as significant through time, the additive one 21 (it misses IL33 and leptin, otherwise it's the same)
# however, generally p values are lower for the additive model which might explain the difference
# the interaction term is significant for 8 analytes, ICAM1, IL10, IL1RA, IL21, MPO, TF, TIMP1, TNFRII

# in general the purely additive model is better
#here are the analytes for which this is not true

subsetter <- AIC_diff<0
res <- subset(list_of_dfs_for_glm, subsetter)
sapply(res, function(x)names(x)[4])






data3 <- spread(long_data, Analyte, Concentration)
data4 <- data3[,rownames(AIC_diff)[-8]]



big_pca <-  prcomp(data4, center = T)
big_pca2 <- cbind(data3[, c(1,2,4,5)], big_pca$x)



all_vols_together_vol_color <- ggplot(big_pca2, aes(x=PC1, y=PC2))+
  geom_text(aes(label=Sample_ID, color=N_infection), size=5)+
  #scale_color_manual(values=my_paired_palette)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,3]*100, "%", sep = ""))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA),
        axis.title = element_text(size=12),
        #axis.text = element_blank(),
        
        plot.title = element_text(size=14, hjust=0.5)
  )


############################################################

list_of_models_plus <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]+x[,3]+x[,1], data=x))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes
names(list_of_models_plus) <- sapply(list_of_models_plus, function(x) unique(names(x$data)[4]))

#                                                         Analyte~Timepoint*N_infection+Volunteer
list_of_models_times <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]*x[,3]+x[,3]+x[,1], data=x))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes
names(list_of_models_times) <- sapply(list_of_models_times, function(x) unique(names(x$data)[4]))

list_of_comparisons <- lapply(names(list_of_models_plus), function(x)anova(list_of_models_times[[x]], list_of_models_plus[[x]], test="Chisq"))
names(list_of_comparisons) <- sapply(list_of_models_times, function(x) unique(names(x$data)[4]))

list_of_res <- vector(length = 39)

list_of_res <- sapply(list_of_comparisons, function(x) x$`Pr(>Chi)`[2])
list_of_res <- p.adjust(list_of_res, method = "fdr")
list_of_sig <- subset(list_of_res, list_of_res<=0.05)

names(list_of_res) <- sapply(list_of_models_times, function(x) unique(names(x$data)[4]))


trimmed_list <- lapply()




