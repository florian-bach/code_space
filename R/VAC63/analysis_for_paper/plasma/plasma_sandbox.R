list_of_models <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]*x[,3]+x[,1], data=x))
# this model showed that IL1RA, MPO were higher in second, TF higher in third



# then we make a list where we extract the glm results as a matrix and add a column with the name of the analyte
list_of_summaries <- lapply(list_of_models, function(x) cbind(summary(x)$coefficients, names(x$data)[4]))

# now we pull the list of matrices together as one big data frame
df_of_model_results <- data.frame(do.call(rbind, list_of_summaries))
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


combo_df$p_adj <- p.adjust(as.numeric(as.character(combo_df$raw_p)), method = "BH")


sig_hits <- subset(combo_df, combo_df$p_adj<=0.05)













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

list_of_modelsX <- list()
for(i in 1:length(list_of_dfs_for_glm)){
  dat <- list_of_dfs_for_glm[[i]]
  print(unique(colnames(dat)[4]))
  biomarker <- colnames(dat)[4]
  colnames(dat)[4] <- "Analyte"
  lm <- nlme::lme(Analyte~timepoint_with_DoD+N_infection, random = ~ 1 | Volunteer_code, data=dat)
  list_of_modelsX[[paste(biomarker)]] <- lm

}
  
resx <- lapply(list_of_modelsX, function(x) anova(x))

ress <- do.call(rbind, resx)
ress <- ress[!grepl("(Intercept)", rownames(ress)),]



ress <- ress[grepl(":", rownames(ress), fixed = T),]


ress <- ress[!grepl("N_infection", rownames(ress), fixed = T),]

ress$p_adj <- p.adjust(as.numeric(as.character(ress$`p-value`)), method = "fdr")
ress <- ress[ress$p_adj<=0.05,]


AIC_diff <- sapply(list_of_models1, AIC)-sapply(list_of_models2, AIC)

# in general the purely additive model is better
#here are the analytes for which this is not true

subsetter <- AIC_diff<0
res <- subset(list_of_dfs_for_glm, subsetter)
sapply(res, function(x)names(x)[4])
