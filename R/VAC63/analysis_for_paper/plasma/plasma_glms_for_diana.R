# making   the glms ####

# the input data should be a data frame with a non-numeric column for each category of metadata (volunteer, timepoint, n_infection)
# ideally, the metadata columns should be either the first or last columns of the df, to make your life easier.. e.g.

# glimpse(data3)
# Rows: 20
# Columns: 41
# $ Volunteer       <chr> "v002", "v002", "v002", "v002", "v003", "v003", "v003", "v003", "v005", "v005", "v005", "v005", "v0…
# $ timepoint       <chr> "C-1", "C+45", "DoD", "T+6", "C-1", "C+45", "DoD", "T+6", "C-1", "C+45", "DoD", "T+6", "C-1", "C+45…
# $ Ang2            <dbl> 2.959485, 3.100371, 2.946339, 3.178689, 2.988001, 2.831268, 3.040207, 3.038620, 3.150142, 3.152900,…
# $ Arginase        <dbl> 4.088420, 4.369939, 3.980276, 3.742254, 4.263162, 4.241198, 4.105749, 4.026615, 3.863442, 3.704579,…
# $ CCL2            <dbl> 2.471189, 2.413769, 2.902661, 2.406148, 2.300813, 2.447917, 2.712178, 2.216377, 2.142921, 2.373611,…




library(dplyr)
library(tidyr)
library(ggplot2)


#library(rgl) # want some 3d?

#data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")

plasma <- read.csv("~/PhD/plasma/vac63/Paired_data_log2transf_proteome.csv", header=T, stringsAsFactors = F)
plasma <- plasma[,-1]
plasma <- plasma %>%
  mutate(timepoint_with_DoD=gsub("C-1", "Baseline", timepoint_with_DoD)) %>%
  mutate(timepoint_with_DoD=gsub("+", "", timepoint_with_DoD, fixed = T)) %>%
  mutate(timepoint_with_DoD=gsub("Diagnosis", "DoD", timepoint_with_DoD)) %>%
  mutate(timepoint_with_DoD=gsub("D6", "T6", timepoint_with_DoD))

plasma$timepoint_with_DoD <- factor(plasma$timepoint_with_DoD, levels=c("Baseline", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C14", "DoD","T6","C28", "C45"))


#subset plasma dataset so that every volunteer has the same set of timepoints ####

# make column with unique identifier for timecourses
plasma <- cbind(list("Sample_ID"=paste(plasma$Volunteer_code, plasma$N_infection, sep="_")), plasma)

split_plasma <- split(plasma, plasma$Sample_ID)

#split dataframe by timecourse and list all timepoints within each timecourse
timepoint_lists <- lapply(split_plasma, function(x)sort(unique(x$timepoint_with_DoD)))

#count occurence of each timepoint and subset data by timepoints that occur as often as there are timecourse
timepoint_counts <- table(unlist(truth_table))
omnipresent_timepoints <- names(subset(timepoint_counts, timepoint_counts==length(unique(plasma$Sample_ID))))
#"Baseline" "C6"       "C8"       "DoD"


#subset plasma dataset so that every volunteer has the same set of timepoints
plasma <- filter(plasma, timepoint_with_DoD %in% c("Baseline", "C6", "C8", "DoD"))


long_data <- gather(plasma, Analyte, Concentration, colnames(plasma)[6:ncol(plasma)])

#long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)

long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)


# exploratory plots ####
# 
ggplot(long_data, aes(x=timepoint_with_DoD, y = Concentration))+
  #geom_point(aes(shape=N_infection))+
  geom_boxplot(aes(fill=N_infection))+
  facet_wrap(~Analyte, scales="free")+
  theme_minimal()


####     pca stuff ####


data3 <- spread(long_data, Analyte, Concentration)
#data3[,3:ncol(data3)] <- (data3[,3:ncol(data3)])



# here's where the metadata columns being at the start comes in handily: we use the column names starting  with the analyte
# columns to make a list of data frame each containing the metadata and one column with one analyte
list_of_dfs_for_glm <- lapply(colnames(data3)[6:ncol(data3)], function(x) data.frame(select(data3, Volunteer_code, timepoint_with_DoD, N_infection, all_of(x))))

# now we use lapply to make a glm for each data frame using the formula glm(analyte~timepoint+volunteer). I was lazy and used
# indexing to do this rather than to name things properly and call them explicitly
# for your data you should test a range of formulas to test the effect of n_infection, but try this very basic one first which
# will pick out the conserved changes that occur in each infection



# colnames(list_of_dfs_for_glm[[1]])
# [1] "Volunteer_code"     "timepoint_with_DoD" "N_infection"        "Ang2"


#                                                         Analyte~Timepoint+Volunteer+N_infection
list_of_models <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]+x[,3]+x[,1], data=x))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes


#                                                         Analyte~Timepoint*N_infection+Volunteer
list_of_models <- lapply(list_of_dfs_for_glm, function(x) glm(x[,4]~x[,2]*x[,3]+x[,3]+x[,1], data=x))
# this model showed that CXCL9, IL18, IL21 were generally higher through time during second infectoin, IFNa was generally lower in third
# this refers to their intercepts, not their slopes




# then we make a list where we extract the glm results as a matrix and add a column with the name of the analyte
list_of_summaries <- lapply(list_of_models_times, function(x) cbind(summary(x)$coefficients, names(x$data)[4]))

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

#evaluate significance of N_infection coefficients
df_of_second<- df_of_model_results[grepl("Second", df_of_model_results$Coefficient, fixed=T),]
df_of_third <- df_of_model_results[grepl("Third", df_of_model_results$Coefficient, fixed=T),]

df_of_n_infection <- rbind(df_of_second, df_of_third)
df_of_n_infection$p_adj <- p.adjust(as.numeric(as.character(df_of_n_infection$raw_p)), method = "BH")

n_infection_sig <- subset(df_of_n_infection, df_of_n_infection$p_adj<=0.05)


# examine significance of timepoint coefficients
df_of_model_results <- df_of_model_results[!grepl("Second", df_of_model_results$Coefficient, fixed=T),]
df_of_model_results <- df_of_model_results[!grepl("Third", df_of_model_results$Coefficient, fixed=T),]



#adjust p value with Bonferroni Hochberg
df_of_model_results$p_adj <- p.adjust(as.numeric(as.character(df_of_model_results$raw_p)), method = "BH")

#get significant hits
sig_hits <- subset(df_of_model_results, df_of_model_results$p_adj<=0.05)

#have a look
sig_hits

### make the figure of significant results ####


#get the order according to p_value
sig_levels <- as.character(sig_hits[order(sig_hits$p_adj),]$Analyte)

# subset big long data frame to only contain significant analytes and make a column with reordered factors, this allows
# appropriate ordering of the facets with facet_wrap() later on
sig_glm_data <- subset(long_data, long_data$Analyte %in% sig_hits$Analyte)
sig_glm_data$AnalyteF <- factor(sig_glm_data$Analyte, levels=unique(sig_levels))


sig_glm_plot <- ggplot(sig_glm_data, aes(x=timepoint_with_DoD, y=Concentration))+
  geom_point(aes(shape=N_infection), position = position_jitter(w = 0.17, h = 0))+
  #geom_point(aes(shape=N_infection))+
  geom_boxplot(aes(fill=timepoint_with_DoD))+
  facet_wrap(~ AnalyteF, scales = "free", ncol=6)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values=list("First"="Red", "Second"="Green", "Third"="Blue"))+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))
  


ggsave(filename = "~/PhD/plasma/vac63/analyte~timepoint+volunter*n_infection.png", sig_glm_plot, width = 16, height=9)






sig_glm_plot <- ggplot(sig_glm_data, aes(x=timepoint_with_DoD, y=Concentration))+
  #geom_point(aes(shape=N_infection), position = position_jitter(w = 0.17, h = 0))+
  geom_point(aes(shape=N_infection))+
  geom_boxplot(aes(fill=N_infection))+
  facet_wrap(~ AnalyteF, scales = "free", ncol=6)+
  #scale_y_log10()+
  theme_bw()+
  #scale_color_manual(values=list("First"="Red", "Second"="Green", "Third"="Blue"))+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))



ggsave(filename = "~/PhD/plasma/vac63/analyte~timepoint+volunter*n_infection_by_n_infection.png", sig_glm_plot, width = 16, height=9)







