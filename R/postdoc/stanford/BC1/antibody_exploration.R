# preamble & read in data ###

library(dplyr)
library(patchwork)

#counts NAs per column
count_na <- function(x){table(is.na(x))}

bc1<- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

# make a dictionary of column names, their labels and codings
data_labs <- lapply(bc1, function(x) attributes(x)$labels)
listy_data_labs <- lapply(names(data_labs), function(x)ifelse(is.null(data_labs[[x]]), "n/a", data_labs[x]))
names(listy_data_labs) <- names(data_labs)

#downsize df
ab_columns <- grep("log", colnames(bc1), value = TRUE)
demo_columns <- c("id",
                  "wealthcat",
                  "totalmalariapreg",
                  "anymalariapreg",
                  "gender",
                  "gestage",
                  "anyparasitemia",
                  "timepoint",
                  "MomFinalRx",
                  "nonmalariafebrile",
                  "ChildFinalRx")

slim_bc <- bc1 %>%
  select(all_of(c(demo_columns, ab_columns)))


# NAs are present. count them!
na_tables <- lapply(slim_bc, count_na)
na_counts <- do.call(rbind, na_tables)
# the output for when its all TRUE or all FALSE misses the count for the other boolean, so let's add that for symmetry
na_counts[,2] <- ifelse(na_counts[,1]==nrow(slim_bc), 0, na_counts[,2])

# NAs mostly come from logGST(146) and logpd(448); get rid of those columns, NA.omit and let's roll
slimmer_bc <- slim_bc %>%
  select(-logpd, logGST) %>%
  na.omit()

big_pca <-  prcomp(slimmer_bc[,grep("log", colnames(slimmer_bc), value = TRUE)], center = T)
pca_plot_data <- as.data.frame(cbind(slimmer_bc[, 1:length(demo_columns)], big_pca$x))

# pca loadings figure

loadings_df <- data.frame(big_pca$rotation)
loadings_df$antibody <- rownames(loadings_df)
loadings_df$antibody <- factor(loadings_df$antibody, levels = loadings_df$antibody[order(loadings_df$PC1)])

pc1_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc1_cols) <- loadings_df$antibody[order(loadings_df$PC1)]
                                        
PC1_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC1)]), y=PC1, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


pc2_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc2_cols) <- loadings_df$antibody[order(loadings_df$PC2)]

loadings_df$antibody <- factor(loadings_df$antibody, levels = loadings_df$antibody[order(loadings_df$PC2)])
PC2_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC2)]), y=PC2, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc2_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


combo_plot <- PC1_plot+PC2_plot



# tl;dr this is the only one where there's seperation
ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(timepoint)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

# other PCAs ####
wealthcat_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(wealthcat)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

  
totalmalariapreg_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(totalmalariapreg)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
anymalariapreg_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(anymalariapreg)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
gender_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(gender)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
  
gestage_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=gestage))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    scale_color_gradient(low="#073376", high="#FFAE42")+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
  
anyparasitemia_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=anyparasitemia))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    scale_color_gradient(low="#073376", high="#FFAE42")+
    ggtitle("Antibody Titre PCA")
  
  
timepoint_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(timepoint)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")

nonmalariafebrile_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=nonmalariafebrile))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

MomFinalRx_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(MomFinalRx)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

ChildFinalRx_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(ChildFinalRx)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")



#long format baby
ginormous_df <- slimmer_bc %>%
  mutate(timepoint=ifelse(timepoint==1, "1st", ifelse(timepoint==2, "2nd", "3rd"))) %>%
  pivot_longer(cols = (length(demo_columns)+1):ncol(slimmer_bc), names_to = "antibody", values_to = "log_conc") %>%
  group_by(antibody)

# data vis ####  

# LOOK at the data

all_ab_plot <- ggplot(ginormous_df, aes(x=timepoint, y=10^(log_conc)))+
  facet_wrap(~antibody, scales = "free")+
  geom_point(aes(color=antibody))+
  geom_boxplot(aes(fill=antibody))+
  scale_y_log10()+
  theme_minimal()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  xlab("")+
  ylab("Concentration")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/all_abs.png", all_ab_plot, height = 10, width=8, bg = "white")


all_ab_plot2 <- ggplot(ginormous_df, aes(x=timepoint, y=10^(log_conc)))+
  facet_wrap(~antibody, scales = "free")+
  geom_violin(aes(fill=antibody), draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_y_log10()+
  theme_minimal()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  xlab("")+
  ylab("Concentration")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/all_ab_plot2.png", all_ab_plot2, height = 10, width=8, bg = "white")




# longitudinal models ####

# split big dataframe into 22 dataframes, one for each antibody
list_of_dfs <- split(ginormous_df, ginormous_df$antibody)

# run a model for each with id as random effect
list_of_models <- lapply(list_of_dfs, function(x) lmer(log_conc~timepoint+(1|id), data=x))

# define conrtrasts
sec_contrast <- t(matrix(c(0,1,0)))
ter_contrast <- t(matrix(c(0,0,1)))
sec_ter_contrast <- t(matrix(c(0,-1,1)))

list_of_tests <- lapply(list_of_models, function(x) multcomp::glht(x, sec_ter_contrast))
list_of_pvalues <- sapply(list_of_tests, function(x) summary(x)$test$pvalues)
#
#
list_of_adj_pvalues <- sort(p.adjust(list_of_pvalues, method = "fdr"))

table(list_of_adj_pvalues<0.05) #all sic at 2nd and 3rd relative to first; 12/22 from 2nd to 3rd
sig_2_3_abs <- subset(list_of_adj_pvalues, list_of_adj_pvalues<0.05)

# make results table
list_of_coefficients <- lapply(list_of_models, function(x)coef(x)$id[1,])
coefficient_df <- do.call(rbind,list_of_coefficients)
coefficient_df$ab <- rownames(coefficient_df)
coefficient_df$raw_p <- list_of_pvalues[coefficient_df$ab]
coefficient_df$p_adj <- list_of_adj_pvalues[coefficient_df$ab]


# vis modelling results
sig_2_3_ab <- ginormous_df %>%
  filter(antibody %in% names(sig_2_3_abs)) %>% 
  ggplot(., aes(x=timepoint, y=10^(log_conc)))+
  facet_wrap(~antibody, scales = "free")+
  geom_violin(aes(fill=antibody), draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_y_log10()+
  theme_minimal()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  xlab("")+
  ylab("Concentration")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/sig_2_3_ab.png", sig_2_3_ab, height = 10, width=8, bg = "white")


all_ab_plot3 <- ginormous_df %>%
  mutate(sig_2_3=ifelse(antibody %in% names(sig_2_3_abs), "sig", "not_sig"))%>%
  ggplot(., aes(x=timepoint, y=10^(log_conc)))+
  facet_wrap(~antibody, scales = "free")+
  geom_violin(aes(fill=antibody, alpha=sig_2_3), draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_y_log10()+
  theme_minimal()+
  scale_alpha_manual(values = c(0.1, 0.9))+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  xlab("")+
  ylab("Concentration")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/all_ab_plot3.png", all_ab_plot3, height = 10, width=8, bg = "white")


# plans ####
# select rows where incident malaria episode occurred within ~30 days?
# select rows where anyparasitaemia occurred within ~30 days


