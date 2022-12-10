# preamble & read in data ###

#palettes
time_palette <- colorspace::sequential_hcl(n=4, "RdPu")[1:3]
pc1_cols <- colorspace::sequential_hcl(23, palette = "Purple Yellow")

#libraries 
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(purrr)
library(lme4)
library(knitr)

# functions

#counts NAs per column
count_na <- function(x){table(is.na(x))}
# mode finder
mode_finder <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

modelable_antigens <- c("Tet Tox", "SBP1", "Rh5", "PfSEA", "PfAMA1", "Hyp2", "HSP40 Ag1", "GST", "GEXP", "CSP GENOVA")


#data
bc1 <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

bc1$age <- bc1$date -bc1$dob

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  select(id)

# make a dictionary of column names, their labels and codings
data_labs <- lapply(bc1, function(x) attributes(x)$labels)
listy_data_labs <- lapply(names(data_labs), function(x)ifelse(is.null(data_labs[[x]]), "n/a", data_labs[x]))
names(listy_data_labs) <- names(data_labs)

#downsize df
ab_columns <- grep("log", colnames(bc1), value = TRUE)

demo_columns <- c("id",
                  "wealthcat",
                  "age",
                  "totalmalariapreg",
                  "anymalariapreg",
                  "gender",
                  "gestage",
                  "anyparasitemia",
                  "timepoint",
                  "MomFinalRx",
                  "nonmalariafebrile",
                  "ChildFinalRx",
                  "malaria",
                  "mal0to6",
                  "febrile0to6",
                  "incidentmalaria",
                  "febrile6to12",
                  "mal6to12")

slim_bc <- bc1 %>%
  dplyr::select(c(demo_columns, ab_columns))


# NAs are present. count them!
na_tables <- lapply(slim_bc, count_na)
na_counts <- do.call(rbind, na_tables)
# the output for when its all TRUE or all FALSE misses the count for the other boolean, so let's add that for symmetry
na_counts[,2] <- ifelse(na_counts[,1]==nrow(slim_bc), 0, na_counts[,2])

# NAs mostly come from logGST(146) and logpd(448); get rid of those columns, NA.omit and let's roll
slimmer_bc <- slim_bc %>%
  dplyr::select(-logpd, -logGST) %>%
  na.omit()


# all data PCA ####

big_pca <-  prcomp(slimmer_bc[,grep("log", colnames(slimmer_bc), value = TRUE)], center = T)
pca_plot_data <- as.data.frame(cbind(slimmer_bc, big_pca$x))


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
age_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=as.numeric(age)/365))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")


wealthcat_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(wealthcat)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

  
totalmalariapreg_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=totalmalariapreg))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
  viridis::scale_color_viridis(option="B")+
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

incident_malaria_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=incidentmalaria))+
  geom_point()+
  facet_wrap(~timepoint)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B", direction = -1)+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")



#

# redo pca seperately for each timepoint ####



list_of_dfs <- split(slimmer_bc, slimmer_bc$timepoint)

list_of_pca <-  lapply(list_of_dfs, function(x) prcomp(x[,grep("log", colnames(x), value = TRUE)], center = T))

list_of_plot_data <- lapply(1:3, function(y) as.data.frame(cbind(list_of_dfs[[y]], list_of_pca[[y]]$x)))


t1 <- ggplot(list_of_plot_data[[1]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[1]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[1]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T1")

t2 <- ggplot(list_of_plot_data[[2]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[2]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[2]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T2")

t3 <- ggplot(list_of_plot_data[[3]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[3]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[3]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T3")

t1+t2+t3



list_of_loadings_df <- lapply(list_of_pca, function(x) data.frame(x$rotation))

list_of_pc123 <- lapply(1:3, function(x) cbind(rownames(list_of_loadings_df[[x]]),
                                               list_of_loadings_df[[x]]$PC1,
                                               list_of_loadings_df[[x]]$PC2,
                                               list_of_loadings_df[[x]]$PC3,
                                               rep(names(list_of_pca[x]), length(list_of_loadings_df[[x]]))))

pc123 <- as.matrix(do.call(rbind, list_of_pc123))
colnames(pc123) <- c("antibody", "PC1", "PC2", "PC3", "timepoint")

class(pc123[,2:5]) <- "numeric"
pc123 <- as.data.frame(pc123)


pc123$PC1 <- as.numeric(pc123$PC1)
pc123$PC2 <- as.numeric(pc123$PC2)
pc123$PC3 <- as.numeric(pc123$PC3)

ggplot(pc123, aes(x=tidytext::reorder_within(antibody, PC1, timepoint, mean), y=PC1, fill=antibody))+
  geom_bar(stat = "identity")+
  facet_wrap(~timepoint, scales = "free_x")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  tidytext::scale_x_reordered()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none",
        axis.text.y = element_blank())



# data vis ####  

# LOOK at the data

# how many kids; 189
# n_distinct(bc1$id)

# how many timepoint sof each
# 1   2   3 
# 178 183 179 
# table(bc1$timepoint)

# how many full timecourses
# 1   2   3 
# 6 15 168 
bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  count(n_time)

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  select(id)

# long format baby
ginormous_df <- slimmer_bc %>%
  mutate(timepoint=ifelse(timepoint==1, "1st", ifelse(timepoint==2, "2nd", "3rd"))) %>%
  pivot_longer(cols = (length(demo_columns)+1):ncol(slimmer_bc), names_to = "antibody", values_to = "log_conc") %>%
  group_by(antibody)

# age at timepoint

ggplot(slimmer_bc, aes(x=timepoint, y=age, group=factor(id), color=factor(id)))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  theme(legend.position = "none")



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

ginormous_df <- slimmer_bc %>%
  mutate(timepoint=ifelse(timepoint==1, "1st", ifelse(timepoint==2, "2nd", "3rd"))) %>%
  pivot_longer(cols = (length(demo_columns)+1):ncol(slimmer_bc), names_to = "antibody", values_to = "log_conc") %>%
  group_by(antibody)

# define conrtrasts
sec_contrast <- t(matrix(c(0,1,0)))
ter_contrast <- t(matrix(c(0,0,1)))
sec_ter_contrast <- t(matrix(c(0,-1,1)))

# brave new world purrrlicious way of doing it:
# make grouped df, nest into one, add column for model & summary, extract values from each to create even more new columns for p values and p_adj values
purrrf <- ginormous_df %>%
  group_by(antibody) %>%
  nest() %>%
  mutate(model=map(data, ~lmer(log_conc~timepoint+(1|id), data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(t2_t1=map(model, ~multcomp::glht(., sec_contrast)),
         t2_t1_p=map_dbl(t2_t1, ~summary(.)$test$pvalues),
         t2_t1_p_adj=map_dbl(t2_t1_p, ~p.adjust(.))) %>%
  mutate(t3_t1=map(model, ~multcomp::glht(., ter_contrast)),
         t3_t1_p=map_dbl(t3_t1, ~summary(.)$test$pvalues),
         t3_t1_p_adj=map_dbl(t3_t1_p, ~p.adjust(.))) %>%
  mutate(t3_t2=map(model, ~multcomp::glht(., sec_ter_contrast)),
         t3_t2_p=map_dbl(t3_t2, ~summary(.)$test$pvalues),
         t3_t2_p_adj=map_dbl(t3_t2_p, ~p.adjust(.)))

#make results table that includes the coef & p_adj values for the sec_ter_contrast
results_table <- purrrf %>%
  mutate(coef=map_dbl(model, ~-coef(.)$id[1,2]+coef(.)$id[1,3]))%>%
  select(antibody, coef, t3_t2_p_adj) %>%
  ungroup()

# more tidy content; doesn't have p values though
# alternative_results_table <- purrrf %>%
#   mutate(tidy = map(model, broom.mixed::tidy)) %>%
#   select(antibody, tidy)%>%
#   unnest(tidy)

sig_2_3_abs <- results_table %>%
  filter(t3_t2_p_adj<0.05) %>%
  select(antibody)


# vis modelling results
sig_2_3_ab <- ginormous_df %>%
  filter(antibody %in% sig_2_3_abs$antibody) %>% 
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


# redo pca seperately for each timepoint ####



list_of_dfs <- split(slimmer_bc, slimmer_bc$timepoint)

list_of_pca <-  lapply(list_of_dfs, function(x) prcomp(x[,grep("log", colnames(x), value = TRUE)], center = T))

list_of_plot_data <- lapply(1:3, function(y) as.data.frame(cbind(list_of_dfs[[y]], list_of_pca[[y]]$x)))


t1 <- ggplot(list_of_plot_data[[1]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[1]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[1]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T1")

t2 <- ggplot(list_of_plot_data[[2]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[2]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[2]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T2")

t3 <- ggplot(list_of_plot_data[[3]], aes(x=PC1, y=PC2, color=logMSP1))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(list_of_pca[[3]])[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(list_of_pca[[3]])[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(option="B")+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("T3")

t1+t2+t3



list_of_loadings_df <- lapply(list_of_pca, function(x) data.frame(x$rotation))

list_of_pc123 <- lapply(1:3, function(x) cbind(rownames(list_of_loadings_df[[x]]),
                                               list_of_loadings_df[[x]]$PC1,
                                               list_of_loadings_df[[x]]$PC2,
                                               list_of_loadings_df[[x]]$PC3,
                                               rep(names(list_of_pca[x]), length(list_of_loadings_df[[x]]))))

pc123 <- as.matrix(do.call(rbind, list_of_pc123))
colnames(pc123) <- c("antibody", "PC1", "PC2", "PC3", "timepoint")

class(pc123[,2:5]) <- "numeric"
pc123 <- as.data.frame(pc123)


pc123$PC1 <- as.numeric(pc123$PC1)
pc123$PC2 <- as.numeric(pc123$PC2)
pc123$PC3 <- as.numeric(pc123$PC3)

ggplot(pc123, aes(x=tidytext::reorder_within(antibody, PC1, timepoint, mean), y=PC1, fill=antibody))+
  geom_bar(stat = "identity")+
  facet_wrap(~timepoint, scales = "free_x")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  tidytext::scale_x_reordered()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none",
        axis.text.y = element_blank())


#predicting MSP1 Ab levels ####
long_demo_df <- ginormous_df %>%
  filter(antibody=="logMSP1")

#nah
wealth_mod <- lmer(log_conc~wealthcat+timepoint+(1|id), data=long_demo_df)
gender_mod <- lmer(log_conc~gender+timepoint+(1|id), data=long_demo_df)
totalmalariapreg_mod <- lmer(log_conc~totalmalariapreg+timepoint+(1|id), data=long_demo_df)
anyparasitemia_mod <- lmer(log_conc~anyparasitemia+timepoint+(1|id), data=long_demo_df)
nonmalariafebrile_mod <- lmer(log_conc~nonmalariafebrile+timepoint+(1|id), data=long_demo_df)

#maybe 
gestage_mod <- lmer(log_conc~gestage+timepoint+(1|id), data=long_demo_df)

ggplot(long_demo_df, aes(x=gestage, y=log_conc, color=as.factor(timepoint)))+
  geom_point()+
  scale_color_manual(values=time_palette)+
  theme_minimal()


#  but wait there's a something fishy happening; plot histograms per antibody and timepoint
ab_time_conc_histo <- ggplot(ginormous_df, aes(x=10^log_conc, fill=antibody))+
  facet_grid(antibody~timepoint)+
  geom_histogram(bins = 40)+
  scale_fill_manual(values=pc1_cols)+
  scale_x_continuous(trans="log10", breaks = 10^seq(-12, 3, by=3))+
  theme_minimal()+
  theme(legend.position = "none",
        strip.text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(size=12),
        strip.text.y = element_text(angle = 0)
        )

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/ab_time_conc_histo.png",ab_time_conc_histo, bg="white", width=9, height=16)


# make a table-ish heatmap of the mode for each antibody-timepoint combination
ab_nest <- ginormous_df %>%
  select(antibody, timepoint, log_conc) %>%
  group_by(antibody, timepoint) %>%
  nest() %>%
  group_by(antibody, timepoint)%>%
  summarise(mode= map_dbl(data, ~mode_finder(.$log_conc)),
            min = map_dbl(data, ~min(.$log_conc)),
            max = map_dbl(data, ~max(.$log_conc))) %>%
  mutate(is_weird=ifelse(round(mode, digits = 6)==-2.892790, "yep", "nope"))
#select(antibody, timepoint, mode) #%>%
# pivot_wider(names_from = "timepoint", values_from = "mode", names_glue = "{timepoint}_{.value}")


ab_mode_heatmap <- ggplot(ab_nest, aes(x=timepoint, y=antibody, label=round(mode, digits = 4)))+
  geom_tile(fill="white")+
  geom_text(aes(color=is_weird))+
  scale_fill_viridis_c(option = "B")+
  scale_color_manual(values=c("black", "red"))+
  theme_minimal()+
  scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/ab_mode_heatmap.png", ab_mode_heatmap, width=4, height=8, bg="white")

# recent malaria episodes ####

bc1_all <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs all visit database.dta")

#downsize df
clinical_columns <- c("fever",
                      "alt",
                      "anymalaria",
                      "parsdens",
                      "fever",
                      "febrile",
                      "feverdur",
                      "wbc",
                      "wbcgrade",
                      "neutro",
                      "neutrograde",
                      "plt",
                      "pltgrade",
                      "MSPinfectiondich",
                      "MPinfectiondich")
demo_columns2 <- c("id",
                  "date",
                  "ageyrs")

slim_bc1_all <- bc1_all %>%
  select(all_of(c(demo_columns2, clinical_columns)))


malaria_episodes <- slim_bc1_all %>%
  filter(anymalaria==1) %>%
  select(id, date, ageyrs)

tp1 <- slimmer_bc %>%
  filter(timepoint==1)%>%
  select(id, date)

tp2 <- slimmer_bc %>%
  filter(timepoint==2)%>%
  select(id, date)

tp3 <- slimmer_bc %>%
  filter(timepoint==3)%>%
  select(id, date)

malaria_episodes$time_to_tp1 <- malaria_episodes$date-tp1$date[match(malaria_episodes$id, tp1$id)]
malaria_episodes$time_to_tp2 <- malaria_episodes$date-tp2$date[match(malaria_episodes$id, tp2$id)]
malaria_episodes$time_to_tp3 <- malaria_episodes$date-tp3$date[match(malaria_episodes$id, tp3$id)]

malaria_episodes_near_sample_dates <- malaria_episodes %>%
  filter(between(as.numeric(time_to_tp3), -60, 0))

# garbage / sandbox


# base r making of models, testing and results tables
# split big dataframe into 22 dataframes, one for each antibody

# ginormous_df <- subset(ginormous_df, ginormous_df$id %in% kids_with_complete_timecourses$id)

# list_of_dfs <- split(ginormous_df, ginormous_df$antibody)

# run a model for each with id as random effect
# list_of_models <- lapply(list_of_dfs, function(x) lmer(log_conc~timepoint+(1|id), data=x))

# list_of_tests <- lapply(list_of_models, function(x) multcomp::glht(x, sec_ter_contrast))
# list_of_pvalues <- sapply(list_of_tests, function(x) summary(x)$test$pvalues)
# #
# #
# list_of_adj_pvalues <- sort(p.adjust(list_of_pvalues, method = "fdr"))
# 
# table(list_of_adj_pvalues<0.05) #all sic at 2nd and 3rd relative to first; 12/22 from 2nd to 3rd
# sig_2_3_abs <- subset(list_of_adj_pvalues, list_of_adj_pvalues<0.05)
# 
# # make results table
# list_of_coefficients <- lapply(list_of_models, function(x)coef(x)$id[1,])
# coefficient_df <- do.call(rbind,list_of_coefficients)
# coefficient_df$ab <- rownames(coefficient_df)
# coefficient_df$raw_p <- list_of_pvalues[coefficient_df$ab]
# coefficient_df$p_adj <- list_of_adj_pvalues[coefficient_df$ab]

# 
# tp1_boolean <- as.numeric(malaria_episodes$time_to_tp1) < 0 & as.numeric(malaria_episodes$time_to_tp1) > -30
# tp2_boolean <- as.numeric(malaria_episodes$time_to_tp2) < 0 & as.numeric(malaria_episodes$time_to_tp2) > -30
# tp3_boolean <- as.numeric(malaria_episodes$time_to_tp3) < 0 & as.numeric(malaria_episodes$time_to_tp3) > -30
# 
# tp1_boolean <- ifelse(is.na(tp1_boolean), FALSE, tp1_boolean)
# tp2_boolean <- ifelse(is.na(tp2_boolean), FALSE, tp2_boolean)
# tp3_boolean <- ifelse(is.na(tp3_boolean), FALSE, tp3_boolean)
# 
#  malaria_episodes_near_sample_dates <- malaria_episodes %>%
#                                         filter(between(as.numeric(time_to_tp3), -60, 0))

                                        
# reducing data to be within reference range ####


#subset data to include ID, timepoint, all antigens, flags and concentrations
raw_data <- bc1[,c(1,67:132, 155)]

#shove all antigens in their own column, same for flags and concentrations
raw_conc_data <- raw_data %>%
  pivot_longer(cols = colnames(raw_data)[seq(3, 66, by=3)], values_to = "conc")%>%
  dplyr::select(conc)

raw_flag_data <- raw_data%>%
  pivot_longer(cols = colnames(raw_data)[seq(4, 67, by=3)], values_to = "flag")%>%
  dplyr::select(flag)

long_raw_df <- raw_data %>%
  pivot_longer(cols = colnames(raw_data)[seq(2, 65, by=3)], values_to = "antigen") %>%
  dplyr::select(id, timepoint, antigen)%>%
  mutate(timepoint=factor(timepoint))%>%
  mutate(conc=raw_conc_data$conc, flag=raw_flag_data$flag)%>%
  mutate(antigen=gsub(".", " ", antigen, fixed = TRUE))%>%
  mutate(antigen=gsub("  ", " ", antigen, fixed = TRUE))%>%
  mutate(id=factor(id))%>%
  mutate(flag_type = case_when(flag==1 ~ "AboveMaxStd",
                               flag==2 ~ "AboveUpperbound",
                               flag==3 ~ "Above_fitted_asymptote",
                               flag==4 ~ "BelowMinStd",
                               is.na(flag) ~ "No Flag"))%>%
  filter(antigen != "")



# AboveMaxStd        AboveUpperbound Above_fitted_asymptote            BelowMinStd 
# 1                      2                      3                      4

flag_dist <- long_raw_df %>%
  count(antigen, flag) %>%
  pivot_wider(names_from = flag, values_from = n, names_prefix = "F") %>%
  group_by(antigen)%>%
  mutate(any_flag = sum(F1,F2,F3,F4, na.rm = TRUE))%>%
  mutate(any_flag_perc= any_flag/sum(any_flag, FNA))

flag_dist$antigen <- factor(flag_dist$antigen, levels=c(flag_dist$antigen[order(flag_dist$any_flag_perc, decreasing = TRUE)]))
 

flagless_df <- long_raw_df %>%
  filter(is.na(flag))


flag_dist_data <- long_raw_df %>%
  group_by(antigen, timepoint, flag_type) %>%
  summarise(flag_n=n())
  
  


flag_dist_plot <- ggplot(flag_dist_data, aes(x=timepoint, y=flag_n, fill=flag_type))+
  geom_bar(stat="identity", position = position_stack())+
  theme_minimal()+
  facet_wrap(~antigen)+
  # scale_fill_manual(values = pc1_cols)+
  theme(axis.title = element_blank(),
        strip.text = element_text(size=8),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/flag_distribution.png", flag_dist_plot, height=4, width=7.5, bg="white")





conc_by_flag <- ggplot(long_raw_df, aes(x=timepoint, y=conc, color=flag_type, group=flag_type))+
  geom_point(position = position_dodge(width=1))+
  facet_wrap(~antigen)+
  theme_minimal()+
  scale_y_log10()+
  theme(axis.title = element_blank(),
        strip.text = element_text(size=8),
        legend.title = element_blank())
  
ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/conc_by_flag.png", conc_by_flag, height=6, width=7.5, bg="white")



flag_dist_plot <- ggplot(flag_dist_data, aes(x=timepoint, y=flag_n, fill=flag_type))+
  geom_bar(stat="identity", position = position_stack())+
  theme_minimal()+
  facet_wrap(~antigen)+
  # scale_fill_manual(values = pc1_cols)+
  theme(axis.title = element_blank(),
        strip.text = element_text(size=8),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/flag_distribution.png", flag_dist_plot, height=4, width=7.5, bg="white")






flagless_ab_histogram <- ggplot(flagless_df, aes(x=conc, color=antigen, fill=antigen))+
  geom_histogram()+
  scale_x_log10()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  facet_wrap(~antigen)+
  theme_minimal()+
  theme(legend.position="none",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/flagless_ab_histogram.png", flagless_ab_histogram, height=6, width=8, bg="white")


number_of_observations_heatmap <- flagless_df %>%
  filter(!is.na(conc))%>%
  group_by(antigen, timepoint)%>%
  summarise(number_of_observations=n())%>%
  mutate(cols=if_else(number_of_observations>100, "black", "white"))%>%
  ggplot(., aes(x=timepoint, y=antigen, label=number_of_observations))+
  geom_tile(aes(fill=number_of_observations))+
  geom_text(aes(color=cols))+
  scale_fill_viridis_c(option = "B")+
  scale_color_manual(values=c("black", "white"))+
  theme_minimal()+
  scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=12),
        legend.position = "none",
        panel.grid = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/number_of_observations_heatmap.png", number_of_observations_heatmap, height=8, width=4, bg="white")


flagless_boxplot_timepoint <- ggplot(flagless_df, aes(x=timepoint, y=conc, color=antigen, fill=antigen))+
  # geom_boxplot()+
  geom_violin()+
  scale_y_log10()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  facet_wrap(~antigen)+
  theme_minimal()+
  theme(legend.position="none",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/flagless_ab_through_time.png", flagless_boxplot_timepoint, height=6, width=8, bg="white")




modelable_antigens_plot <- flagless_df%>%
  filter(antigen %in% modelable_antigens)%>%
  group_by(antigen, timepoint)%>%
  add_count(name = "n_observations")%>%
  ggplot(., aes(x=timepoint, y=conc, color=antigen, fill=antigen))+
  #geom_boxplot()+
  geom_violin()+
  geom_text(aes(y=100, label=n_observations))+
  scale_y_log10()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  facet_wrap(~antigen)+
  theme_minimal()+
  theme(legend.position="none",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/modelable_antigens_through_time.png", modelable_antigens_plot, height=6, width=8, bg="white")





modelable_antigens_plot <- flagless_df%>%
  filter(antigen %in% modelable_antigens)%>%
  group_by(antigen, timepoint)%>%
  add_count(name = "n_observations")%>%
  ggplot(., aes(x=timepoint, y=conc, color=antigen, fill=antigen))+
  #geom_boxplot()+
  geom_violin()+
  geom_text(aes(y=100, label=n_observations))+
  scale_y_log10()+
  scale_color_manual(values=pc1_cols)+
  scale_fill_manual(values=pc1_cols)+
  facet_wrap(~antigen)+
  theme_minimal()+
  theme(legend.position="none",
        axis.title = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/modelable_antigens_through_time.png", modelable_antigens_plot, height=6, width=8, bg="white")




# define conrtrasts
sec_contrast <- t(matrix(c(0,1,0)))
ter_contrast <- t(matrix(c(0,0,1)))
sec_ter_contrast <- t(matrix(c(0,-1,1)))

# brave new world purrrlicious way of doing it:
# make grouped df, nest into one, add column for model & summary, extract values from each to create even more new columns for p values and p_adj values
purrrf <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  group_by(antigen) %>%
  filter(!is.na(conc))%>%
  nest() %>%
  mutate(model=map(data, ~lmer(log10(conc)~timepoint+(1|id), data=., REML = FALSE)))%>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(t2_t1=map(model, ~multcomp::glht(., sec_contrast)),
         t2_t1_p=map_dbl(t2_t1, ~summary(.)$test$pvalues),
         t2_t1_p_adj=map_dbl(t2_t1_p, ~p.adjust(.))) %>%
  mutate(t3_t1=map(model, ~multcomp::glht(., ter_contrast)),
         t3_t1_p=map_dbl(t3_t1, ~summary(.)$test$pvalues),
         t3_t1_p_adj=map_dbl(t3_t1_p, ~p.adjust(.))) %>%
  mutate(t3_t2=map(model, ~multcomp::glht(., sec_ter_contrast)),
         t3_t2_p=map_dbl(t3_t2, ~summary(.)$test$pvalues),
         t3_t2_p_adj=map_dbl(t3_t2_p, ~p.adjust(.)))

results <- purrrf %>%
  mutate(t3_t2_coef=10^(map_dbl(model, ~-coef(.)$id[1,2]+coef(.)$id[1,3])))%>%
  mutate(t3_t2_coef=round(t3_t2_coef, digits = 2 ), t3_t2_p_adj=round(t3_t2_p_adj, digits = 2))%>%
  select(antigen, t3_t2_coef, t3_t2_p_adj) %>%
  ungroup()


kbl(results, booktabs = T, linesep = "", escape=FALSE) %>% 
  kable_paper(full_width = F) %>%
  column_spec(3, color = ifelse(results$t3_t2_p_adj < 0.05, "red", "black"))%>%
  save_kable(file = "~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/regression_results.html", self_contained = T)



results_table$t3_t2_p_adj <- ifelse(
  results_table$t3_t2_p_adj < 0.05,
    cell_spec(results_table$t3_t2_p_adj, color = "red"),
    cell_spec(results_table$t3_t2_p_adj, color = "black")
  )

kbl(results_table) %>%
  kable_paper() %>%
  save_kable(file = "~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/regression_results.html", self_contained = T)


wide_flagless <- flagless_df %>%
  filter(antigen %in% c("CSP GENOVA", "GEXP", "PfSEA", "Tet Tox ", "Rh5"))%>%
  select(-flag)%>%
  filter(!is.na(conc))%>%
  pivot_wider(names_from = antigen, values_from = conc) %>%
  filter(timepoint==1)

no_na_flagless <- na.omit(wide_flagless)



t1_pca <-  prcomp(no_na_flagless[,c("CSP GENOVA", "GEXP", "PfSEA", "Tet Tox ", "Rh5")], center = T)
t1_pca_plot_data <- as.data.frame(cbind(no_na_flagless, t1_pca$x))


loadings_df <- data.frame(t1_pca$rotation)
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

t1_demo_data <- bc1 %>%
  filter(timepoint==1 & id %in% t1_pca_plot_data$id)%>%
  mutate(sample_id=paste(id, timepoint, sep="_"))%>%
  select(c(colnames(bc1)[2:55]))

t1_pca_plot_data <- cbind(t1_demo_data, t1_pca_plot_data)

#cont

ggplot(t1_pca_plot_data, aes(x=PC1, y=PC2, color=covarate))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(t1_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(t1_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")
  


# individual-level plots ####

hi <-  bc1 %>%
  mutate(sample_id=paste(id, timepoint, sep="_"))


flagless_df <- long_raw_df %>%
  filter(is.na(flag))


flagless_df$sample_id <- paste(flagless_df$id, flagless_df$timepoint, sep="_")

demo_columns2 <- c("malaria","mal0to6","febrile0to6", "incidentmalaria", "febrile6to12", "mal6to12", "mbloodLAMP")
flagless_df[, demo_columns2] <- hi[match(flagless_df$sample_id, hi$sample_id), demo_columns2]

mal0_6_plot <- flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(mal0to6)%>%
  ggplot(aes(x=id, y=conc, color=mal0to6, group=id))+
  geom_point()+
  scale_y_log10()+
  # geom_line()+
  facet_grid(timepoint~antigen)+
  theme_minimal()+
  theme(axis.text.x = element_blank())+
  viridis::scale_color_viridis(option="B", direction = -1)

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/mal0_6_plot.png", mal0_6_plot, width = 10, height = 6, bg="white")

mal6_12_plot <-flagless_df %>%
  filter(antigen %in% modelable_antigens)%>%
  arrange(mal6to12)%>%
  ggplot(aes(x=id, y=conc, color=mal6to12, group=id))+
  geom_point()+
  scale_y_log10()+
  # geom_line()+
  facet_grid(timepoint~antigen)+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank()
  )+
  viridis::scale_color_viridis(option="B", direction = -1)

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/mal6_12_plot.png", mal6_12_plot, width = 10, height = 6, bg="white")




