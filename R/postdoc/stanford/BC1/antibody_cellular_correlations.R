#libraries 
library(patchwork)
library(tidyr)
library(ggplot2)
library(purrr)
library(MASS)
library(knitr)
library(dplyr)


#palettes
time_palette <- colorspace::sequential_hcl(n=4, "RdPu")[1:3]
pc1_cols <- colorspace::sequential_hcl(23, palette = "Purple Yellow")
incidence_cols <- colorspace::sequential_hcl(11, palette = "Purple Yellow")
n_infection_cols <- c("white", colorspace::sequential_hcl(n=5, palette = "Lajolla")[-1])


`%notin%` <- Negate(`%in%`)

fdr_cutoff <- 0.1

# from looking at the frequencies of observations only these antibodies can be modelled longitudinally; NB some drop off at timepoint 3 so only 1 and 2 can be compared properly
modelable_antigens <- c("Tet Tox", "SBP1", "Rh5", "PfSEA", "PfAMA1", "Hyp2", "HSP40 Ag1", "GST", "GEXP", "CSP GENOVA")


# read in antibody data
bc1 <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

ab_columns <- grep("log", colnames(bc1), value = TRUE)


long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "anyHPfinal")))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))%>%
  mutate(antigen=gsub("log", "", antigen, fixed = TRUE))%>%
  mutate(antigen=gsub("_", " ", antigen, fixed = TRUE))%>%
  mutate(MomFinalRxx=if_else(MomFinalRx==1, "3 Dose SP",
                             ifelse(MomFinalRx==2, "3 Dose DP",
                                    if_else(MomFinalRx==3, "Monthly DP", "NA")))
         
  )%>%
  mutate(MomFinalRxx=factor(MomFinalRxx, levels = c("3 Dose SP", "3 Dose DP", "Monthly DP")))%>%
  mutate(anyHPfinalx=if_else(anyHPfinal==1, "Placental Malaria",
                             if_else(anyHPfinal==0, "No Pathology", "Results missing")))



# read in new visits database to look at correlations with malaria incidence
# clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_ALL.dta")
clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_REV.dta")

# it's a big table, so let's only include kids for whom we have any antibody measurements
# clin_data <- clin_data %>%
#   filter(id %in% long_raw_dfff$id)

# make a data frame where we add a bunch of columns that contain how many (symptomatic) infections were experienced by the child in the indicated time window

infs <- clin_data %>%
  group_by(id) %>%
  dplyr::select(id, age, anyinfection, sxifinfected) %>%
  mutate(inf_0_6   = sum(if_else(age<6, anyinfection, 0), na.rm = TRUE),
         inf_6_12  = sum(if_else(age>6  & age<12, anyinfection, 0), na.rm = TRUE),
         inf_12_18 = sum(if_else(age>12 & age<18, anyinfection, 0), na.rm = TRUE),
         inf_12_24 = sum(if_else(age>12 & age<24, anyinfection, 0), na.rm = TRUE),
         inf_0_12  = sum(if_else(age<12, anyinfection, 0), na.rm = TRUE),
         inf_0_24  = sum(if_else(age<24, anyinfection, 0), na.rm = TRUE),
         symp_0_6   = sum(if_else(age<6, sxifinfected, 0), na.rm = TRUE),
         symp_6_12  = sum(if_else(age>6  & age<12, sxifinfected, 0), na.rm = TRUE),
         symp_12_18 = sum(if_else(age>12 & age<18, sxifinfected, 0), na.rm = TRUE),
         symp_12_24 = sum(if_else(age>12 & age<24, sxifinfected, 0), na.rm = TRUE),
         symp_0_12  = sum(if_else(age<12, sxifinfected, 0), na.rm = TRUE),
         symp_0_24  = sum(if_else(age<24, sxifinfected, 0), na.rm = TRUE),
         any_inf_0_6 = ifelse(inf_0_6==0, 0, 1),
         any_inf_6_12 = ifelse(inf_6_12==0, 0, 1),
         any_inf_12_18 = ifelse(inf_12_18==0, 0, 1),
         any_symp_0_6 = ifelse(symp_0_6==0, 0, 1),
         any_symp_6_12 = ifelse(symp_6_12==0, 0, 1),
         any_symp_12_18 = ifelse(symp_12_18==0, 0, 1)
  ) %>%
  dplyr::select(-anyinfection, -sxifinfected, -age) %>%
  distinct()


# combine antibody data with malaria incidence data
# the -1 removes the id column from the infs df, otherwise it's duplicated         
combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])

# for these kids we only have antibody measurements at birth and they're not in the clinical database so we'll cut them here
combo_data <- filter(combo_data, id %notin% c(11130, 11084, 11037))

raw_combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  dplyr::select(id)

combo_data <- filter(combo_data, id %in% kids_with_complete_timecourses$id)


#read in tfh data 

tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
abc <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/atypical_b.csv")

#fix & harmonise colnames, add id column that matches rest of data
cleaned_tfh <- tfh
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
cleaned_tfh$idt <- paste("1", substr(cleaned_tfh$idt, 15, 18), sep='')

#fix & harmonise colnames, add id column that matches rest of data
cleaned_abc <- abc
colnames(cleaned_abc) <- c("idb", "abc_q4", "abc_q3","abc_q2","abc_q1","cd10_neg_b", "parent2b", "abc_drop")
cleaned_abc$idb <- paste("1", substr(cleaned_abc$idb, 17, 20), sep='')
# fcs filename has a typo need to fix this
cleaned_abc$idb <- gsub("11108", "11083", cleaned_abc$idb)

# all good, ready for merge
combo_cells_data <- cbind(cleaned_tfh, cleaned_abc[match(cleaned_tfh$idt, cleaned_abc$idb),])
# all(combo_cells_data$idt == combo_cells_data$idb)
# [1] TRUE
combo_cells_data <- combo_cells_data %>%
  mutate("id"=as.numeric(idt))|>
  dplyr::select(-idt, -idb, -parent2t, -parent2b)

long_combo_combo <- full_join(combo_data, combo_cells_data, by="id")
write.csv(long_combo_combo, "~/postdoc/stanford/clinical_data/BC1/cell_inf_ab_database.csv", quote=FALSE, row.names=FALSE)

wide_combo_12months <- combo_data |> 
  dplyr::select(id, timepoint, anyHPfinal, antigen, conc)|>
  filter(timepoint==3)|>
  pivot_wider(names_from = antigen, values_from = conc)

twelve_cor_combo <- inner_join(wide_combo_12months, combo_cells_data, by="id")


# big matrix

cor_data <- twelve_cor_combo |>
  dplyr::select(-id, -timepoint, -anyHPfinal, -tfh_drop, -abc_drop)%>%
  na.omit()
  

spearman <- cor(cor_data, method = "spearman")
twelve_dist <- dist(spearman, method = "euclid", diag = FALSE, upper = FALSE, p = 2)
twelve_hclust <- hclust(twelve_dist)

cor_mat <- as.matrix(spearman)

#order matrix according to clustering results
# cor_mat <- as.matrix(
#   cor_data[
#     rownames(cor_data)[rev(twelve_hclust$order)],colnames(cor_data)[rev(twelve_hclust$order)]
#     ]
#   )


# col_fun_pearson <- circlize::colorRamp2(c(min(cor_mat), 0, max(cor_mat)), c("#0859C6", "black", "#FFA500"))

col_fun_pearson <- circlize::colorRamp2(c(-1, 0, 1), c("#0859C6", "black", "#FFA500"))

#doesnt work for pdf..
# colnames(pearson_matrix) <- gsub("IFNy", "IFNγ", colnames(pearson_matrix))
# rownames(pearson_matrix) <-  gsub("IFNy", "IFNγ",rownames(pearson_matrix))

library(ComplexHeatmap)

pearson_heatmap <- Heatmap(matrix = cor_mat,
                           cluster_rows = TRUE,
                           cluster_columns=TRUE,
                           show_row_dend = FALSE,
                           show_column_dend = TRUE,
                           show_heatmap_legend = TRUE,
                           name = "Spearmon rho",
                           #name = "Pearson r",
                           #cluster_columns = FALSE,
                           column_names_gp = gpar(fontsize = 9),
                           row_names_gp = gpar(fontsize = 9),
                           row_names_side = "left",
                           col = col_fun_pearson,
                           column_names_rot = 45,
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 5))
                           })


png("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/big_correlation_heatmap.png", width = 8, height=8, units="in", res=700)
draw(pearson_heatmap)
dev.off()

# statsy correlations


cleaned_broom_data <- twelve_cor_combo %>%
  #filter(tfh_drop=="No", abc_drop=="No") %>%
  dplyr::select(-tfh_drop, -abc_drop) %>%
  pivot_longer(cols = contains("_"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = colnames(.)[4:24], names_to = "antigen", values_to = "ab_conc")







cleaned_broomer <- cleaned_broom_data %>%
  group_by(cell_pop, antigen)%>%
  do(broom::tidy(cor.test(.$cell_freq, .$ab_conc, method="spearman")))%>%
  group_by(cell_pop)%>%
  mutate("p_adj"=p.adjust(p.value))%>%
  ungroup()

sig_cleaned_broomer <- filter(cleaned_broomer, p.value<0.05)



corr_figure <- cleaned_broom_data %>%
  filter(cell_pop=="cd10_neg_b", antigen %in% sig_cleaned_broomer$antigen)%>% 
  mutate("immature_bcell"=100-cell_freq)%>%
  ggplot(., aes(x=immature_bcell, y=ab_conc))+
  ggtitle("almost significant correlations (Rho~0.33, p_adj=0.138)")+
  geom_point(aes(color=antigen))+
  xlab("percentage of immature B cells")+
  ylab("antibody intensity")+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_smooth(method="lm")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/bcell_ab_correlations.png", corr_figure, height=4, width=8, dpi=444, bg="white")


#best hit
cleaned_broom_data %>%
  filter(cell_pop=="tfh_q2", antigen =="MSP1")%>%
  ggplot(., aes(x=cell_freq, y=ab_conc))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


# homework: do proper models, not just correlations; redo purr models of incidence using cellular frequencies


# only 
inf_cell_combo <- left_join(combo_cells_data, infs, by="id")

long_inf_cell_combo <- inf_cell_combo%>%
  dplyr::select(-tfh_drop, -abc_drop) %>%
  #filter(tfh_drop=="No", abc_drop=="No") %>%
  pivot_longer(cols = matches('tfh|abc|cd10'), names_to = "cell_pop", values_to = "cell_freq")

combo_combo_purrrf <- long_inf_cell_combo%>%
  na.omit()%>%
  group_by(cell_pop) %>%
  nest() %>%
  mutate(model_12_18=map(data, ~glm.nb(inf_12_18 ~ cell_freq, data=.)))%>%
  mutate(model_12_18_symp=map(data, ~glm.nb(symp_12_18 ~ cell_freq, data=.)))%>%
  mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~cell_freq, data=., family = "binomial", weights = inf_12_18)))%>%
  #mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ cell_freq, data=.)))%>%
  mutate(model_6_12 =map(data, ~glm(any_symp_12_18~cell_freq, data=., family = "poisson")))%>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ cell_freq, data=.)))%>%
  mutate(model_any_symp_12_18=map(data, ~glm(any_symp_12_18~cell_freq, data=., family = "binomial")))%>%
  mutate(model_any_12_18=map(data, ~glm(any_inf_12_18~cell_freq, data=., family = "binomial")))%>%
  
  mutate(summary_12_18=map(model_12_18, ~summary(.))) %>%
  mutate(summary_12_18_symp=map(model_12_18_symp, ~summary(.))) %>%
  mutate(summary_12_18_symp_prob=map(model_12_18_symp_prob, ~summary(.))) %>%
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  mutate(summary_any_symp_12_18=map(model_any_symp_12_18, ~summary(.))) %>%
  mutate(summary_any_12_18=map(model_any_12_18, ~summary(.))) %>%
  
  mutate(summary_12_18_p=map_dbl(summary_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_p=map_dbl(summary_12_18_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_prob_p=map_dbl(summary_12_18_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_symp_12_18_p=map_dbl(summary_any_symp_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_12_18_p=map_dbl(summary_any_12_18, ~unlist(.$coefficients[8])))%>%
  
  mutate(model_12_24=map(data, ~glm.nb(inf_12_24 ~ cell_freq, data=.)))%>%
  mutate(model_12_24_symp=map(data, ~glm.nb(symp_12_24 ~ cell_freq, data=.)))%>%
  mutate(model_12_24_symp_prob=map(data, ~glm(symp_12_24/inf_12_24~cell_freq, data=., family = "binomial", weights = inf_12_24)))%>%
  
  mutate(summary_12_24=map(model_12_24, ~summary(.))) %>%
  mutate(summary_12_24_symp=map(model_12_24_symp, ~summary(.))) %>%
  mutate(summary_12_24_symp_prob=map(model_12_24_symp_prob, ~summary(.))) %>%
  
  mutate(summary_12_24_p=map_dbl(summary_12_24, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_p=map_dbl(summary_12_24_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_prob_p=map_dbl(summary_12_24_symp_prob, ~unlist(.$coefficients[8])))%>%
  ungroup()%>%
  mutate(summary_12_18_padj= p.adjust(summary_12_18_p))%>%
  mutate(summary_12_18_symp_padj=p.adjust(summary_12_18_symp_p))%>%
  mutate(summary_12_18_symp_prob_padj=p.adjust(summary_12_18_symp_prob_p))%>%
  mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
  mutate(summary_12_24_padj=p.adjust(summary_12_24_p))%>%
  mutate(summary_12_24_symp_padj=p.adjust(summary_12_24_symp_p))%>%
  mutate(summary_12_24_symp_prob_padj=p.adjust(summary_12_24_symp_prob_p))%>%
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
  mutate(summary_any_symp_12_18_padj=p.adjust(summary_any_symp_12_18_p))%>%
  mutate(summary_any_12_18_padj=p.adjust(summary_any_12_18_p))



fdr_cutoff <- 0.1


twelve_18_sig <- combo_combo_purrrf %>%
  filter(summary_12_18_padj<fdr_cutoff)%>%
  dplyr::select(cell_pop, summary_12_18_padj)

twelve_18_symp_sig <- combo_combo_purrrf %>%
  filter(summary_12_18_symp_padj<fdr_cutoff)%>%
  dplyr::select(cell_pop, summary_12_18_symp_padj)

twelve_18_symp_prob_sig <- combo_combo_purrrf %>%
  filter(summary_12_18_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(cell_pop, summary_12_18_symp_prob_padj)

twelve_6_12_sig <- combo_combo_purrrf %>%
  filter(summary_6_12_padj<fdr_cutoff)%>%
  dplyr::select(cell_pop, summary_6_12_padj)


padj_combo_combo_purrrf <-  combo_combo_purrrf %>%
  pivot_longer(cols = ends_with('_p'), names_to = "model", values_to = "raw_p")%>%
  dplyr::select(cell_pop, model, raw_p)%>%
  mutate("padj"=p.adjust(raw_p))%>%
  filter(raw_p<0.1)



long_inf_cell_combo%>%
  filter(cell_pop=="tfh_q3")%>%
  na.omit()%>%
  ggplot(., aes(x=any_symp_12_18, y=cell_freq))+
  geom_point(alpha=0.2, na.rm = TRUE)+
  # xlab("Number of Infections in Months 12-18")+
  # ylab("Cell Frequency at 12 Months")+
  ggpubr::stat_cor(method = "spearman", label.x = 1.5, label.y = 1, na.rm = TRUE)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")#+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)

try <- long_combo_combo%>%
  filter(timepoint==3)%>%
  filter(tfh_drop=="No", abc_drop=="No") %>%
  dplyr::select(-antigen, -conc, -tfh_drop, -abc_drop) %>%
  distinct(id, timepoint, .keep_all = TRUE)%>%
  pivot_longer(cols = matches('tfh|abc|cd10'), names_to = "cell_pop", values_to = "cell_freq")




# summary level correlations
all_ab_broom_data <- cleaned_broom_data%>%
  group_by(id)%>%
  summarise("all_abs"=log10(sum(10^ab_conc)), cell_pop, cell_freq)%>%
  distinct()%>%
  ungroup()

all_ab_broomer <- all_ab_broom_data%>%
  group_by(cell_pop)%>%
  do(broom::tidy(cor.test(.$cell_freq, .$all_abs, method="spearman")))%>%
  ungroup()%>%
  mutate("p_adj"=p.adjust(p.value, method="fdr" ))


sig_all_ab_cleaned_broomer <- filter(all_ab_broomer, p_adj<0.1)


all_ab_broom_data %>%
  filter(cell_pop %in% "cd10_neg_b")%>% 
  ggplot(., aes(x=cell_freq, y=all_abs))+
  geom_point(aes(color=cell_pop))+
  geom_smooth(method="lm")+
  facet_wrap(~cell_pop, scales="free")+
  theme_minimal()



all_ab_broom_data <- cleaned_broom_data%>%
  group_by(id)%>%
  summarise("all_abs"=log10(sum(10^ab_conc)), cell_pop, cell_freq)%>%
  distinct()%>%
  ungroup()




all_ab_inf_purf <- combo_data %>%
  filter(inf_12_18!=0, timepoint==3)%>%
  group_by(id)%>%
  mutate("all_abs"=log10(sum(10^conc, na.rm = TRUE)))%>%
  dplyr::select(-antigen, -conc)%>%
  distinct()%>%
  ungroup()%>%
  nest(data=everything()) %>%
  mutate(model_12_18=map(data, ~glm.nb(inf_12_18 ~ all_abs, data=.)))%>%
  mutate(model_12_18_symp=map(data, ~glm.nb(symp_12_18 ~ all_abs, data=.)))%>%
  mutate(model_12_18_symp_prob=map(data, ~glm(symp_12_18/inf_12_18~all_abs, data=., family = "binomial", weights = inf_12_18)))%>%
  mutate(model_6_12=map(data, ~glm.nb(inf_6_12 ~ all_abs, data=.)))%>%
  mutate(model_0_12=map(data, ~glm.nb(inf_0_12 ~ all_abs, data=.)))%>%
  mutate(model_any_symp_12_18=map(data, ~glm(any_symp_12_18~all_abs, data=., family = "binomial")))%>%
  mutate(model_any_12_18=map(data, ~glm(any_inf_12_18~all_abs, data=., family = "binomial")))%>%
  
  mutate(summary_12_18=map(model_12_18, ~summary(.))) %>%
  mutate(summary_12_18_symp=map(model_12_18_symp, ~summary(.))) %>%
  mutate(summary_12_18_symp_prob=map(model_12_18_symp_prob, ~summary(.))) %>%
  mutate(summary_6_12=map(model_6_12, ~summary(.))) %>%
  mutate(summary_0_12=map(model_0_12, ~summary(.))) %>%
  mutate(summary_any_symp_12_18=map(model_any_symp_12_18, ~summary(.))) %>%
  mutate(summary_any_12_18=map(model_any_12_18, ~summary(.))) %>%
  
  mutate(summary_12_18_p=map_dbl(summary_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_p=map_dbl(summary_12_18_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_18_symp_prob_p=map_dbl(summary_12_18_symp_prob, ~unlist(.$coefficients[8])))%>%
  mutate(summary_6_12_p=map_dbl(summary_6_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_0_12_p=map_dbl(summary_0_12, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_symp_12_18_p=map_dbl(summary_any_symp_12_18, ~unlist(.$coefficients[8])))%>%
  mutate(summary_any_12_18_p=map_dbl(summary_any_12_18, ~unlist(.$coefficients[8])))%>%
  
  mutate(model_12_24=map(data, ~glm.nb(inf_12_24 ~ all_abs, data=.)))%>%
  mutate(model_12_24_symp=map(data, ~glm.nb(symp_12_24 ~ all_abs, data=.)))%>%
  mutate(model_12_24_symp_prob=map(data, ~glm(symp_12_24/inf_12_24~all_abs, data=., family = "binomial", weights = inf_12_24)))%>%
  
  mutate(summary_12_24=map(model_12_24, ~summary(.))) %>%
  mutate(summary_12_24_symp=map(model_12_24_symp, ~summary(.))) %>%
  mutate(summary_12_24_symp_prob=map(model_12_24_symp_prob, ~summary(.))) %>%
  
  mutate(summary_12_24_p=map_dbl(summary_12_24, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_p=map_dbl(summary_12_24_symp, ~unlist(.$coefficients[8])))%>%
  mutate(summary_12_24_symp_prob_p=map_dbl(summary_12_24_symp_prob, ~unlist(.$coefficients[8])))%>%
  ungroup()%>%
  mutate(summary_12_18_padj= p.adjust(summary_12_18_p))%>%
  mutate(summary_12_18_symp_padj=p.adjust(summary_12_18_symp_p))%>%
  mutate(summary_12_18_symp_prob_padj=p.adjust(summary_12_18_symp_prob_p))%>%
  mutate(summary_6_12_padj=p.adjust(summary_6_12_p))%>%
  mutate(summary_12_24_padj=p.adjust(summary_12_24_p))%>%
  mutate(summary_12_24_symp_padj=p.adjust(summary_12_24_symp_p))%>%
  mutate(summary_12_24_symp_prob_padj=p.adjust(summary_12_24_symp_prob_p))%>%
  mutate(summary_0_12_padj=p.adjust(summary_0_12_p))%>%
  mutate(summary_any_symp_12_18_padj=p.adjust(summary_any_symp_12_18_p))%>%
  mutate(summary_any_12_18_padj=p.adjust(summary_any_12_18_p))






twelve_18_sig <- twelve_purf %>%
  filter(summary_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_padj)

twelve_18_symp_sig <- twelve_purf %>%
  filter(summary_12_18_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_symp_padj)

twelve_18_symp_prob_sig <- twelve_purf %>%
  filter(summary_12_18_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_18_symp_prob_padj)

twelve_6_12_sig <- twelve_purf %>%
  filter(summary_6_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_6_12_padj)


twelve_24_sig <- twelve_purf %>%
  filter(summary_12_24_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_padj)

twelve_24_symp_sig <- twelve_purf %>%
  filter(summary_12_24_symp_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_symp_padj)

twelve_24_symp_prob_sig <- twelve_purf %>%
  filter(summary_12_24_symp_prob_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_12_24_symp_prob_padj)


twelve_0_12_sig <- twelve_purf %>%
  filter(summary_0_12_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_0_12_padj)

twelve_18_any_sig <- twelve_purf %>%
  filter(summary_any_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_12_18_padj)

twelve_18_any_symp_sig <- twelve_purf %>%
  filter(summary_any_symp_12_18_padj<fdr_cutoff)%>%
  dplyr::select(antigen, summary_any_symp_12_18_padj)


# make complex figure with all significant associations before multiple testing correction

# x_and_ys <- cbind(sig_cleaned_broomer$cell_pop, paste("log", gsub(" ", "_", sig_cleaned_broomer$antigen), sep=""))

list_of_plots <- list(matrix(nrow = 14))

qual_palette <- viridis::inferno(14)

x_and_ys <- cbind(sig_cleaned_broomer$cell_pop, sig_cleaned_broomer$antigen)

for(i in 1:14){
  
  plot <- cleaned_broom_data %>%
    filter(cell_pop==x_and_ys[i,1], antigen==x_and_ys[i,2])%>%
    ggplot(., aes(x=cell_freq, y=ab_conc))+
    geom_point(color=pc1_cols[i])+
    geom_smooth(method="lm")+
    xlab(x_and_ys[i,1])+
    ylab(x_and_ys[i,2])+
    theme_minimal()

list_of_plots[[i]] <- plot

}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow = 3)

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/big_indie_ab_cell_correlation_plot.png", big_plot, height = 5.5, width=8, bg="white")
