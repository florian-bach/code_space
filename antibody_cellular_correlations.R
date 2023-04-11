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
clin_data <- clin_data %>%
  filter(id %in% long_raw_dfff$id)

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

# for these kids we only have antibody measurements at birth and they're not in the clinical database so we'll cut them her
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
                           column_names_gp = gpar(fontsize = 6),
                           row_names_gp = gpar(fontsize = 6),
                           row_names_side = "left",
                           col = col_fun_pearson,
                           column_names_rot = 45,
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                           })


# statsy correlations


cleaned_broom_data <- twelve_cor_combo %>%
  filter(tfh_drop=="No", abc_drop=="No") %>%
  dplyr::select(-tfh_drop, -abc_drop) %>%
  pivot_longer(cols = contains("_"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = colnames(.)[4:24], names_to = "antigen", values_to = "ab_conc")


cleaned_broomer <- cleaned_broom_data %>%
  group_by(cell_pop, antigen)%>%
  do(broom::tidy(cor.test(.$cell_freq, .$ab_conc, method="spearman")))%>%
  mutate("p_adj"=p.adjust(p.value))%>%
  ungroup()

sig_cleaned_broomer <- filter(cleaned_broomer, p_adj<0.1)



cleaned_broom_data %>%
  filter(cell_pop=="cd10_neg_b", antigen %in% sig_cleaned_broomer$antigen)%>%
  ggplot(., aes(x=cell_freq, y=ab_conc))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

#best hit
cleaned_broom_data %>%
  filter(cell_pop=="tfh_q2", antigen =="MSP1")%>%
  ggplot(., aes(x=cell_freq, y=ab_conc))+
  geom_point(aes(color=antigen))+
  geom_smooth(method="lm")+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


# homework: do proper models, not just correlations; redo purr models of incidence using cellular frequencies
