# libraries & preamble ####
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

# antibody stuff ####
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



#clinical data ####
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


# no correlations between CXCR3+ memory CD4 T cells and antibodies
# no correlations between  CCR6+ memory CD4 T cells and antibodies
# no correlations between CXCR3+ memory CD4 T cells and infections (at least batch3 and batch4)
# no correlations between  CCR6+ memory CD4 T cells and infections (at least batch3 and batch4)

# import tfh data ####
# import final tfh data
tfh_batch1 <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/kareena_gating/12,2f,02,2f,16 tfh/batch1_stats_with_memory.csv")
tfh_batch2 <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/kareena_gating/14,2f,01,2f,16 tfh/batch2_stats_with_memory.csv")
tfh_batch3 <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/kareena_gating/16,2f,03,2f,16 tfh/flo_batch3_table_with_memory.csv")
tfh_batch4 <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/kareena_gating/15,2f,12,2f,15/batch4_stats_with_memory.csv")

colnames(tfh_batch1) <- c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive", "X.1")
colnames(tfh_batch2) <- c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive","X.1")
colnames(tfh_batch3) <- c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "CD4_Th1", "CD4_Th17", "Th_memory", "Th_naive", "X.1")
colnames(tfh_batch4) <- c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "CD4_Th1", "CD4_Th17", "Th_memory", "Th_naive", "X.1")


tfh_combo_batch <- rbind(tfh_batch1[,c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive")],
                     tfh_batch2[,c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive")],
                     tfh_batch3[,c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive")],
                     tfh_batch4[,c("X", "Tfh_count", "Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "Th_memory", "Th_naive")])

tfh_combo_batch <- tfh_combo_batch %>%
  filter(X!="SD", X!="Mean", !grepl("Control.fcs", tfh_combo_batch$X))%>%
  mutate(file_name=X, .keep = "unused")%>%
  mutate(id=as.numeric(paste(1, substr(.$file_name, 15, 18), sep="")))




# import abc data ####
abc_batch1 <- read.csv("~/postdoc/stanford/clinical_data/BC1/b_cell_data/box/11,2f,02,2f,16 AtMBC/batch1_stats.csv")
abc_batch2 <- read.csv("~/postdoc/stanford/clinical_data/BC1/b_cell_data/box/14,2f,01,2f,16 AtMBC/batch2_stats.csv")
abc_batch3 <- read.csv("~/postdoc/stanford/clinical_data/BC1/b_cell_data/box/17,2f,03,2f,16 AtMBC/batch3_stats.csv")

#messed up the order in flowjo, let's fix it for consistency
colnames(abc_batch1) <- c("X", "b_count", "mature_b_count", "activated_b", "memory_b", "naive_b", "atypical_b", "X.1")
colnames(abc_batch2) <- c("X", "b_count", "mature_b_count", "activated_b", "memory_b", "atypical_b", "naive_b", "X.1")
colnames(abc_batch3) <- c("X", "b_count", "mature_b_count", "activated_b", "naive_b", "memory_b", "atypical_b", "X.1")

#rbind, get order of columns smoothed out
abc_combo_batch <- rbind(abc_batch1[,c("X", "b_count", "mature_b_count", "activated_b", "memory_b", "naive_b", "atypical_b")],
                         abc_batch2[,c("X", "b_count", "mature_b_count", "activated_b", "memory_b", "naive_b", "atypical_b")],
                         abc_batch3[,c("X", "b_count", "mature_b_count", "activated_b", "memory_b", "naive_b", "atypical_b")])

abc_combo_batch <- abc_combo_batch %>%
  filter(X!="SD", X!="Mean", X!="B CELL PANEL_NF.fcs", !grepl("Control.fcs", abc_combo_batch$X))%>%
  mutate(file_name=X, .keep = "unused")%>%
  mutate(id=as.numeric(paste(1, substr(.$file_name, 17, 20), sep="")))



# combine tfh and abc data ####
half_combo <- inner_join(combo_data, tfh_combo_batch, by="id")
true_combo <- inner_join(half_combo, abc_combo_batch, by="id")


long_true_combo <- true_combo %>%
  mutate("immature_b_perc"=1-(mature_b_count/b_count), "any_0_12"=if_else(inf_0_12>0, "yes", "no"))%>%
  #filter(Tfh_count>500)%>%
  pivot_longer(cols = c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2", "activated_b", "memory_b", "naive_b", "atypical_b", "Th_memory", "Th_naive"), names_to = "cell_pop", values_to = "cell_freq")

# all_ab_broomer <- long_true_combo%>%
#   group_by(cell_pop)%>%
#   do(broom::tidy(cor.test(.$cell_freq, .$all_abs, method="spearman")))%>%
#   ungroup()%>%
#   mutate("p_adj"=p.adjust(p.value, method="fdr" ))


big_plot_12_24 <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive",
                         "immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"))%>%
  ggplot(., aes(x=factor(inf_12_24), y=cell_freq/100, fill=factor(inf_12_24)))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free", ncol=4)+
  # scale_fill_manual(values = n_infection_cols)+
  theme_minimal()


big_plot_0_12 <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive",
                         "immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"))%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100, fill=factor(inf_0_12)))+
  geom_boxplot()+
  scale_y_continuous(labels = scales::label_percent())+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free", ncol=4)+
  scale_fill_manual(values = n_infection_cols)+
  theme_minimal()


big_plot_any_0_12 <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2",
                         "immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"))%>%
  ggplot(., aes(x=factor(any_0_12), y=cell_freq/100, fill=factor(any_0_12)))+
  geom_boxplot()+
  scale_y_continuous(labels = scales::label_percent())+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free", ncol=4)+
  scale_fill_manual(values = n_infection_cols[4:5])+
  theme_minimal()


big_plot <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2"))%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100, fill=factor(inf_0_12)))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free")+
  # scale_fill_manual(values = n_infection_cols)+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/spot_check_plot.png", big_plot, height = 24, width=16, bg="white")

#needs attention
tfh_pop_ab_purf <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_Th1", "Tfh_Th1_Th17",  "Tfh_Th17", "Tfh_Th2"), Tfh_count>80, timepoint==3)%>%
  group_by(cell_pop, antigen)%>%
  mutate(subset_count=round(cell_freq/100*Tfh_count), id=factor(id))%>%
  mutate(cell_freq_estimate=subset_count/Tfh_count)%>%
  nest() %>%
  mutate(cell_ab_model=map(data, ~MASS::glmmPQL(data=., cell_freq_estimate ~ conc, random = ~1 | id, family=binomial, weights = Tfh_count)))%>%
  mutate(cell_ab_model_summary=map(cell_ab_model, ~summary(.)))%>%
  mutate(cell_ab_model_summary_p=map_dbl(cell_ab_model_summary, ~.$tTable[10]))%>%
  group_by(cell_pop)%>%
  mutate(cell_ab_model_summary_padj= p.adjust(cell_ab_model_summary_p))

tfh_sigs <- filter(cell_pop_ab_purf, cell_ab_model_summary_padj<0.1)
tfh_sigs

#no padj sig, not even a p
abc_pop_ab_purf <- long_true_combo %>%
  filter(cell_pop %in% c("activated_b", "memory_b", "naive_b", "atypical_b"), timepoint==3)%>%
  mutate(subset_count=round(cell_freq/100*mature_b_count), id=factor(id))%>%
  mutate(cell_freq_estimate=subset_count/mature_b_count)%>%
  group_by(cell_pop, antigen)%>%
  nest() %>%
  # mutate(cell_ab_model=map(data, ~MASS::glmmPQL(data=., cell_freq_estimate ~ conc, random = ~1 | id, family=binomial, weights = mature_b_count)))%>%
  mutate(cell_ab_model=map(data, ~MASS::glmmPQL(data=., cell_freq/100 ~ conc, random = ~1 | id, family=binomial)))%>%
  mutate(cell_ab_model_summary=map(cell_ab_model, ~summary(.)))%>%
  mutate(cell_ab_model_summary_p=map_dbl(cell_ab_model_summary, ~.$tTable[10]))%>%
  group_by(cell_pop)%>%
  mutate(cell_ab_model_summary_padj= p.adjust(cell_ab_model_summary_p))

abc_sigs <- filter(abc_pop_ab_purf, cell_ab_model_summary_padj<0.1)
abc_sigs

long_true_combo %>%
  filter(antigen=="GLURP", cell_pop=="Tfh_Th17", timepoint==3)%>%
  ggplot(., aes(x=cell_freq, y=conc))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_minimal()

long_true_combo %>%
  filter(antigen=="EBA181", cell_pop=="Tfh_Th1_Th17", timepoint==3)%>%
  ggplot(., aes(x=cell_freq, y=conc))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_minimal()



cell_pop_incidence_purf <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_Th1", "Tfh_Th1_Th17",  "Tfh_Th17", "Tfh_Th2"), Tfh_count>100, timepoint==3)%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  group_by(cell_pop, Incidence_Type)%>%
  mutate(subset_count=round(cell_freq/100*Tfh_count), id=factor(id))%>%
  mutate(cell_freq_estimate=subset_count/Tfh_count)%>%
  nest() %>%
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., cell_freq_estimate ~ Incidence_Value, random = ~1 | id, family=binomial, weights = Tfh_count)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  group_by(cell_pop)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

sigs <- filter(cell_pop_incidence_purf, cell_incidence_model_summary_padj<0.1)
sigs

hmm <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_Th1", "Tfh_Th1_Th17",  "Tfh_Th17", "Tfh_Th2"),
         Tfh_count>100,
         timepoint==3)%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_6_12", "inf_12_18", "symp_0_12", "symp_6_12", "symp_12_18"))%>%
  ggplot(., aes(x=factor(Incidence_Value), y=cell_freq, fill=factor(Incidence_Value)))+
  geom_boxplot()+
  facet_wrap(cell_pop~Incidence_Type, scales="free")+
  # geom_smooth(method = "lm")+
  scale_fill_manual(values = n_infection_cols)+
  theme_minimal()

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/tfh_incidence_plot.png", hmm, width = 8, height = 16, bg="white")


hmm2 <- long_true_combo %>%
  filter(cell_pop %in% c("Tfh_Th1", "Tfh_Th1_Th17",  "Tfh_Th17", "Tfh_Th2", "activated_b", "memory_b", "naive_b", "atypical_b"),
         Tfh_count>100,
         timepoint==3)%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_6_12", "inf_12_18", "symp_0_12", "symp_6_12", "symp_12_18"))%>%
  ggplot(., aes(x=factor(Incidence_Value), y=cell_freq, fill=factor(Incidence_Value)))+
  geom_boxplot()+
  facet_wrap(cell_pop~Incidence_Type, scales="free")+
  # geom_smooth(method = "lm")+
  scale_fill_manual(values = n_infection_cols)+
  theme_minimal()
ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/tfh_abc_incidence_plot.png", hmm2, width = 11, height = 16, bg="white")


# for(i in 1:nrow(sigs)){
#   
#   plot_data <- long_twelve %>%
#     filter(!is.na(antigen), timepoint==3) %>%
#     filter(recode_cell_pop==x_and_ys[i,1])%>%
#     filter(antigen==x_and_ys[i,2])
#   
#   plot <- ggplot(plot_data, aes(x=freq, y=10^ab_conc))+
#     geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
#     geom_smooth(method="lm")+
#     scale_y_log10()+
#     ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
#     xlab(x_and_ys[i,1])+
#     ylab(x_and_ys[i,2])+
#     theme_minimal()
#   
#   list_of_plots[[i]] <- plot
#   
# }
# 
# big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_cleaned_broomer)/5))
# 
# ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_big_indie_ab_cell_correlation_plot", fdr_cutoff, "raw_p.png", sep="_"), big_plot, height = 6, width=15, bg="white")



flipped_cell_pop_incidence_purf <- long_true_combo %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"), timepoint==3)%>%
  group_by(cell_pop, Incidence_Type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., Incidence_Value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., Incidence_Value ~ cell_freq + (1 | id), family="poisson")))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  group_by(Incidence_Type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

sigs <- flipped_cell_pop_incidence_purf %>%
  dplyr::select(cell_pop, Incidence_Type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

sigs

longer_true_combo <- long_true_combo %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")

list_of_plots <- list()

for(i in 1:nrow(sigs)){
  
  plot_data <- longer_true_combo %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(timepoint==3)%>%
    filter(cell_pop==paste(sigs[i,1]) & Incidence_Type==paste(sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=Incidence_Value, y=cell_freq, fill=factor(Incidence_Value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=incidence_cols)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(sigs[i,1])+
    xlab(sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_plots[[i]] <- plot
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sigs)/5))

ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_big_indie_cell_incidence_plot.png"), big_plot, height = 6, width=6, bg="white")



# solo tfh ####


tfh_clinab <- right_join(combo_data, tfh_combo_batch, by="id")

long_tfh_clinab <- tfh_clinab %>%
  pivot_longer(cols = c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive"), names_to = "cell_pop", values_to = "cell_freq")

tfh_0_12 <- long_tfh_clinab %>%
  filter(!is.na(inf_0_12))%>%
  # filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive",
  #                        "immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"))%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100, fill=factor(inf_0_12)))+
  geom_boxplot()+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free", ncol=4)+
  scale_fill_manual(values = n_infection_cols)+
  theme_minimal()



tfh_incidence_purf <- long_tfh_clinab %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"), timepoint==3)%>%
  group_by(cell_pop, Incidence_Type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., Incidence_Value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  
  mutate(cell_incidence_model=map(data, ~MASS::glm.nb(Incidence_Value ~ cell_freq + id, data=.)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  group_by(Incidence_Type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

tfh_incidence_sigs <- tfh_incidence_purf %>%
  dplyr::select(cell_pop, Incidence_Type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)


tfh_ab_purf <- long_tfh_clinab %>%
  filter(timepoint==3)%>%
  group_by(cell_pop, antigen)%>%
  mutate(subset_count=round(cell_freq/100*Tfh_count), id=factor(id))%>%
  mutate(cell_freq_estimate=subset_count/Tfh_count)%>%
  nest() %>%
  mutate(cell_ab_model=map(data, ~MASS::glmmPQL(data=., cell_freq_estimate ~ conc, random = ~1 | id, family=binomial, weights = Tfh_count)))%>%
  mutate(cell_ab_model_summary=map(cell_ab_model, ~summary(.)))%>%
  mutate(cell_ab_model_summary_p=map_dbl(cell_ab_model_summary, ~.$tTable[10]))%>%
  group_by(cell_pop)%>%
  mutate(cell_ab_model_summary_padj= p.adjust(cell_ab_model_summary_p))
# solo abc ####

abc_combo_batch


abc_clinab <- right_join(combo_data, abc_combo_batch, by="id")

long_abc_clinab <- abc_clinab %>%
  mutate("immature_b_perc"=(1-(mature_b_count/b_count))*100)%>%
  pivot_longer(cols = c("immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"), names_to = "cell_pop", values_to = "cell_freq")



abc_incidence_purf <- long_abc_clinab %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"), timepoint==3)%>%
  group_by(cell_pop, Incidence_Type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  mutate(cell_incidence_model=map(data, ~MASS::glm.nb(Incidence_Value ~ cell_freq + id, data=.)))%>%
  # poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., Incidence_Value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  #poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  
  group_by(Incidence_Type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

sigs <- abc_incidence_purf %>%
  dplyr::select(cell_pop, Incidence_Type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

# cell_pop        Incidence_Type cell_incidence_model_summary_p cell_incidence_model_summary_padj
# <chr>           <chr>                                   <dbl>                             <dbl>
#   1 immature_b_perc inf_0_12                             9.73e-18                          4.86e-17
# 2 immature_b_perc symp_0_12                            5.19e-19                          2.60e-18
# 3 memory_b        inf_0_12                             5.63e-10                          2.25e- 9

abc_0_12 <- long_abc_clinab %>%
  filter(!is.na(inf_0_12))%>%
  # filter(cell_pop %in% c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive",
  #                        "immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"))%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100, fill=factor(inf_0_12)))+
  geom_boxplot()+
  geom_point()+
  scale_y_continuous(labels = scales::label_percent())+
  #geom_smooth(method="lm")+
  #ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  facet_wrap(~cell_pop, scales = "free", ncol=4)+
  scale_fill_manual(values = n_infection_cols)+
  theme_minimal()
