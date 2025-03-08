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
names(pc1_cols) <- sort(c("CSP", "EBA140", "EBA75",  "EBA181", "Etramp4", "Etramp5", "GEXP", "H103", "HSP40", "Hyp2", "MSP2 CH150", "MSP2 DD2", "AMA1", "GLURP","MSP1", "SEA", "Rh2", "Rh4 2",  "Rh5", "SBP1", "TT"))
incidence_cols <- colorspace::sequential_hcl(11, palette = "Purple Yellow")
n_infection_cols <- c("white", colorspace::sequential_hcl(n=5, palette = "Lajolla")[-1])


`%notin%` <- Negate(`%in%`)

fdr_cutoff <- 0.1

dropper <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/kareena_gating/drop_quality_control.csv")
dropper$id <- as.numeric(paste(1, substr(dropper$Sample., 15, 18), sep=""))
ids_to_keep <- subset(dropper, dropper$Flo.Drop=="no", select = "id")
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
  # filter(X!="SD", X!="Mean", !grepl("Control.fcs", tfh_combo_batch$X))%>%
  mutate(file_name=X, .keep = "unused")%>%
  mutate(id=as.numeric(paste(1, substr(.$file_name, 15, 18), sep="")))

tfh_combo_batch <- subset(tfh_combo_batch, tfh_combo_batch$id %in% ids_to_keep$id)


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
  mutate(id=factor(paste(1, substr(.$file_name, 17, 20), sep="")))

abc_combo_batch <- subset(abc_combo_batch, abc_combo_batch$id %in% ids_to_keep$id)

# import ab data ####
bc1 <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

ab_columns <- grep("log", colnames(bc1), value = TRUE)


long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "anyHPfinal")))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))%>%
  mutate(antigen=gsub("log", "", antigen, fixed = TRUE))%>%
  mutate(antigen=gsub("_", " ", antigen, fixed = TRUE))%>%
  mutate(id=factor(id))%>%
  mutate(MomFinalRxx=if_else(MomFinalRx==1, "3 Dose SP",
                             ifelse(MomFinalRx==2, "3 Dose DP",
                                    if_else(MomFinalRx==3, "Monthly DP", "NA")))

  )%>%
  mutate(MomFinalRxx=factor(MomFinalRxx, levels = c("3 Dose SP", "3 Dose DP", "Monthly DP")))%>%
  mutate(anyHPfinalx=if_else(anyHPfinal==1, "Placental Malaria",
                             if_else(anyHPfinal==0, "No Pathology", "Results missing")))

long_raw_dfff %>%
  filter(!is.na(anyHPfinal), timepoint==1)%>%
ggplot(aes(x=factor(anyHPfinal),
                          y=conc,
                          fill=factor(anyHPfinal)
                          ))+
  geom_boxplot()+
  geom_point(shape=21)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()

# ab histopathology ?????
# 
# histo_purf <- long_raw_dfff %>%
#   filter(!is.na(anyHPfinalx), !is.na(conc))%>%
#   group_by(antigen, timepoint)%>%
#   nest() %>%
#   mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., conc ~ timepoint + anyHPfinalx, random = ~1 | id, family=gaussian)))%>%
#   # mutate(cell_incidence_model=map(data, ~lme4::lmer(data=., conc ~ anyHPfinalx + (1 | id), REML=FALSE)))%>%
# 
#   mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
#   mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
#   # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[4]))
#   ungroup()%>%
#   mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))
# 
# 
# histo_purf <- long_raw_dfff %>%
#   filter(!is.na(anyHPfinalx), !is.na(conc))%>%
#   group_by(antigen, timepoint)%>%
#   nest() %>%
#   mutate(cell_incidence_model=map(data, ~lm(conc ~ timepoint + anyHPfinalx, data=.)))%>%
#   # mutate(cell_incidence_model=map(data, ~lme4::lmer(data=., conc ~ anyHPfinalx + (1 | id), REML=FALSE)))%>%
# 
#   mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
#   mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
#   # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[4]))
#   ungroup()%>%
#   group_by(timepoint)%>%
#   mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))
# 
# abc_incidence_sigs <- abc_incidence_purf %>%
#   dplyr::select(cell_pop, incidence_type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
#   filter(cell_incidence_model_summary_padj<0.05)
# 

histopath_wilcox <- long_raw_dfff %>%
  filter(!is.na(anyHPfinalx), !is.na(conc), timepoint==1)%>%
  group_by(antigen)%>%
  pivot_wider(names_from = anyHPfinalx, values_from = conc)%>%
  nest()%>%
  mutate(wilcox=map(data, ~wilcox.test(.$`No Pathology`, .$`Placental Malaria`)))%>%
  mutate(wilcox2=map(data, ~wilcox.test(.$`No Pathology`, .$`Placental Malaria`)))%>%
  mutate(raw_p=map_dbl(wilcox, ~.$p.value))%>%
  ungroup()%>%
  mutate(padj=p.adjust(raw_p, method="fdr"))

sig_histopath <- filter(histopath_wilcox, padj<0.1)

histopathology <- long_raw_dfff %>%
  filter(!is.na(anyHPfinalx), !is.na(conc), timepoint==1, antigen %in% sig_histopath$antigen)%>%
  mutate(anyHPfinalx=case_match(anyHPfinalx,
                                 "No Pathology" ~ "Absent",
                                 "Placental Malaria"~ "Present"))%>%
  ggplot(., aes(x=anyHPfinalx, y=conc))+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=anyHPfinalx))+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 1)+
  # ggtitle("Histopathology")+
  xlab("\nPlacental Malaria")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", hjust=0.5),
        strip.text = element_text(size=13.5, color = "black"))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/histopathology_wilcox.png", histopathology, height = 4, width=8, bg="white", limitsize = FALSE)



# import clinical data ####
# read in new visits database to look at correlations with malaria incidence
# clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_ALL.dta")
clin_data <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/BC-1 childs routine visit database FINAL_REV.dta")

# it's a big table, so let's only include kids for whom we have any antibody measurements
# clin_data <- clin_data %>%
#   filter(id %in% long_raw_dfff$id)

# make a data frame where we add a bunch of columns that contain how many (symptomatic) infections were experienced by the child in the indicated time window

infs <- clin_data %>%
  group_by(id) %>%
  dplyr::select(id, age, anyinfection, sxifinfected, ChildFinalRx) %>%
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
  mutate(ChildFinalRx = case_match(ChildFinalRx, 1~"trimonthly", 2~"monthly"),
         id=factor(id))%>%
  dplyr::select(-anyinfection, -sxifinfected, -age) %>%
  distinct()


# abc incidence stuff ####
abc_clin <- inner_join(infs, abc_combo_batch, by="id")%>%
  mutate("immature_b_perc"=(1-(mature_b_count/b_count))*100)%>%
  dplyr::select(-matches("any"))%>%
  pivot_longer(cols = c("immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")


abc_clin%>%
filter(incidence_type %in% c("inf_0_12", "inf_12_18","inf_12_24", "symp_12_18"))%>%
  arrange(cell_pop)%>%
ggplot(aes(x=factor(incidence_value), y=cell_freq))+
  geom_point(aes(fill=factor(ChildFinalRx)), position = position_dodge(width=0.75), shape=23)+
  geom_boxplot(aes(fill=factor(ChildFinalRx)))+
  facet_wrap(incidence_type~cell_pop, scales="free", ncol=5)+
  # scale_fill_manual(values=incidence_cols)+
  theme_minimal()+
  theme(legend.position = "right")

abc_clin %>%
  group_by(incidence_type, incidence_value)%>%
  filter(!duplicated(id))%>%
  summarise("infs"=n())%>%
  print(n = 62)



abc_incidence_purf <- abc_clin %>%
  # there is very little transmission in the b cell cohort so it's hard to evaulate stuff
  # filter(incidence_type %notin% c("symp_0_6", "symp_6_12", "inf_6_12", "symp_0_12"))%>%
  group_by(cell_pop, incidence_type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(incidence_value ~ cell_freq + id, data=.)))%>%
  
  # MASS poisson model
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  # lme4 possoin model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ cell_freq * ChildFinalRx + (1 | id), family="poisson")))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  # MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  # lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(incidence_type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="fdr"))

abc_incidence_sigs <- abc_incidence_purf %>%
  dplyr::select(cell_pop, incidence_type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

# visualise modelling results

list_of_abc_plots <- list()

for(i in 1:nrow(abc_incidence_sigs)){
  
  plot_data <- abc_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==paste(abc_incidence_sigs[i,1]) & incidence_type==paste(abc_incidence_sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=cell_freq, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=incidence_cols)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(abc_incidence_sigs[i,1])+
    xlab(abc_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_abc_plots[[i]] <- plot
}

sig_abc_plot <- cowplot::plot_grid(plotlist = list_of_abc_plots, nrow = 1)
ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_abc_only_sig_incidence_plot.png"), sig_abc_plot, height = 3, width=4, bg="white")

# tfh incidence stuff ####

tfh_clin <- tfh_combo_batch%>%
  mutate(id=factor(id))%>%
  left_join(., infs, by="id")%>%
  pivot_longer(cols = c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")


tfh_clin%>%
  filter(incidence_type %in% c("symp_12_24"))%>%
  ggplot(aes(x=factor(incidence_value), y=cell_freq,))+
  geom_boxplot(aes(fill=ChildFinalRx))+
  geom_point(aes(fill=ChildFinalRx), shape=21, position=position_dodge(width = 0.75))+
  facet_wrap(incidence_type~cell_pop, scales="free")+
  scale_fill_manual(values=incidence_cols[c(1,5)])+
  theme_minimal()+
  theme(legend.position = "right")

tfh_clin %>%
  group_by(incidence_type, incidence_value)%>%
  filter(!duplicated(id))%>%
  summarise("infs"=n())%>%
  print(n = 62)



tfh_incidence_purf <- tfh_clin %>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24", "inf_12_18", "symp_0_12", "symp_12_24", "symp_12_18"))%>%
  group_by(cell_pop, incidence_type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(incidence_value ~ cell_freq + id, data=.)))%>%
  
  # poisson model
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  #poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  ungroup()%>%
  group_by(incidence_type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))

tfh_incidence_sigs <- tfh_incidence_purf %>%
  dplyr::select(cell_pop, incidence_type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

# visualise modelling results

list_of_tfh_plots <- list()

for(i in 1:nrow(tfh_incidence_sigs)){
  
  plot_data <- tfh_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==paste(tfh_incidence_sigs[i,1]) & incidence_type==paste(tfh_incidence_sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=cell_freq, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=incidence_cols)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(tfh_incidence_sigs[i,1])+
    xlab(tfh_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_tfh_plots[[i]] <- plot
}

sig_tfh_plot <- cowplot::plot_grid(plotlist = list_of_tfh_plots, nrow = 2)

ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_tfh_only_sig_incidence_plot.png"), sig_tfh_plot, height = 6, width=6, bg="white")

# antibody tfh correlations ####

tfh_clin$id <- factor(tfh_clin$id)

thf_ab <- long_raw_dfff%>%
  filter(timepoint==3)%>%
  inner_join(., tfh_clin, by="id")

tfh_ab_broomer <- thf_ab%>%
  group_by(cell_pop, antigen)%>%
  filter(!duplicated(id))%>%
  do(broom::tidy(cor.test(.$cell_freq, .$conc, method="spearman")))%>%
  # group_by(cell_pop)%>%
  ungroup()%>%
  mutate("p_adj"=p.adjust(p.value, method="fdr" ))

tfh_ab_sigs <- tfh_ab_broomer%>%
  filter(p.value<0.05, cell_pop=="Th_memory")



list_of_tfh_ab_plots <- list()

for(i in 1:8){
  
  plot_data <- thf_ab %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==paste(tfh_ab_sigs[i,1]) & antigen==paste(tfh_ab_sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=cell_freq, y=conc, fill=antigen))+
    geom_point(shape=21)+
    geom_smooth(method="lm")+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=pc1_cols)+
    xlab("")+
    ylab("")+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2, )+
    # xlab(tfh_ab_sigs[i,1])+
    ggtitle(tfh_ab_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5))
  list_of_tfh_ab_plots[[i]] <- plot
}

tfh_ab_plot <- cowplot::plot_grid(plotlist = list_of_tfh_ab_plots, nrow = 2)

# ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_tfh_ab_sig_plot.png"), tfh_ab_plot, height = 6, width=4, bg="white")
ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/tfh_ab_raw_p_memory_sig_plot.png"), tfh_ab_plot, height = 5, width=8, bg="white")

# antibody bcell correlations ####
abc_ab <- long_raw_dfff%>%
  filter(timepoint==3)%>%
  inner_join(., abc_clin, by="id")


abc_ab_broomer <- abc_ab%>%
  group_by(cell_pop, antigen)%>%
  filter(!duplicated(id))%>%
  do(broom::tidy(cor.test(.$cell_freq, .$conc, method="spearman")))%>%
  group_by(cell_pop)%>%
  mutate("p_adj"=p.adjust(p.value, method="fdr" ))

abc_ab_sigs <- abc_ab_broomer%>%
  filter(p_adj<0.11)



list_of_abc_ab_plots <- list()

for(i in 1:nrow(abc_ab_sigs)){
  
  plot_data <- abc_ab %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==paste(abc_ab_sigs[i,1]) & antigen==paste(abc_ab_sigs[i,2]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=cell_freq, y=conc, fill=antigen))+
    geom_point(shape=21)+
    geom_smooth(method="lm")+
    # ggtitle(paste(unique(sigs[i,2])))+
    scale_fill_manual(values=pc1_cols)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    xlab(abc_ab_sigs[i,1])+
    ylab(abc_ab_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_abc_ab_plots[[i]] <- plot
}

abc_ab_plot <- cowplot::plot_grid(plotlist = list_of_abc_ab_plots, nrow = 3)

ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/new_gating_abc_ab_sig_plot.png"), abc_ab_plot, height = 8, width=10, bg="white")

# antibody incidence stuff ####

# regular concentration~incidence ##

ab_clin <- long_raw_dfff%>%
  # mutate(id=as.numeric(as.character(id)))%>%
  # filter(timepoint==3)%>%
  inner_join(., infs, by="id")%>%
  dplyr::select(-matches("any"))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")%>%
  pivot_wider(names_from = timepoint, values_from = conc, names_prefix = "t_")%>%
  mutate("increases_6_12"=ifelse(10^t_3>10^t_2, 1, 0))%>%
  pivot_longer(cols = c("t_1", "t_2", "t_3"), names_to = "timepoint", values_to = "conc", names_prefix = "t_")



ab_clin%>%
  # filter(incidence_type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"))%>%
  ggplot(aes(x=factor(incidence_value), y=conc, fill=factor(incidence_value)))+
  geom_point()+
  facet_wrap(incidence_type~antigen, scales="free")+
  scale_fill_manual(values=incidence_cols)+
  theme_minimal()+
  theme(legend.position = "none")

# ab_clin %>%
#   group_by(incidence_type, incidence_value)%>%
#   filter(!duplicated(id))%>%
#   summarise("infs"=n())%>%
#   print(n = 62)





ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  nest() %>%
  # negative binomial model
  # mutate(nb_cell_incidence_model=map(data, ~lme4::glmer.nb(incidence_value ~ conc + (1|id), data=.)))%>%
  # mutate(nb_cell_incidence_model_summary=map(nb_cell_incidence_model, ~summary(.)))%>%
  # mutate(nb_cell_incidence_model_summary_p=map_dbl(nb_cell_incidence_model_summary, ~coef(.)[11]))%>%
  # 
  # MASS poisson model
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ log10(conc), random = ~1 | id, family=poisson)))%>%
  # lme4 poisson model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ log10(conc) + (1 | id), family="poisson")))%>%
  
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  #lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  # mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))



ab_poisson_incidence_sigs <- ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

View(ab_poisson_incidence_sigs)

# visualise modelling results

list_of_ab_poisson_plots <- list()

for(i in 1:nrow(ab_poisson_incidence_sigs)){
  
  plot_data <- ab_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(antigen==paste(ab_poisson_incidence_sigs[i,1]) & incidence_type==paste(ab_poisson_incidence_sigs[i,2]) & timepoint==paste(ab_poisson_incidence_sigs[i,3]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=conc, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    ggtitle(paste(unique(ab_poisson_incidence_sigs[i,3])))+
    scale_fill_manual(values=incidence_cols)+
    # facet_wrap(~timepoint)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(ab_poisson_incidence_sigs[i,1])+
    xlab(ab_poisson_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_ab_poisson_plots[[i]] <- plot
}

sig_ab_poisson_plot <- cowplot::plot_grid(plotlist = list_of_ab_poisson_plots)
ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/ab_only_poisson_sig_incidence_plot.png"), sig_ab_poisson_plot, height = 12, width=12, bg="white")



ab_nb_incidence_sigs <- ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_p, cell_incidence_model_summary_padj, nb_cell_incidence_model_summary_p, nb_cell_incidence_model_summary_padj)%>%
  filter(nb_cell_incidence_model_summary_padj<0.05)


# visualise modelling results

list_of_ab_nb_plots <- list()

for(i in 1:nrow(ab_nb_incidence_sigs)){
  
  plot_data <- ab_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(antigen==paste(ab_nb_incidence_sigs[i,1]) & incidence_type==paste(ab_nb_incidence_sigs[i,2]) & timepoint==paste(ab_nb_incidence_sigs[i,3]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=conc, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    ggtitle(paste(unique(ab_nb_incidence_sigs[i,3])))+
    scale_fill_manual(values=incidence_cols)+
    # facet_wrap(~timepoint)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(ab_nb_incidence_sigs[i,1])+
    xlab(ab_nb_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_ab_nb_plots[[i]] <- plot
}

sig_ab_nb_plot <- cowplot::plot_grid(plotlist = list_of_ab_nb_plots, nrow = 5)
ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/ab_only_NB_sig_incidence_plot.png"), sig_ab_nb_plot, height = 12, width=12, bg="white")


# increasers and non-increasers ~ incidence 
bino_ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type%notin%c("symp_0_12"), antigen != "TT", timepoint=="3")%>%
  nest() %>%
  # poisson model
  # tried flipping model formula to do binomial regression, only two sig, included here, overall bad performance
  mutate(conc_increase_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ factor(increases_6_12), random = ~1 | id, family=poisson)))%>%
  mutate(conc_increase_incidence_model_summary2=map(conc_increase_incidence_model, ~summary(.)))%>%
  mutate(conc_increase_incidence_model_summary2_p=map_dbl(conc_increase_incidence_model_summary2, ~.$tTable[10]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(conc_increase_incidence_model_summary2_padj= p.adjust(conc_increase_incidence_model_summary2_p, method="BH"))


bino_ab_poisson_incidence_sigs <- bino_ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, conc_increase_incidence_model_summary2_p, conc_increase_incidence_model_summary2_padj)%>%
  filter(conc_increase_incidence_model_summary2_padj<0.1, incidence_type %in% c("inf_0_6", "symp_0_6") )

unique(bino_ab_poisson_incidence_sigs$antigen)




list_of_bino_ab_poisson_plots <- list()

for(i in 1:nrow(bino_ab_poisson_incidence_sigs)){
  
  plot_data <- ab_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(antigen==paste(bino_ab_poisson_incidence_sigs[i,1]) & incidence_type==paste(bino_ab_poisson_incidence_sigs[i,2]) & timepoint==paste(bino_ab_poisson_incidence_sigs[i,3]))%>%
    filter(!duplicated(id))%>%
    group_by(antigen, incidence_type, incidence_value)%>%
    filter(!is.na(conc))%>%
    mutate("n_obs"=n())%>%
    mutate("increase_perc"=sum(increases_6_12, na.rm = TRUE)/n_obs)%>%
    filter(!duplicated(increase_perc))
  
  
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=increase_perc, fill=factor(incidence_value)))+
    geom_bar(stat="identity")+
    # geom_point()+
    ggtitle(paste(unique(bino_ab_poisson_incidence_sigs[i,3])))+
    scale_fill_manual(values=incidence_cols)+
    geom_text(aes(label=paste0("frac(",n_obs*increase_perc, ",", n_obs,")")),parse = TRUE, vjust= -0.2, size=3)+
    # facet_wrap(~timepoint)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    scale_y_continuous(limits = c(0,1.25), breaks = seq(0,1,by=0.25), labels = scales::percent_format())+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(paste("increased ", bino_ab_poisson_incidence_sigs[i,1], sep=""))+
    xlab(bino_ab_poisson_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_bino_ab_poisson_plots[[i]] <- plot
}

sig_bino_ab_poisson_plot <- cowplot::plot_grid(plotlist = list_of_bino_ab_poisson_plots, nrow = 5)
ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/bino_ab_only_poisson_sig_incidence_plot.png"), sig_bino_ab_poisson_plot, height = 12, width=12, bg="white")
# memory incidence stuff


tplot <- tfh_clin%>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24"), cell_pop%in%c("Th_memory"))%>%
  ggplot(aes(x=factor(incidence_value), y=cell_freq,fill=factor(incidence_value)))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(incidence_type~cell_pop, scales="free")+
  scale_fill_manual(values=incidence_cols)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())



bplot <- abc_clin%>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24"), cell_pop %in% c("memory_b"))%>%
  ggplot(aes(x=factor(incidence_value), y=cell_freq,fill=factor(incidence_value)))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(incidence_type~cell_pop, scales="free")+
  scale_fill_manual(values=incidence_cols)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())


patch <- bplot / tplot

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/memory_incidence_cont.png", patch, height = 8, width=6, dpi=444, bg="white")



tplot2 <- tfh_clin%>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24"), cell_pop%in%c("Th_memory"))%>%
  mutate(incidence_dich = if_else(incidence_value>0, "some infection", "no infection"))%>%
  ggplot(aes(x=factor(incidence_dich), y=cell_freq,fill=factor(incidence_dich)))+
  geom_boxplot()+
  facet_wrap(incidence_type~cell_pop, scales="free")+
  scale_fill_manual(values=incidence_cols[c(1,5)])+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())



bplot2 <- abc_clin%>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24"), cell_pop %in% c("memory_b"))%>%
  mutate(incidence_dich = if_else(incidence_value>0, "some infection", "no infection"))%>%
  ggplot(aes(x=factor(incidence_dich), y=cell_freq,fill=factor(incidence_dich)))+
  geom_boxplot()+
  facet_wrap(incidence_type~cell_pop, scales="free")+
  scale_fill_manual(values=incidence_cols[c(1,5)])+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

patch2 <- bplot2 / tplot2

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/memory_incidence_dich.png", patch2, height = 8, width=6, dpi=444, bg="white")

big_thing <- full_join(abc_combo_batch, tfh_combo_batch, by = "id")

big_thing %>%
  # pivot_wider(names_from = cell_pop, values_from = cell_freq)%>%
  ggplot(., aes(x=memory_b, y=Th_memory))+
  geom_point()+
  ggpmisc::stat_poly_eq()+
  geom_smooth(method="lm")+
  theme_minimal()
