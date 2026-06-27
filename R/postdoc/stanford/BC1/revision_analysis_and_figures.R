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
disease_pal <- rev(colorspace::sequential_hcl("Reds", n=4)[1:3])


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
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "momid", "anyHPfinal", "gestage", "gender", "anymalariapreg", "logpd")))%>%
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

ab_clin <- long_raw_dfff%>%
  # mutate(id=factor(id))%>%
  # filter(timepoint==3)%>%
  inner_join(., infs, by="id")%>%
  dplyr::select(-matches("any"))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")%>%
  pivot_wider(names_from = timepoint, values_from = conc, names_prefix = "t_")%>%
  mutate("increases_6_12"=ifelse(10^t_3>10^t_2, 1, 0))%>%
  pivot_longer(cols = c("t_1", "t_2", "t_3"), names_to = "timepoint", values_to = "conc", names_prefix = "t_")%>%
  filter(id %notin% c(12015, 12028, 12337))

combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])


# recompute infection incidence regressions ####

ab_future_incidence_purf <- ab_clin %>%
  filter(timepoint==1 |
           (timepoint == 2 & incidence_type %notin% c("inf_0_6", "symp_0_6", "inf_0_12", "symp_0_12")) |
           (timepoint == 3 & incidence_type %in% c("symp_12_18", "symp_12_24", "inf_12_18",  "inf_12_24")))%>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  # mutate(nb_cell_incidence_model=map(data, ~lme4::glmer.nb(incidence_value ~ conc + (1|id), data=.)))%>%
  # mutate(nb_cell_incidence_model_summary=map(nb_cell_incidence_model, ~summary(.)))%>%
  # mutate(nb_cell_incidence_model_summary_p=map_dbl(nb_cell_incidence_model_summary, ~coef(.)[11]))%>%
  # 
  # MASS poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ log10(conc), random = ~1 | id, family=poisson)))%>%
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(data=., incidence_value ~ log10(conc))))%>%
  # mutate(cell_incidence_model=map(data, ~glm(data=., incidence_value ~ conc), family="poisson"))%>%
  mutate(cell_incidence_model=map(data, ~glmmTMB::glmmTMB(data=., incidence_value ~ conc, family="nbinom1")))%>%
  # mutate(cell_incidence_model=map(data, ~lm(data=., log10(conc)~incidence_value)))%>%
  # lme4 poisson model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ log10(conc) + (1 | id), family="poisson")))%>%
  
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #MASS poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  #lme4 poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)$cond[8]))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  # mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))

# no significant results, neither with poisson, nor with negative binomial
ab_future_incidence_sigs <- ab_future_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)



ab_past_incidence_purf <- ab_clin %>%
  filter(!(antigen=="MSP1"&timepoint==3),
    (timepoint == 2 & incidence_type %in% c("inf_0_6", "symp_0_6")) |
      (timepoint == 3 & incidence_type %notin% c("symp_12_18", "symp_12_24", "inf_12_18", "inf_12_24"))
  ) %>%
  group_by(antigen, incidence_type, timepoint) %>%
  filter(!duplicated(id)) %>%
  nest()%>%
  # negative binomial model
  # mutate(nb_cell_incidence_model=map(data, ~lme4::glmer.nb(incidence_value ~ conc + (1|id), data=.)))%>%
  # mutate(nb_cell_incidence_model_summary=map(nb_cell_incidence_model, ~summary(.)))%>%
  # mutate(nb_cell_incidence_model_summary_p=map_dbl(nb_cell_incidence_model_summary, ~coef(.)[11]))%>%
  # 
  # MASS poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ log10(conc), random = ~1 | id, family=poisson)))%>%
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(data=., incidence_value ~ log10(conc))))%>%
  mutate(cell_incidence_model=map(data, ~lm(data=., conc ~ incidence_value)))%>%
  # mutate(cell_incidence_model=map(data, ~lm(data=., log10(conc)~incidence_value)))%>%
  # lme4 poisson model
  # mutate(cell_incidence_model=map(data, ~lme4::glmer(data=., incidence_value ~ log10(conc) + (1 | id), family="poisson")))%>%
  
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  #lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint, incidence_type)%>%
  # mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))



ab_past_incidence_sigs <- ab_past_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.15)




ab_past_incidence_sigs%>%
  filter(incidence_type=="inf_0_6", timepoint==2, cell_incidence_model_summary_padj<0.1)%>%
  pull(antigen) -> antigens_exposure_6 

(exposure_6 <- ab_clin %>%
  mutate("timepointf" = recode(timepoint, 
                               "1"="Cord Blood",
                               "2"="6 Months",
                               "3"="12 Months"))%>%
  filter(timepointf=="6 Months", antigen %in% antigens_exposure_6, incidence_type=="inf_0_6")%>%
    mutate(antigen=toupper(antigen))%>%
  ggplot(., aes(x=factor(incidence_value), y=10^conc))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(incidence_value)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Etramp5, GLURP,  Hyp2,  MSP1, SBP1, are significantly increased at six months\nin children who get infected early in life")+
  ylab("Concentration at 6 Months [AU]")+
  xlab("Number of Parasitemic Months in Months 0-6")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 0), limits=10^c(-5, 0))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols))


# ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/exposure_6.png", exposure_6, height=3, width=8, bg="white", dpi=444)
ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/exposure_6_revision.pdf", exposure_6, height=3, width=6, bg="white", dpi=444)


ab_past_incidence_sigs%>%
  filter(incidence_type=="symp_6_12", timepoint==3, cell_incidence_model_summary_padj<0.1)%>%
  pull(antigen) -> antigens_exposure_12 


(exposure_12 <- ab_clin %>%
  mutate("timepointf" = recode(timepoint, 
                               "1"="Cord Blood",
                               "2"="6 Months",
                               "3"="12 Months"))%>%
  filter(timepointf=="12 Months", antigen %in% antigens_exposure_12, incidence_type=="symp_6_12")%>%
    mutate(antigen=toupper(antigen))%>%
  ggplot(., aes(x=factor(incidence_value), y=10^conc))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(incidence_value)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Only MSP1, MSP2 & Etramp5 are significantly increased\nat twelve months in children who get infected in months 6-12")+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Malaria Episodes in Months 6-12")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 0), limits=10^c(-5, NA))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    plot.title = element_text(size=9),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols))


# ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/exposure_12.png", exposure_12, height=3, width=4.8, bg="white", dpi=444)
ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/exposure_12_revision.pdf", exposure_12, height=3, width=4.5, bg="white", dpi=444)


## breadth score ####
breadth_incidence_future_purf <- breadth_and_infs %>%
  filter(!is.na(timepoint))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")%>%
  filter(timepoint==1 |
           (timepoint == 2 & incidence_type %notin% c("inf_0_6", "symp_0_6", "inf_0_12", "symp_0_12")) |
           (timepoint == 3 & incidence_type %in% c("symp_12_18", "symp_12_24", "inf_12_18",  "inf_12_24")))%>%
  group_by(incidence_type, timepoint)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  mutate(cell_incidence_model=map(data, ~glmmTMB::glmmTMB(data=., incidence_value ~ breadth_score, family="nbinom1")))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)$cond[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))


#only exposure, no prediction
breadth_incidence_future_sigs <- breadth_incidence_future_purf %>%
  dplyr::select(incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)



breadth_incidence_past_purf <- breadth_and_infs %>%
  filter(!is.na(timepoint))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")%>%
  filter(
    (timepoint == 2 & incidence_type %in% c("inf_0_6", "symp_0_6")) |
      (timepoint == 3 & incidence_type %notin% c("symp_12_18", "symp_12_24", "inf_12_18", "inf_12_24"))
  ) %>%
  group_by(incidence_type, timepoint) %>%
  filter(!duplicated(id)) %>%
  nest()%>%
  mutate(cell_incidence_model=map(data, ~glm(data=., breadth_score~incidence_value)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))


#only exposure, no prediction
breadth_incidence_past_sigs <- breadth_incidence_past_purf %>%
  dplyr::select(incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)


# recompute cell frequency stats ####
## b cells ####
abc_clin <- infs%>%
  # mutate(id=factor(id))%>%
  inner_join(., abc_combo_batch, by="id")%>%
  mutate("immature_b_perc"=(1-(mature_b_count/b_count))*100)%>%
  dplyr::select(-matches("any"))%>%
  pivot_longer(cols = c("immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")


abc_incidence_purf <- abc_clin %>%
  # there is very little transmission in the b cell cohort so it's hard to evaulate stuff
  # filter(incidence_type %notin% c("symp_0_6", "symp_6_12", "inf_6_12", "symp_0_12"))%>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24", "symp_12_24"))%>%
  # filter(cell_pop!="immature_b_perc")%>%
  group_by(cell_pop, incidence_type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(incidence_value ~ cell_freq + id, data=.)))%>%
  
  # MASS poisson model
  # mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., incidence_value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  # lme4 possoin model
  mutate(cell_incidence_model=map(data, ~glm(data=., cell_freq~incidence_value)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  # MASS poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  # lme4 poisson p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(incidence_type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="fdr"))


abc_incidence_sigs <- abc_incidence_purf %>%
  dplyr::select(cell_pop, incidence_type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)

## t cells ####


tfh_clin <- tfh_combo_batch%>%
  mutate(id=factor(id))%>%
  left_join(., infs, by="id")%>%
  pivot_longer(cols = c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive"), names_to = "cell_pop", values_to = "cell_freq")%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")



tfh_incidence_purf <- tfh_clin %>%
  filter(incidence_type %in% c("inf_0_12", "inf_12_24", "inf_12_18", "symp_0_12", "symp_12_24", "symp_12_18"))%>%
  group_by(cell_pop, incidence_type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # negative binomial model
  # mutate(cell_incidence_model=map(data, ~MASS::glm.nb(incidence_value ~ cell_freq, data=.)))%>%
  
  # poisson model
  mutate(cell_incidence_model=map(data, ~glm(data=., cell_freq~incidence_value)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[11]))%>%
  #poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(incidence_type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="BH"))


tfh_incidence_sigs <- tfh_incidence_purf %>%
  dplyr::select(cell_pop, incidence_type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)


# revised supplementary ####
## maternal chemoprevention ####
(
mom_rx_plot <- combo_data %>%
    mutate("timepointf" = recode(timepoint, 
                                 "1"="Cord Blood",
                                 "2"="6 Months",
                                 "3"="12 Months"))%>%
    filter(!is.na(gestage), timepointf=="Cord Blood")%>%
    ggplot(., aes(x=factor(MomFinalRxx), y=10^conc))+
    geom_point(alpha=0.2, shape=21)+
    geom_boxplot(aes(fill=factor(MomFinalRxx)), outlier.shape = NA)+
    # ggpubr::stat_compare_means(size=2, comparisons = list(c("3 Dose SP", "3 Dose DP"),
    #                                                       c("3 Dose SP", "Monthly DP"), 
    #                                                       c("3 Dose DP", "Monthly DP")))+
    facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=3)+
    xlab("")+
    ylab("Concentration in Cord Blood [AU]")+
    scale_y_log10(labels = scales::label_log())+
    theme_minimal()+
    theme(#panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust=1),
      strip.text = element_text(size=7))+
    scale_fill_viridis_d(option = "G"))

ggsave("/Users/fbach/postdoc/stanford/manuscripts/feli_antibody/author_comments/final_figures/S1_revised.svg", mom_rx_plot, height=10, width = 8)

## maternal malaria ####
(cord_blood_matmal <- combo_data %>%
   mutate(matmal=if_else(anyHPfinal==1, "placental malaria", if_else(anymalariapreg==1, "non-placental malaria", "no malaria")))%>%
   mutate(antigen=toupper(antigen),
          "timepointf" = recode(timepoint, 
                                "1"="Cord Blood",
                                "2"="6 Months",
                                "3"="12 Months"))%>%
   filter(
     timepointf == "Cord Blood",
     !is.na(matmal)) %>%
   mutate(
     matmal = factor(matmal),
     log_conc = log10(conc)
   ) %>%
   ggplot(aes(x = matmal, y = 10^conc, fill = matmal)) +
   geom_point(alpha = 0.2, shape = 21) +
   geom_boxplot(outliers = F) +
   facet_wrap(~antigen, scales="free", labeller = labeller(antigen = label_wrap_gen(width = 6))) +
   scale_y_log10(labels=scales::label_log()) +
   coord_cartesian(ylim = c(-5, NA))+
   ylab("Cord Blood Concentration [AU]")+
   theme_minimal()+
   theme( legend.position="bottom",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.x = element_text(hjust = 0.5, angle = 90, vjust = 0.5),
          axis.title.x = element_blank(),
          strip.text = element_text())+
   scale_fill_manual(values=disease_pal))

ggsave("/Users/fbach/postdoc/stanford/manuscripts/feli_antibody/author_comments/final_figures/S3_revised.svg", cord_blood_matmal, height=10, width=8, bg="white", dpi=444)

## gestational age ####
### main paper ####
gestages_cord <- combo_data %>%
  mutate(antigen=toupper(antigen),
         "timepointf" = recode(timepoint, 
                               "1"="Cord Blood",
                               "2"="6 Months",
                               "3"="12 Months"))%>%
  filter(!is.na(gestage), antigen %in% c("GEXP", "HYP2", "GLURP", "AMA1", "HSP40"), timepointf=="Cord Blood")%>%
  mutate(gestagef=factor(if_else(gestage<28, "<28", 
                                 if_else(gestage>=28 & gestage<32, "28-32", 
                                         if_else(gestage>=32 & gestage<37, "32-37", 
                                                 if_else(gestage>=37, ">37", "whoops")))), levels=c("<28", "28-32", "32-37", ">37")))%>%
  ggplot(., aes(x=gestage, y=10^conc))+
  geom_point(alpha=0.35)+
  geom_smooth(method="lm", aes(color=antigen))+
  # geom_boxplot(aes(fill=gestagef), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "fixed", nrow=1)+
  # ggtitle("Higher Gestational Age Facilitates Transfer of Protective Antibodies")+
  xlab("\n Gestational Age in Weeks")+
  ylab("Cord Blood Concentration [AU]")+
  scale_x_continuous(breaks = seq(32, 42, by=2))+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  ggpubr::stat_cor(method = "spearman", size=2.1)+
  # coord_cartesian(ylim = c(-5, 2))+
  theme_minimal()+
  theme(#panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
    strip.text = element_text())+
  scale_color_manual(values=pc1_cols)

# ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/gestages_cord_scatter.png", gestages_cord, height = 2.7, width=5, bg="white", limitsize = FALSE)
ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/final_version/gestages_cord_scatter.pdf", gestages_cord, height = 2.7, width=5, bg="white", limitsize = FALSE)

### supplementary ####
(all_gestages_corr <- combo_data %>%
    mutate("timepointf" = recode(timepoint, 
                                 "1"="Cord Blood",
                                 "2"="6 Months",
                                 "3"="12 Months"))%>%
    filter(!is.na(gestage), timepointf=="Cord Blood")%>%
    ggplot(., aes(x=gestage, y=10^conc))+
    geom_point(alpha=0.35)+
    geom_smooth(method="lm", aes(color=antigen))+
    facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "fixed", nrow=3)+
    xlab("\n Gestational Age in Weeks")+
    ylab("Cord Blood Concentration [AU]")+
    ggpubr::stat_cor(method = "spearman", size=2, vjust=2, color="red")+
    scale_x_continuous(breaks = seq(32, 42, by=2))+
    scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
    theme_minimal()+
    theme(#panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust=1),
      strip.text = element_text())+
    scale_color_manual(values=pc1_cols))

ggsave("/Users/fbach/postdoc/stanford/manuscripts/feli_antibody/author_comments/final_figures/S4_revised.svg", all_gestages_corr, height=10, width=8, bg="white", dpi=444)

#all no inf ####

all_no_inf_matmal <- combo_data %>%
  mutate("timepointf" = recode(timepoint, 
                               "1"="Cord Blood",
                               "2"="6 Months",
                               "3"="12 Months"))%>%
  mutate(antigen=toupper(antigen))%>%
  filter(inf_0_12==0)%>%
  ggplot(., aes(x=factor(timepointf, levels=c("Cord Blood", "6 Months", "12 Months")), y=10^conc))+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=timepointf), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 3)+
  ylab("Concentration")+
  scale_y_log10(labels=scales::label_log())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1),
    axis.title.x = element_blank(),
    strip.text = element_text(size=7))+
  scale_fill_manual(values=rev(unname(pc1_cols[c(7, 13,19)])))


ggsave("/Users/fbach/postdoc/stanford/manuscripts/feli_antibody/author_comments/final_figures/S5_revised.svg", all_no_inf_matmal, height=10, width=8, bg="white", dpi=444)
