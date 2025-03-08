# preamble & data generation ####
#palettes
time_palette <- colorspace::sequential_hcl(n=4, "RdPu")[1:3]
pc1_cols <- colorspace::sequential_hcl(23, palette = "Purple Yellow")
names(pc1_cols) <- c("CSP", "EBA140", "EBA75", "EBA181", "Etramp4", "Etramp5", "GEXP", "H103", "HSP40", "Hyp2", "MSP2 CH150", "MSP2 DD2", "AMA1", "GLURP", "MSP1", "SEA", "Rh2", "Rh4 2", "Rh5", "SBP1", "TT")  
incidence_cols <- colorspace::sequential_hcl(11, palette = "Purple Yellow")
n_infection_cols <- c("white", colorspace::sequential_hcl(n=5, palette = "Lajolla")[-1])
gestage_pal <- colorspace::sequential_hcl(n=5, "RdPu")[1:3]
disease_pal <- rev(colorspace::sequential_hcl("Reds", n=4)[1:3])

#libraries 
library(patchwork)
library(tidyr)
library(ggplot2)
library(purrr)
library(MASS)
library(knitr)
library(dplyr)
library(ComplexHeatmap)


`%notin%` <- Negate(`%in%`)

fdr_cutoff <- 0.1

# from looking at the frequencies of observations only these antibodies can be modelled longitudinally; NB some drop off at timepoint 3 so only 1 and 2 can be compared properly
modelable_antigens <- c("Tet Tox", "SBP1", "Rh5", "PfSEA", "PfAMA1", "Hyp2", "HSP40 Ag1", "GST", "GEXP", "CSP GENOVA")


# read in antibody data
bc1 <- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

ab_columns <- grep("log", colnames(bc1), value = TRUE)[-1]


long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "anyHPfinal", "gestage", "gender", "anymalariapreg", "logpd")))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  mutate(conc=10^conc)%>%
  filter(antigen %notin% c("logGST"))%>%
  mutate(antigen=gsub("log", "", antigen,  fixed = TRUE))%>%
  mutate(antigen=gsub("_", " ", antigen, fixed = TRUE))%>%
  mutate(antigen=toupper(antigen))%>%
  mutate(antigen=case_match(antigen,
                        "RH4 2"~"RH4",
                        "MSP2 CH150"~"MSP2-CH150",
                        "MSP2 DD2"~"MSP2-DD2",
                        .default=antigen))%>%
  mutate(MomFinalRx=if_else(MomFinalRx==1, "3 Dose SP",
                            ifelse(MomFinalRx==2, "3 Dose DP",
                                   if_else(MomFinalRx==3, "Monthly DP", "NA")))
         
  )%>%
  mutate(MomFinalRx=factor(MomFinalRx, levels = c("3 Dose SP", "3 Dose DP", "Monthly DP")))%>%
  mutate(anyHPfinalx=if_else(anyHPfinal==1, "Placental Malaria",
                             if_else(anyHPfinal==0, "No Pathology", "Results missing")))%>%
  mutate(matmal=if_else(anyHPfinal==1, "placental malaria", if_else(anymalariapreg==1, "non-placental malaria", "no malaria")))



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
         inf_6_18 = sum(if_else(age>6 & age<18, anyinfection, 0), na.rm = TRUE),
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
combo_data <- combo_data %>%
  filter(id %notin% c(11130, 11084, 11037))%>%
  mutate("timepointf" = recode(timepoint, 
                               "1"="Cord Blood",
                               "2"="6 Months",
                               "3"="12 Months"))%>%
  mutate(timepointf=factor(timepointf, levels=c("Cord Blood", "6 Months", "12 Months")))

raw_combo_data <- cbind(long_raw_dfff, infs[match(long_raw_dfff$id, infs$id),-1])

kids_with_complete_timecourses <- bc1 %>%
  group_by(id)%>%
  summarise("n_time"=n()) %>%
  filter(n_time==3)%>%
  dplyr::select(id)



# add cellular frequencies ##
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


# make antibody cutoff data frame

raw_data <- bc1[,c(1,67:132, 155)]

#shove all antigens in their own column, same for flags and concentrations
raw_conc_data <- raw_data %>%
  pivot_longer(cols = colnames(raw_data)[seq(3, 66, by=3)], values_to = "conc")%>%
  dplyr::select(conc)

raw_flag_data <- raw_data%>%
  pivot_longer(cols = colnames(raw_data)[seq(4, 67, by=3)], values_to = "flag")%>%
  dplyr::select(flag)

# make long data frame combining both conc and flag data
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
  filter(antigen != "GST", antigen!="")

cutoff_df <- long_raw_df %>%
  filter(flag_type=="BelowMinStd", antigen != "GST")%>%
  group_by(antigen)%>%
  summarise("max_below_standard"=max(conc))%>%
  mutate(antigen2=antigen,
         antigen=unique(combo_data$antigen))




# wide_combo_12months <- combo_data |> 
#   dplyr::select(id, timepoint, anyHPfinal, antigen, conc)|>
#   filter(timepoint==3)|>
#   pivot_wider(names_from = antigen, values_from = conc)
# 
# twelve_cor_combo <- inner_join(wide_combo_12months, combo_cells_data, by="id")

# twelve_cor_combo <- very_long_combo_combo %>%
#   filter(timepoint==3)%>%
#   select(id, timepoint, anyHPfinal, antigen, conc, cell_pop, cell_freq)%>%
#   pivot_wider(names_from = cell_pop, values_from = cell_freq)%>%
#   pivot_wider(names_from = antigen, values_from = conc)


# these have padj < 0.1;  > raw_p range 0.009 - 0.015
placental_antigens <- c("GEXP", "TT", "AMA1", "HYP2", "SEA")

# fig 0. venn diagram of # of individuals with sample type ####

library(ggvenn)

feli_data <- list("Malaria Serology"=combo_data$id,
                  "T Cell Cytometry"=tfh_combo_batch$id,
                  "B Cell Cytometry"=abc_combo_batch$id)
# make data a list (of assyas) of lists (of individuals) 
venn_diagram <- ggvenn(
  feli_data, 
  show_stats = "c",
  text_color = "black",
  text_size = 6,
  fill_color = c("#CD534CFF", "#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
ggsave("~/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/sample_set_diagram.png", venn_diagram, width=5, height=5, bg='white')

    # fig 1. gestational age, maternal malaria, 0 infection kids ####
# Ab responses stratified by gestational age,
# Ab responses in cord blood stratified by pregnancy malaria outcome 
# Ab responses of no-inf kids and effects later on

gestages_cord <- combo_data %>%
  filter(!is.na(gestage), !is.na(matmal), antigen %in% c("GEXP", "HYP2", "GLURP", "AMA1", "HSP40"), timepointf=="Cord Blood")%>%
  mutate(gestagef=factor(if_else(gestage<28, "<28", 
                                 if_else(gestage>=28 & gestage<32, "28-32", 
                                         if_else(gestage>=32 & gestage<37, "32-37", 
                                                 if_else(gestage>=37, ">37", "whoops")))), levels=c("<28", "28-32", "32-37", ">37")))%>%
  ggplot(., aes(x=gestage, y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c("GEXP", "HYP2", "GLURP", "AMA1", "HSP40")), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.35)+
  geom_smooth(method="lm", aes(color=antigen))+
  # geom_boxplot(aes(fill=gestagef), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "fixed", nrow=1)+
  # ggtitle("Higher Gestational Age Facilitates Transfer of Protective Antibodies")+
  xlab("\n Gestational Age in Weeks")+
  ylab("Cord Blood Concentration [AU]")+
  ggpubr::stat_cor(method = "spearman", size=2.3)+
  scale_x_continuous(breaks = seq(32, 42, by=2))+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(#panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_color_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/gestages_cord_scatter.png", gestages_cord, height = 2.7, width=5, bg="white", limitsize = FALSE)


cord_blood_matmal <- combo_data %>%
  filter(timepointf=="Cord Blood", !is.na(matmal), antigen %in% placental_antigens, is.finite(conc))%>%
  ggplot(., aes(x=matmal, y=conc))+
  # scale_x_discrete(labels = scales::label_wrap(11))+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  # annotation_logticks(sides = "l", size = 0.1)+
  geom_hline(data=filter(cutoff_df, antigen %in% placental_antigens), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=matmal), outlier.shape = NA,)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "fixed", nrow=1)+
  # ggtitle("Maternal Pf Infection Negatively Impacts Antibody Transfer")+
  ylab("Cord Blood Concentration [AU]")+
  theme_minimal()+
    theme(
        legend.position="bottom",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(hjust = 0.5, angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.text = element_text())+
  scale_fill_manual(values=disease_pal)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/cord_blood_matmal.png", cord_blood_matmal, height=2.7, width=5, bg="white", dpi=444)



# fig 2. without infection antibodies quickly wane below linear range of assay ####
no_inf_matmal <- combo_data %>%
  filter(inf_0_12==0, antigen %in% c(placental_antigens), is.finite(conc))%>%
  ggplot(., aes(x=timepointf, y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c(placental_antigens)), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_line(aes(group=id), alpha=0.1)+
  geom_boxplot(aes(fill=timepointf), outlier.shape = NA)+
  # stat_summary(fun = median, fun.min = median, fun.max = median,
  #              geom = "crossbar", width = 0.5, color="darkred")+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("In The Absence of Infection, Maternal Antibodies Wane Below Limit of Quantification\nWithin 6 Months")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  ylab("Concentration [AU]")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1),
    axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=rev(unname(pc1_cols[c(7, 13,19)])))

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/no_inf_six_twelve.png", no_inf_matmal, height=2.7, width=5, bg="white", dpi=444)
{# yes_inf_matmal <- combo_data %>%
#   filter(inf_0_12!=0, antigen %in% placental_antigens, is.finite(conc))%>%
#   ggplot(., aes(x=timepointf, y=conc))+
#   geom_hline(data=filter(cutoff_df, antigen %in% placental_antigens), aes(yintercept = max_below_standard), linetype="dashed")+
#   # geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
#   geom_boxplot(aes(fill=antigen), outlier.shape = NA)+
#   facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
#   # ggtitle("In The Absence of Infection, Maternal Antibodies Wane Below Limit of Quantification\nWithin 6 Months")+
#   scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(
#     legend.position="none",
#     legend.title = element_blank(),
#     axis.text.x = element_text(angle=90, hjust=1),
#     axis.title.x = element_blank(),
#     strip.text = element_text())+
#   scale_fill_manual(values=pc1_cols)
  }


# fig 3. antibody markers of exposure and/or protection ####
protection <- combo_data %>%
  filter(timepointf=="Cord Blood", antigen %in% c("HSP40", "SEA"))%>%
  filter(antigen %in% c("HSP40", "SEA"))%>%
  ggplot(., aes(x=factor(symp_0_12), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c("HSP40", "SEA")), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_0_12)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Cord Blood Antibodies HSP40 & SEA are Significantly Higher in Infants who Remain\nUninfected in the First Six Months of Life")+
  ylab("Cord Blood Concentration [AU]")+
  xlab("Number of Malaria Episodes\nin the First Year of Life")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.title.x = element_text(hjust=0),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/protection_cord.png", protection, height=3, width=3, bg="white", dpi=444)


 antigens_exposure_6 <- c("ETRAMP5", "GLURP",  "HYP2",  "MSP1", "SBP1")

exposure_6 <- combo_data %>%
  filter(timepointf=="6 Months", antigen %in% antigens_exposure_6)%>%
  ggplot(., aes(x=factor(inf_0_6), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_exposure_6), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_0_6)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Etramp5, GLURP,  Hyp2,  MSP1, SBP1, are significantly increased at six months\nin children who get infected early in life")+
  ylab("Concentration at 6 Months [AU]")+
  xlab("Number of Infections in Months 0-6")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 0), limits=10^c(-5, 0))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/exposure_6.png", exposure_6, height=3, width=8, bg="white", dpi=444)


antigens_exposure_12 <- c("MSP1", "MSP2-CH150", "ETRAMP5")

exposure_12 <- combo_data %>%
  filter(timepointf=="12 Months", antigen %in% antigens_exposure_12)%>%
  ggplot(., aes(x=factor(inf_6_12), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_exposure_12), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_6_12)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Only MSP1, MSP2 & Etramp5 are significantly increased\nat twelve months in children who get infected in months 6-12")+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Infections in Months 6-12")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 0), limits=10^c(-5, 0))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    plot.title = element_text(size=9),
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/exposure_12.png", exposure_12, height=3, width=4.8, bg="white", dpi=444)

# fig 4. early life malaria causes ab defects ####
# antigens_defect_12 <- c("Etramp5", "H103", "Hyp2", "GLURP", "MSP1", "Rh4 2", "SBP1")
antigens_defect_symp_12 <- c("CSP", "EBA75", "H103", "AMA1", "GLURP")

defect_12_symp <- combo_data %>%
  filter(timepointf=="12 Months", antigen %in% antigens_defect_symp_12)%>%
  # mutate(symp_0_6=case_when(symp_0_6=="0"~"none",
  #                               symp_0_6=="1"~"some"))%>%
  ggplot(., aes(x=factor(symp_0_6), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_defect_symp_12), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_0_6)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Etramp5, H103, Hyp2, GLURP, MSP1, Rh4 2, SBP1 are lower at 12 months\nin children who were infected in the first six months of life")+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Malaria Episodes in Months 0-6")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-8, 0))+
  theme_minimal()+
  theme(
      legend.position="none",
      legend.title = element_blank(),
      # axis.text.x = element_text(angle = 90, hjust=1),
      strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/defect_12_symp.png", defect_12_symp, height=3, width=8, bg="white", dpi=444)



antigens_defect_inf_12 <- c("RH4", "EBA140", "H103")

defect_12_inf <- combo_data %>%
  filter(timepointf=="12 Months", antigen %in% antigens_defect_inf_12)%>%
  # mutate(symp_0_6=case_when(symp_0_6=="0"~"none",
  #                               symp_0_6=="1"~"some"))%>%
  ggplot(., aes(x=factor(inf_0_6), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_defect_inf_12), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_0_6)), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  # ggtitle("Etramp5, H103, Hyp2, GLURP, MSP1, Rh4 2, SBP1 are lower at 12 months\nin children who were infected in the first six months of life")+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Infections in Months 0-6")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-8, 0))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/defect_12_inf.png", defect_12_inf, height=3, width=4.8, bg="white", dpi=444)





# fig 5 impact of early life infection on b & t cell subsets ####

tfh_clinab <- inner_join(combo_data, tfh_combo_batch, by="id")

long_tfh_clinab <- tfh_clinab %>%
  pivot_longer(cols = c("Tfh_perc", "Tfh_Th1", "Tfh_Th1_Th17", "Tfh_Th17", "Tfh_Th2","Th_memory", "Th_naive"), names_to = "cell_pop", values_to = "cell_freq")%>%
  mutate(cell_popf=case_when(cell_pop=="Tfh_perc"~"Tfh % of CD4",
                             cell_pop=="Tfh_Th1"~"Th1-like of Tfh",
                             cell_pop=="Tfh_Th2"~"Th2-like of Tfh",
                             cell_pop=="Tfh_Th17"~"Th17-like of Tfh",
                             cell_pop=="Tfh_Th1_Th17"~"Th1/Th17-like of Tfh",
                             cell_pop=="Th_naive"~"naive T Cells",
                             cell_pop=="Th_memory"~"memory T Cells"))

tfh_incidence_purf <- long_tfh_clinab %>%
  select(-starts_with("any_"), -c(symp_6_12, inf_6_12))%>%
  filter(timepoint==3)%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  # filter(Incidence_Type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"), timepoint==3)%>%
  group_by(cell_pop, Incidence_Type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # mutate(cell_incidence_model=map(data, ~glm(data=., Incidence_Value ~ cell_freq, family="poisson")))%>%
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., Incidence_Value~cell_freq , random = ~1 | id,family=poisson)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  # group_by(Incidence_Type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="fdr"))

tfh_incidence_sigs <- tfh_incidence_purf %>%
  dplyr::select(cell_pop, Incidence_Type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)





# tfh_0_6 <- long_tfh_clinab %>%
#   filter(!duplicated(id), !is.na(inf_0_6))%>%
#   ggplot(., aes(x=factor(any_inf_0_6), y=cell_freq/100, fill=factor(any_inf_0_6)))+
#   geom_boxplot()+
#   geom_point()+
#   scale_y_continuous(labels = scales::label_percent())+
#   facet_wrap(~cell_pop, scales = "free", ncol=4)+
#   scale_fill_manual(values = n_infection_cols)+
#   theme_minimal()

big_n_infection_cols <- c(n_infection_cols[1:4],  colorspace::sequential_hcl(n=20, palette = "Lajolla")[17:20])

tfh_0_12 <- long_tfh_clinab %>%
  filter(!is.na(symp_0_12), cell_pop %in% c("Tfh_Th1_Th17", "Tfh_Th1", "Th_memory","Th_naive"))%>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_0_12), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_0_12)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  # ggtitle("early life malaria is associated with decreased memory differentiation in T Cells \nand a slight shift in Tfh differentiation")+
  ylab("% of Parent Population at 12 Months")+
  xlab("Number of Malaria Epiosdes in The First Year of Life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=big_n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/t_cell_incidence_012.png", tfh_0_12, height=3, width=8, bg="white", dpi=444)




#b cell stuff
abc_clinab <- right_join(combo_data, abc_combo_batch, by="id")

long_abc_clinab <- abc_clinab %>%
  mutate("immature_b_perc"=(1-(mature_b_count/b_count))*100)%>%
  pivot_longer(cols = c("immature_b_perc", "activated_b", "memory_b", "naive_b", "atypical_b"), names_to = "cell_pop", values_to = "cell_freq")%>%
  mutate(cell_popf=case_when(cell_pop=="activated_b"~"activated B Cells",
                             cell_pop=="atypical_b"~"atypical B Cells",
                             cell_pop=="memory_b"~"memory B Cells",
                             cell_pop=="naive_b"~"naive B Cells",
                             cell_pop=="immature_b_perc"~"CD10+ B Cells"))




abc_incidence_purf <- long_abc_clinab %>%
  filter(timepoint==3)%>%
  select(-starts_with("any_"), -c(symp_6_12, inf_6_18))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(Incidence_Type %in% c("inf_0_12", "inf_12_24", "symp_0_12", "symp_12_24"))%>%
  group_by(cell_pop, Incidence_Type)%>%
  filter(!duplicated(id))%>%
  nest() %>%
  # glm
  # mutate(cell_incidence_model=map(data, ~glm(Incidence_Value ~ cell_freq, data=.)))%>%
  # poisson model
  mutate(cell_incidence_model=map(data, ~MASS::glmmPQL(data=., Incidence_Value ~ cell_freq, random = ~1 | id, family=poisson)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  #negative binomial p
  # mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  #poisson p
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~.$tTable[10]))%>%
  ungroup()%>%
  # group_by(Incidence_Type)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p, method="fdr"))

abc_incidence_sigs <- abc_incidence_purf %>%
  dplyr::select(cell_pop, Incidence_Type, cell_incidence_model_summary_p, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)



abc_0_12 <- long_abc_clinab %>%
  filter(!is.na(symp_0_12), cell_pop %in% c("memory_b", "naive_b", "immature_b_perc", "atypical_b"))%>%
  # filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_0_12), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(symp_0_12)), outlier.shape = NA)+
  facet_wrap(~factor(cell_popf, levels=c("atypical B Cells", "activated B Cells", "naive B Cells", "memory B Cells", "CD10+ B Cells")), labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  # ggtitle("early life malaria is associated with increased\nmemory differentiation in B cells")+
  ylab("% of Parent Population at 12 Months")+
  xlab("Number of Malaria Episodes in The First Year of Life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/b_cell_incidence_0_12.png", abc_0_12, height=3, width=8, bg="white", dpi=444)


abc_0_6 <- long_abc_clinab %>%
  filter(!is.na(inf__12), cell_pop %in% c("atypical_b"))%>%
  # filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(inf_6_12), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_6_12)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  # ggtitle("early life malaria is associated with increased\nmemory differentiation in B cells")+
  ylab("% of Parent Population at 12 Months")+
  xlab("Number of Infections in The First Year of Life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/b_cell_incidence_0_12.png", abc_0_12, height=3, width=4, bg="white", dpi=444)





 # fig 6 antibody ~ cell interactions ####
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
  mutate(cell_ab_model_summary_padj= p.adjust(cell_ab_model_summary_p, method="fdr"))

tfh_ab_sigs <- tfh_ab_purf%>%
  filter(cell_ab_model_summary_padj<0.1)

tfh_ab_broomer <- long_tfh_clinab%>%
  group_by(cell_pop, antigen)%>%
  filter(!duplicated(id))%>%
  do(broom::tidy(cor.test(.$cell_freq, .$conc, method="spearman")))%>%
  ungroup()%>%
  mutate("p_adj"=p.adjust(p.value, method="fdr"))


# needs work


abc_ab_purf <- long_abc_clinab %>%
  filter(timepoint==3, antigen %notin%c("EBA75"))%>%
  group_by(cell_pop, antigen)%>%
  mutate(subset_count=round(cell_freq/100*b_count), id=factor(id))%>%
  mutate(cell_freq_estimate=subset_count/b_count)%>%
  nest() %>%
  mutate(cell_ab_model=map(data, ~MASS::glmmPQL(data=., cell_freq_estimate ~ conc, random = ~1 | id, family=binomial, weights = b_count)))%>%
  mutate(cell_ab_model_summary=map(cell_ab_model, ~summary(.)))%>%
  mutate(cell_ab_model_summary_p=map_dbl(cell_ab_model_summary, ~.$tTable[10]))%>%
  group_by(cell_pop)%>%
  mutate(cell_ab_model_summary_padj= p.adjust(cell_ab_model_summary_p, method="fdr"))

abc_ab_sigs <- abc_ab_purf%>%
  filter(cell_ab_model_summary_padj<0.1)

all_ab_broomer <- long_true_combo%>%
  group_by(cell_pop, antigen)%>%
  do(broom::tidy(cor.test(.$cell_freq, .$conc, method="spearman")))%>%
  ungroup()%>%
  mutate("p_adj"=p.adjust(p.value, method="fdr" ))

# supp. fig. 1 ALL DATA gestational age, maternal malaria, 0 infection kids ####

gestage_purf <- combo_data %>%
  filter(!is.na(gestage), timepointf=="Cord Blood")%>%
  mutate(gestagef=factor(if_else(gestage<28, "<28", 
                                 if_else(gestage>=28 & gestage<32, "28-32", 
                                         if_else(gestage>=32 & gestage<37, "32-37", 
                                                 if_else(gestage>=37, ">37", "whoops")))), levels=c("<28", "28-32", "32-37", ">37")))%>%
  # filter(gestagef %in%c(">37","32-32", "28-32"))%>%
  group_by(antigen)%>%
  nest()%>%
  mutate(model=map(data, ~lm(data=., log10(conc)~gestage)))%>%
  mutate(model_summary=map(model, ~summary(.)))%>%
  mutate(model_summary_p=map_dbl(model_summary, ~.$coefficients[8]))%>%
  ungroup()%>%
  mutate(model_summary_padj= p.adjust(model_summary_p, method="fdr"))

  
filter(gestage_purf, model_summary_padj<0.1)
  

all_gestages_cord <- combo_data %>%
  filter(!is.na(gestage), timepointf=="Cord Blood")%>%
  mutate(gestagef=factor(if_else(gestage<28, "<28", 
                                 if_else(gestage>=28 & gestage<32, "28-32", 
                                         if_else(gestage>=32 & gestage<37, "32-37", 
                                                 if_else(gestage>=37, ">37", "whoops")))), levels=c("<28", "28-32", "32-37", ">37")))%>%
  ggplot(., aes(x=gestage, y=log10(conc)))+
  geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  # geom_boxplot(aes(fill=gestagef, group=gestagef), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=3)+
  ggtitle("Higher Gestational Age Facilitates Transfer of Antibodies")+
  xlab("\n Gestational Age in Weeks")+
  ylab("Concentration")+
  theme_minimal()+
  theme(#panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=gestage_pal)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_gestages_cord.png", all_gestages_cord, height = 7.5, width=8, bg="white", limitsize = FALSE)



all_cord_blood_matmal <- combo_data %>%
  filter(timepointf=="Cord Blood", !is.na(matmal))%>%
  ggplot(., aes(x=matmal, y=log10(conc)))+
  geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=matmal), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=3)+
  ggtitle("Maternal Pf Infection Negatively Impacts Antibody Transfer")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=disease_pal)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_cord_blood_matmal.png", all_cord_blood_matmal, height=7.5, width=8, bg="white", dpi=444)



all_no_inf_matmal <- combo_data %>%
  filter(inf_0_12==0)%>%
  ggplot(., aes(x=timepointf, y=log10(conc)))+
  geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=timepointf), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 3)+
  ggtitle("Maternal Antibodies Wane Within 6 Months")+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1),
    axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=rev(unname(pc1_cols[c(7, 13,19)])))


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_no_inf_six_twelve.png", all_no_inf_matmal, height=7.5, width=8, bg="white", dpi=444)


# supp fig. 2 antibody vs. exposures


all_zero_6 <- combo_data %>%
  filter(!is.na(inf_0_6), timepointf=="6 Months")%>%
  ggplot(., aes(x=factor(inf_0_6), y=log10(conc)))+
  geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_0_6)), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 3)+
  # ggtitle("Maternal Antibodies Wane Within 6 Months, Removing the Impact Maternal Malaria")+
  ylab("Concentration at 6 Months [AU]")+
  xlab("Number of Infections Months 0-6")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle=90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_inf06.png", all_zero_6, height=7.5, width=8, bg="white", dpi=444)



all_6_12 <- combo_data %>%
  filter(!is.na(inf_6_12), timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(inf_6_12), y=log10(conc)))+
  geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_6_12)), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 3)+
  # ggtitle("Maternal Antibodies Wane Within 6 Months, Removing the Impact Maternal Malaria")+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Infections Months 6-12")+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle=90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_inf612.png", all_6_12, height=7.5, width=8, bg="white", dpi=444)




# c_all <- combo_data %>%
#   filter(timepointf == "12 Months", inf_6_12==0)%>%
#   mutate(any_inf_0_6=case_when(any_inf_0_6=="0"~"none",
#                                any_inf_0_6=="1"~"one or more"))%>%
#   ggplot(., aes(x=factor(any_inf_0_6), y=conc))+
#   geom_hline(data=filter(cutoff_df), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
#   geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
#   geom_boxplot(aes(fill=factor(any_inf_0_6, levels = c("none", "one or more"))), outlier.shape = NA)+
#   facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
#   # ggtitle("parasitemia in The First Six Months of a Defect in\nAntimalarial Antibodies in Infants at 12 Months of Age")+
#   xlab("Number of Infections Months 0 to 6")+
#   ylab("Concentration")+
#   theme_minimal()+
#   theme(
#     legend.position="none",
#     legend.title = element_blank(),
#     # axis.text.x = element_text(angle=90, hjust=1),
#     strip.text = element_text())+
#   scale_fill_manual(values=n_infection_cols)
# 
# ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/c_all.png", c_all, height=4, width=16, bg="white", dpi=444)

# supplementary Figure X: early life incidence predicts later life incidence ####
early_late_abc_incidence <- long_abc_clinab %>%
  filter(!is.na(inf_0_12), !is.na(inf_12_24))%>%
  filter(!duplicated(id))%>%
  ggplot(., aes(x=inf_0_12, y=inf_12_18))+
  geom_point(shape=21, position = position_jitter(width = 0.1, height = 0.1))+
  ggtitle("incidence in B cell data")+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman", label.y = 3.5, na.rm = TRUE, size=2)+
  ylab("Number of Infections Months 12-18")+
  xlab("Number of Infections Months 0-12")+
  # scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  coord_cartesian(xlim = c(0,4), ylim=c(0,4))+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

early_late_tfh_incidence <- long_tfh_clinab %>%
  filter(!is.na(inf_0_12), !is.na(inf_12_24))%>%
  filter(!duplicated(id))%>%
  ggplot(., aes(x=inf_0_12, y=inf_12_18))+
  geom_point(shape=21, position = position_jitter(width = 0.1, height = 0.1))+
  ggtitle("incidence in T cell data")+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman", label.y = 3.5, na.rm = TRUE, size=2)+
  ylab("Number of Infections Months 12-18")+
  xlab("Number of Infections Months 0-12")+
  # scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  coord_cartesian(xlim = c(0,4), ylim=c(0,4))+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


early_late_incidence_plot <- early_late_tfh_incidence + early_late_abc_incidence

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/early_late_incidence_plot.png", early_late_incidence_plot, height=3, width=6, bg="white", dpi=444)


long_abc_clinab %>%
  filter(!is.na(inf_0_12), !is.na(inf_12_18))%>%
  filter(!duplicated(id))%>%
  do(broom::tidy(cor.test(.$inf_0_12, .$inf_12_24, method="pearson")))


long_tfh_clinab %>%
  filter(!is.na(inf_0_12), !is.na(inf_12_24))%>%
  filter(!duplicated(id))%>%
  ggplot(., aes(x=inf_0_12, y=inf_12_18))+
  geom_point(shape=21, position = position_jitter(width = 0.1, height = 0.1))+
  ggtitle("incidence in T cell data")+
  geom_smooth(method="lm")+
  # ylab("% of parent population")+
  # xlab("Number of Infections in the first year of life")+
  # scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  coord_cartesian(xlim = c(0,4), ylim=c(0,4))+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)


long_tfh_clinab %>%
  filter(!is.na(inf_0_12), !is.na(inf_12_18))%>%
  filter(!duplicated(id))%>%
  do(broom::tidy(cor.test(.$inf_0_12, .$inf_12_24, method="pearson")))

# supplementary figure X: T memory association looks reasonable even when controlling for early life exposure,
# supplemenraey Figure X: no protective effect of HSP40, SEA antibodies measured at 6 or 12 months, looking 12 months into the future ####

n_infection_cols_supp <- c("white", colorspace::sequential_hcl(n=8, palette = "Lajolla")[-1], "black")

six_prot <- combo_data %>%
  filter(antigen %in% c("SEA", "HSP40"), is.finite(conc), timepoint==2)%>%
  ggplot(., aes(x=factor(inf_6_18), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c("SEA", "HSP40")), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_6_18)), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 1)+
  ylab("Concentration at 6 Months [AU]")+
  xlab("Number of Infections Months 6-18")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle=90, hjust=1),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols_supp)

twelve_prot <- combo_data %>%
  filter(antigen %in% c("SEA", "HSP40"), is.finite(conc), timepoint==2)%>%
  ggplot(., aes(x=factor(inf_12_24), y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% c("SEA", "HSP40")), aes(yintercept = max_below_standard), linetype="dashed")+
  geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=factor(inf_12_24)), outlier.shape = NA)+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow = 1)+
  ylab("Concentration at 12 Months [AU]")+
  xlab("Number of Infections Months 12-24")+
  scale_y_log10(labels=scales::label_log(), breaks=10^seq(-5, 2), limits=10^c(-5, 2.5))+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle=90, hjust=1),
    # axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols_supp)

no_protection_supp_figure <- six_prot / twelve_prot

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/supp_no_protection_figure.png", no_protection_supp_figure, height=5, width=5, bg="white", dpi=444)

# supplementary Figure X: probability of increasing Ab concentration from month 6 to 12 based on early life infections

ab_clin <- long_raw_dfff%>%
  # filter(timepoint==3)%>%
  inner_join(., infs, by="id")%>%
  dplyr::select(-matches("any"))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "incidence_type", values_to = "incidence_value")%>%
  pivot_wider(names_from = timepoint, values_from = conc, names_prefix = "t_")%>%
  mutate("increases_6_12"=ifelse(10^t_3>10^t_2, 1, 0))%>%
  pivot_longer(cols = c("t_1", "t_2", "t_3"), names_to = "timepoint", values_to = "conc", names_prefix = "t_")


cont_defect_data <- ab_clin %>%
  # filter(!is.na(antigen), timepoint==3) %>%
  filter(incidence_type=="inf_0_6", timepoint==3)%>%
  group_by(antigen, incidence_type, incidence_value)%>%
  filter(!is.na(conc))%>%
  mutate("n_obs"=n())%>%
  mutate("increase_perc"=sum(increases_6_12, na.rm = TRUE)/n_obs)%>%
  filter(!duplicated(id), !duplicated(increase_perc))

# unique(bino_ab_poisson_incidence_sigs$antigen)
# [1] "CSP"     "EBA140"  "EBA181"  "Etramp5" "GEXP"    "H103"    "HSP40"   "Hyp2"    "GLURP"   "MSP1"   
# [11] "SEA"     "Rh2"     "Rh5"     "SBP1"  
defect_abs_inf <- c("CSP", "EBA140", "Etramp5", "GEXP", "HSP40",  "HYP2", "GLURP", "MSP1", "SEA", "RH2", "SBP1")

cont_defect_plot_inf <- cont_defect_data %>%
  filter(antigen %in%defect_abs_inf)%>%
  ggplot(., aes(x=factor(incidence_value), y=increase_perc, fill=factor(incidence_value)))+
  geom_bar(stat="identity", color="black")+
  geom_text(aes(label=paste0("frac(",n_obs*increase_perc, ",", n_obs,")")),parse = TRUE, vjust= -0.2, size=3)+
  facet_wrap(~antigen, nrow=2)+
  scale_y_continuous(limits = c(0,1.25), breaks = seq(0,1,by=0.25), labels = scales::percent_format())+
  scale_fill_manual(values=n_infection_cols)+
  # ggtitle("Children who were infected in the first six months of life are less likely to increase\nantibody responses from 6 to 12 months")+
  ylab("Fraction of Children Who Increase Concentration From Month 6 to 12")+
  xlab("Infections in Months 0-6")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/cont_defect_plot_inf.png", cont_defect_plot_inf, height=6, width=8, bg="white", dpi=444)



defect_abs_symp <- c("EBA181", "Hyp2", "GLURP", "MSP1", "Rh2",  "Rh5")

cont_defect_plot_symp <- cont_defect_data %>%
  filter(antigen %in%defect_abs_symp)%>%
  ggplot(., aes(x=factor(incidence_value), y=increase_perc, fill=factor(incidence_value)))+
  geom_bar(stat="identity", color="black")+
  geom_text(aes(label=paste0("frac(",n_obs*increase_perc, ",", n_obs,")")),parse = TRUE, vjust= -0.2, size=3)+
  facet_wrap(~antigen, nrow=1)+
  scale_y_continuous(limits = c(0,1.25), breaks = seq(0,1,by=0.25), labels = scales::percent_format())+
  scale_fill_manual(values=n_infection_cols)+
  ylab("Fraction of Children Who Increase Concentration From Month 6 to 12")+
  xlab("Malaria Episodes in Months 0-6")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/cont_defect_plot_symp.png", cont_defect_plot_symp, height=3, width=8, bg="white", dpi=444)


# binary outcome rather than continuous


binary_defect_data <- ab_clin %>%
  # filter(!is.na(antigen), timepoint==3) %>%
  filter(incidence_type=="inf_0_6", timepoint==3)%>%
  mutate(any_inf_0_6=if_else(incidence_value>0, "some", "none"))%>%
  group_by(antigen, any_inf_0_6)%>%
  filter(!is.na(conc))%>%
  mutate("n_obs"=n())%>%
  mutate("increase_perc"=sum(increases_6_12, na.rm = TRUE)/n_obs)%>%
  filter(!duplicated(id), !duplicated(increase_perc))

defect_abs <- c("CSP", "EBA140", "Etramp5", "GEXP", "H103", "HSP40",  "Hyp2", "GLURP", "MSP1", "Rh2", "SBP1")


binary_defect_plot <- binary_defect_data %>%
  filter(antigen %in%defect_abs)%>%
  ggplot(., aes(x=factor(any_inf_0_6), y=increase_perc, fill=factor(incidence_value)))+
  geom_bar(stat="identity", color="black")+
  geom_text(aes(label=paste0("frac(",n_obs*increase_perc, ",", n_obs,")")),parse = TRUE, vjust= -0.2, size=3)+
  facet_wrap(~antigen, nrow=2)+
  scale_y_continuous(limits = c(0,1.25), breaks = seq(0,1,by=0.25), labels = scales::percent_format())+
  scale_fill_manual(values=n_infection_cols)+
  ylab("Concentration")+
  xlab("Infections in Months 0-6")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1))

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/binary_defect_plot.png", binary_defect_plot, height=6, width=8, bg="white", dpi=444)

# B memory looks a bit more dodge ####

early_vs_late_exposure_t_memory <- long_tfh_clinab %>%
  filter(!is.na(symp_12_18), cell_pop %in% c("Th_memory","Th_naive"))%>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_12_18), y=cell_freq/100))+
  stat_summary(aes(group=factor(inf_0_12)), position = position_dodge(width=0.75), fun.y = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1)+
  geom_point(shape=21, aes(fill=factor(inf_0_12)), position = position_dodge(width=0.75))+
  # geom_boxplot(aes(fill=factor(inf_0_12)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  ggtitle("lower T cell memory differentiation is associated with increased\nfuture malaria even when controlling for early life exposure")+
  ylab("% of CD4 T cells")+
  xlab("Number of symptomatic episodes in months 12-18")+
  scale_y_continuous(labels = scales::label_percent())+
  guides(fill=guide_legend(title = "# of Infections\nMonths 0-12"))+
  theme_minimal()+
  theme( legend.title = element_text(),
         # axis.text.x = element_text(angle = 90, hjust=1),
         strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)+
  scale_color_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/early_vs_late_exposure_t_memory.png", early_vs_late_exposure_t_memory, height=3, width=6, bg="white", dpi=444)


early_vs_late_exposure_b_memory <- long_abc_clinab %>%
  filter(!is.na(symp_12_18), cell_pop %in% c("memory_b","naive_b"))%>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(symp_12_18), y=cell_freq/100))+
  stat_summary(aes(group=factor(inf_0_12)), position = position_dodge(width=0.75), fun.y = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1)+
  geom_point(shape=21, aes(fill=factor(inf_0_12)), position = position_dodge(width=0.75))+
  # geom_boxplot(aes(fill=factor(inf_0_12)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  ggtitle("higher B cell memory differentiation is associated with increased\nfuture malaria even when controlling for early life exposure")+
  ylab("% of B cells")+
  xlab("Number of symptomatic episodes in months 12-18")+
  scale_y_continuous(labels = scales::label_percent())+
  guides(fill=guide_legend(title = "# of Infections\nMonths 0-12"))+
  theme_minimal()+
  theme( legend.title = element_text(),
         # axis.text.x = element_text(angle = 90, hjust=1),
         strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)+
  scale_color_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/early_vs_late_exposure_b_memory.png", early_vs_late_exposure_b_memory, height=3, width=6, bg="white", dpi=444)

# antibody levels through time juxtaposed between kids that have any malaria in first year of life and none ####
combo_data %>%
  mutate("any_inf_0_12"=if_else(inf_0_12==0, "none", "some"))%>%
  filter(antigen %in% antigens_exposure_6, is.finite(conc))%>%
  ggplot(., aes(x=timepointf, y=conc))+
  geom_hline(data=filter(cutoff_df, antigen %in% antigens_exposure_6), aes(yintercept = log10(max_below_standard)), linetype="dashed")+
  # geom_point(alpha=0.2, shape=21, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=any_inf_0_12), outlier.shape = NA)+
  facet_grid(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("")+
  scale_y_log10()+
  ylab("Concentration")+
  theme_minimal()+
  theme(
    # legend.position="none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1),
    axis.title.x = element_blank(),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

# [DATA CONTAINS ONLY 55 INDIVIDUALS FROM WHOM BOTH B AND T CELL DATA IS AVAILABLE] ####
long_true_combo <- long_true_combo %>%
  mutate(cell_popf=case_when(cell_pop=="activated_b"~"activated B Cells",
                           cell_pop=="atypical_b"~"atypical B Cells",
                           cell_pop=="memory_b"~"memory B Cells",
                           cell_pop=="naive_b"~"naive B Cells",
                           cell_pop=="Tfh_perc"~"Tfh % of CD4",
                           cell_pop=="Tfh_Th1"~"Th1-like of Tfh",
                           cell_pop=="Tfh_Th2"~"Th2-like of Tfh",
                           cell_pop=="Tfh_Th17"~"Th17-like of Tfh",
                           cell_pop=="Tfh_Th1_Th17"~"Th1/Th17-like of Tfh",
                           cell_pop=="Th_naive"~"naive T Cells",
                           cell_pop=="Th_memory"~"memory T Cells"))%>%
  mutate(cell_popf = factor(cell_popf, levels=c("naive B Cells", "memory B Cells", "activated B Cells", "atypical B Cells",
                                                "Th1-like of Tfh", "Th2-like of Tfh", "Th17-like of Tfh", "Th1/Th17-like of Tfh",
                                                "Tfh % of CD4", "naive T Cells", "memory T Cells")))


all_malaria_cells_012 <- long_true_combo %>%
  filter(timepointf=="12 Months")%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_0_12)), outlier.shape = NA)+
  facet_wrap(~cell_popf, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", ncol = 4)+
  ggtitle("")+
  ylab("% of parent population")+
  xlab("Number of Infections in the first year of life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols[])

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/all_malaria_cells_012.png", all_malaria_cells_012, height=8, width=7, bg="white", dpi=444)

sig <- sigs %>%
  filter(Incidence_Type =="inf_0_12")%>%
  select(cell_pop)

sig_malaria_cells_012 <- long_true_combo %>%
  filter(timepointf=="12 Months", cell_pop %in% c("memory_b", "Th_memory", "Tfh_Th1", "Tfh_Th1_Th17"))%>%
  ggplot(., aes(x=factor(inf_0_12), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_0_12)), outlier.shape = NA)+
  facet_wrap(~factor(cell_popf, levels=c("memory B Cells", "memory T Cells", "Th1-like of Tfh", "Th1/Th17-like of Tfh")), labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", nrow=1)+
  ggtitle("early life malaria is associated with increased memory differentiation in B cells,\nbut not T cells, and a slight shift in Tfh differentiation")+
  ylab("% of parent population")+
  xlab("Number of Infections in the first year of life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols[])

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/sig_malaria_cells_012.png", sig_malaria_cells_012, height=3, width=8, bg="white", dpi=444)




any_012_malaria_cells <- long_true_combo %>%
  filter(timepointf=="12 Months")%>%
  mutate(any_inf_0_12=case_when(any_0_12=="no"~"none",
                                any_0_12=="yes"~"some"))%>%
  ggplot(., aes(x=factor(any_inf_0_12), y=cell_freq/100))+
  # geom_point(alpha=0.2, shape=21)+
  geom_violin(aes(fill=factor(any_inf_0_12)))+
  facet_wrap(~cell_pop, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", ncol = 4)+
  ggtitle("")+
  ylab("% of parent population")+
  xlab("Number of Infections in the first year of life")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/remix/figures_for_paper/any_012_malaria_cells.png", any_012_malaria_cells, height=8, width=6, bg="white", dpi=444)


big_early_life_malaria_cells06 <- long_true_combo %>%
  filter(timepointf=="12 Months")%>%
  # mutate(any_inf_0_6=case_when(any_inf_0_6=="0"~"none",
  #                              any_inf_0_6=="1"~"some"))%>%
  ggplot(., aes(x=factor(inf_0_6), y=cell_freq/100))+
  geom_point(alpha=0.2, shape=21)+
  geom_boxplot(aes(fill=factor(inf_0_6)), outlier.shape = NA)+
  facet_wrap(~cell_pop, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free", ncol = 4)+
  ggtitle("")+
  ylab("% of Parent Population")+
  xlab("Number of Infections in Months 0-6")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust=1),
    strip.text = element_text())+
  scale_fill_manual(values=n_infection_cols[])
