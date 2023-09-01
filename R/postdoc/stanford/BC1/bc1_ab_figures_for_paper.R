# preamble & data generation ####
#palettes
time_palette <- colorspace::sequential_hcl(n=4, "RdPu")[1:3]
pc1_cols <- colorspace::sequential_hcl(23, palette = "Purple Yellow")
incidence_cols <- colorspace::sequential_hcl(11, palette = "Purple Yellow")
n_infection_cols <- c("white", colorspace::sequential_hcl(n=5, palette = "Lajolla")[-1])
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

ab_columns <- grep("log", colnames(bc1), value = TRUE)


long_raw_dfff <- bc1 %>%
  dplyr::select(all_of(c("id", "timepoint", ab_columns, "MomFinalRx", "anyHPfinal", "gestage", "gender")))%>%
  pivot_longer(cols=all_of(ab_columns), names_to = "antigen", values_to = "conc")%>%
  filter(antigen %notin% c("logpd", "logGST"))%>%
  mutate(antigen=gsub("log", "", antigen, fixed = TRUE))%>%
  mutate(antigen=gsub("_", " ", antigen, fixed = TRUE))%>%
  mutate(MomFinalRx=if_else(MomFinalRx==1, "3 Dose SP",
                             ifelse(MomFinalRx==2, "3 Dose DP",
                                    if_else(MomFinalRx==3, "Monthly DP", "NA")))
         
  )%>%
  mutate(MomFinalRx=factor(MomFinalRx, levels = c("3 Dose SP", "3 Dose DP", "Monthly DP")))%>%
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


tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
abc <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/atypical_b.csv")

#fix & harmonise colnames, add id column that matches rest of data
cleaned_tfh <- tfh
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
cleaned_tfh$idt <- paste("1", substr(cleaned_tfh$idt, 15, 18), sep='')

#fix & harmonise colnames, add id column that matches rest of data
cleaned_abc <- abc
colnames(cleaned_abc) <- c("idb", "abc_q4", "abc_q3","abc_q2","abc_q1","cd10_neg_b", "parent2b", "abc_drop")
cleaned_abc$cd10_pos_b <- 100-cleaned_abc$cd10_neg_b
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

very_long_combo_combo <- long_combo_combo %>%
  filter(tfh_drop=="No", abc_drop=="No")%>%
  dplyr::select(-tfh_drop, -abc_drop) %>%
  pivot_longer(cols = matches('tfh|abc|cd10'), names_to = "cell_pop", values_to = "cell_freq") %>%
  mutate(cell_pop = recode(cell_pop, 
                           "abc_q1"="Activated B Cells",
                           "abc_q2"="Memory B Cells",
                           "abc_q3"="Naive B Cells",
                           "abc_q4"="Atypical B Cells",
                           #"cd10_neg_b"="Mature B Cells",
                           "cd10_pos_b"="Immature B Cells",
                           "tfh_q1"="Th17",
                           "tfh_q2"="Th1_17",
                           "tfh_q3"="Th1",
                           "tfh_q4"="Th2",
                           "tfh_of_cd4"="Tfh of CD4 T Cells"))%>%
  filter(cell_pop != "cd10_neg_b")

# wide_combo_12months <- combo_data |> 
#   dplyr::select(id, timepoint, anyHPfinal, antigen, conc)|>
#   filter(timepoint==3)|>
#   pivot_wider(names_from = antigen, values_from = conc)
# 
# twelve_cor_combo <- inner_join(wide_combo_12months, combo_cells_data, by="id")

twelve_cor_combo <- very_long_combo_combo %>%
  filter(timepoint==3)%>%
  select(id, timepoint, anyHPfinal, antigen, conc, cell_pop, cell_freq)%>%
  pivot_wider(names_from = cell_pop, values_from = cell_freq)%>%
  pivot_wider(names_from = antigen, values_from = conc)

# antibody temporal trends ####

antibodies_during_first_year_of_life <- ggplot(combo_data, aes(x=timepointf, y=conc, fill=antigen))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(alpha=0.1, shape=21)+
  geom_line(aes(group=id), alpha=0.1)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  facet_wrap(~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow = 3, scales = "free")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_during_first_year_of_life.png", antibodies_during_first_year_of_life, height = 12, width=12, bg="white")


# antibody~incidence figures ####

antibodies_incidence_cord <- combo_data %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==1, Incidence_Type %in% c("inf_0_6", "symp_0_6", "any_inf_0_6", "any_symp_0_6"))%>%
  ggplot(., aes(x=Incidence_Value, y=conc))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=Incidence_Value))+
  facet_grid(Incidence_Type~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Cord Blood")+
  ylab("Concentration")+
  xlab("Number of Infections")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_incidence_cord.png", antibodies_incidence_cord, height = 8, width=20, bg="white", limitsize = FALSE)






antibodies_incidence_six <- combo_data %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==2, Incidence_Type %in% c("inf_0_6", "symp_0_6", "any_inf_0_6", "any_symp_0_6"))%>%
  ggplot(., aes(x=Incidence_Value, y=conc))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=Incidence_Value))+
  facet_grid(Incidence_Type~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Incidence Before Six Months")+
  ylab("Concentration")+
  xlab("Number of Infections")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_incidence_six.png", antibodies_incidence_six, height = 8, width=20, bg="white", limitsize = FALSE)




antibodies_incidence_six_post <- combo_data %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==2, Incidence_Type %in% c("inf_6_12", "symp_6_12", "any_inf_6_12", "any_symp_6_12"))%>%
  ggplot(., aes(x=Incidence_Value, y=conc))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_boxplot(aes(fill=antigen, group=Incidence_Value))+
  scale_x_continuous(breaks = c(0:3))+
  facet_grid(Incidence_Type~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  ggtitle("Incidence After Six Months")+
  ylab("Concentration")+
  xlab("Number of Infections")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_incidence_six_post.png", antibodies_incidence_six_post, height = 8, width=20, bg="white", limitsize = FALSE)



antibodies_incidence_twelve <- combo_data %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==2, Incidence_Type %in% c("inf_6_12", "symp_6_12", "any_inf_6_12", "any_symp_6_12"))%>%
  ggplot(., aes(x=Incidence_Value, y=conc))+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=Incidence_Value))+
  facet_grid(Incidence_Type~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Incidence Before Twelve Months")+
  xlab("Number of Infections")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_incidence_twelve.png", antibodies_incidence_twelve, height = 8, width=20, bg="white", limitsize = FALSE)



antibodies_incidence_twelve_post <- combo_data %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==3, Incidence_Type %in% c("inf_12_18", "symp_12_18", "any_inf_12_18", "any_symp_12_18"))%>%
  ggplot(., aes(x=Incidence_Value, y=conc))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=Incidence_Value))+
  facet_grid(Incidence_Type~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Incidence After Twelve Months")+
  xlab("Number of Infections")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/antibodies_incidence_twelve_post.png", antibodies_incidence_twelve_post, height = 8, width=20, bg="white", limitsize = FALSE)


# antibody~placental malaria, preterm birth, maternal chemoprevention ####



histopathology <- combo_data %>%
  filter(!is.na(anyHPfinalx))%>%
  ggplot(., aes(x=anyHPfinalx, y=conc))+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=anyHPfinalx))+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Histopathology")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/histopathology.png", histopathology, height = 8, width=20, bg="white", limitsize = FALSE)





MomFinalRx <- combo_data %>%
  ggplot(., aes(x=MomFinalRx, y=conc))+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=MomFinalRx))+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Maternal Chemoprevention")+
  xlab("Number of Infections")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/MomFinalRx.png", MomFinalRx, height = 8, width=20, bg="white", limitsize = FALSE)





gestages_bins <- combo_data %>%
  filter(!is.na(gestage))%>%
  mutate(gestagef=factor(if_else(gestage<28, "<28 weeks", 
                          if_else(gestage>=28 & gestage<32, "28-32 weeks", 
                                  if_else(gestage>=32 & gestage<37, "32-37 weeks", 
                                          if_else(gestage>=37, ">37 weeks", "whoops")))), levels=c("<28 weeks", "28-32 weeks", "32-37 weeks", ">37 weeks")))%>%
  ggplot(., aes(x=gestagef, y=conc))+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=antigen, group=gestagef))+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Gestational Age")+
  xlab("Number of Infections")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/gestages_bins.png", gestages_bins, height = 8, width=20, bg="white", limitsize = FALSE)



gestages_cont <- combo_data %>%
  filter(!is.na(gestage))%>%
  ggplot(., aes(x=gestage, y=conc))+
  geom_smooth(method="lm")+
  geom_point(aes(fill=antigen), alpha=0.1, shape=21)+
  facet_grid(timepointf~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Gestational Age")+
  xlab("\nGestational Age (Weeks)")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/gestages_cont.png", gestages_cont, height = 8, width=20, bg="white", limitsize = FALSE)


# antibodies ~ cell frequencies ####
cor_data <- twelve_cor_combo |>
  dplyr::select(-id, -timepoint, -anyHPfinal)%>%
  na.omit()

spearman <- cor(cor_data, method = "spearman")
twelve_dist <- dist(spearman, method = "euclid", diag = FALSE, upper = FALSE, p = 2)
twelve_hclust <- hclust(twelve_dist)

cor_mat <- as.matrix(spearman)


# col_fun_spearment <- circlize::colorRamp2(c(min(cor_mat), 0, max(cor_mat)), c("#0859C6", "black", "#FFA500"))

col_fun_spearment <- circlize::colorRamp2(c(-1, 0, 1), c("#0859C6", "black", "#FFA500"))

#doesnt work for pdf..
# colnames(spearment_matrix) <- gsub("IFNy", "IFNγ", colnames(spearment_matrix))
# rownames(spearment_matrix) <-  gsub("IFNy", "IFNγ",rownames(spearment_matrix))

spearment_heatmap <- Heatmap(matrix = cor_mat,
                           cluster_rows = TRUE,
                           cluster_columns=TRUE,
                           show_row_dend = FALSE,
                           show_column_dend = TRUE,
                           show_heatmap_legend = TRUE,
                           name = "Spearmon rho",
                           #name = "spearment r",
                           #cluster_columns = FALSE,
                           column_names_gp = gpar(fontsize = 9),
                           row_names_gp = gpar(fontsize = 9),
                           row_names_side = "left",
                           col = col_fun_spearment,
                           column_names_rot = 45,
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 5, col="white"))
                           })


png("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/ab_cell_correlation_heatmap.png", width = 8, height=8, units="in", res=700)
draw(spearment_heatmap)
dev.off()


spearment_heatmap_no_lab <- Heatmap(matrix = cor_mat,
                             cluster_rows = TRUE,
                             cluster_columns=TRUE,
                             show_row_dend = FALSE,
                             show_column_dend = TRUE,
                             show_heatmap_legend = TRUE,
                             name = "Spearmon rho",
                             #name = "spearment r",
                             #cluster_columns = FALSE,
                             column_names_gp = gpar(fontsize = 9),
                             row_names_gp = gpar(fontsize = 9),
                             row_names_side = "left",
                             col = col_fun_spearment,
                             column_names_rot = 45
                             )


png("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/ab_cell_correlation_heatmap_no_lab.png", width = 8, height=8, units="in", res=700)
draw(spearment_heatmap_no_lab)
dev.off()



cleaned_broomer <- very_long_combo_combo %>%
  filter(!is.na(antigen), !is.na(cell_pop), timepoint==3)%>%
  group_by(cell_pop, antigen)%>%
  do(broom::tidy(cor.test(.$cell_freq, .$conc, method="spearman")))%>%
  group_by(cell_pop)%>%
  mutate("p_adj"=p.adjust(p.value))%>%
  ungroup()


fdr_cutoff <- 0.05
sig_cleaned_broomer <- filter(cleaned_broomer, p.value<fdr_cutoff)


list_of_plots <- list(matrix(nrow = nrow(sig_cleaned_broomer)))

qual_palette <- viridis::inferno(nrow(sig_cleaned_broomer))

x_and_ys <- cbind(sig_cleaned_broomer$cell_pop, sig_cleaned_broomer$antigen)

for(i in 1:nrow(sig_cleaned_broomer)){
  
  plot <- very_long_combo_combo %>%
    filter(!is.na(antigen), timepoint==3) %>%
    filter(cell_pop==x_and_ys[i,1], antigen==x_and_ys[i,2])%>%
    ggplot(., aes(x=cell_freq, y=conc))+
    geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    geom_smooth(method="lm")+
    ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    xlab(x_and_ys[i,1])+
    ylab(x_and_ys[i,2])+
    theme_minimal()
  
  list_of_plots[[i]] <- plot
  
}

big_plot <- cowplot::plot_grid(plotlist = list_of_plots, nrow =  round(nrow(sig_cleaned_broomer)/5))

ggsave(filename = paste("~/postdoc/stanford/clinical_data/BC1/figures_for_paper/old_gating_big_indie_ab_cell_correlation_plot", fdr_cutoff, "raw_p.png", sep="_"), big_plot, height = 6, width=15, bg="white")


# big plots of all cells vs all antibodies

very_long_combo_combo%>%
  filter(cell_pop=="Th2", timepoint==3)%>%
  ggplot(., aes(x=cell_freq, y=conc))+
  geom_point(aes(fill=cell_pop), alpha=1, shape=21)+
  scale_fill_manual(values = pc1_cols)+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
  facet_wrap(~antigen, scales="free")+
  theme_minimal()


# bootstrapping spearman correlations

bla <- RVAideMemoire::spearman.cor.multcomp(very_long_combo_combo$cell_freq, very_long_combo_combo$conc, fact=interaction(very_long_combo_combo$cell_pop, very_long_combo_combo$antigen), alpha = 0.05, nrep = 1000)

bootstrap <- bla$tab

# cell frequencies at 1 year of age ####

all_cells <- very_long_combo_combo %>%
  ggplot(., aes(x="", y=cell_freq, fill=cell_pop))+
  geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  facet_wrap(~cell_pop, labeller = labeller(cell_pop = label_wrap_gen(width = 8)), nrow = 2, scales = "free")+
  ylab("Percentage of Parent")+
  ggtitle("All Cells")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/all_cells.png", all_cells, height = 8, width=8, bg="white")



all_cells_no_drop <- very_long_combo_combo %>%
  ggplot(., aes(x="", y=cell_freq, fill=cell_pop))+
  geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  facet_wrap(~cell_pop, labeller = labeller(cell_pop = label_wrap_gen(width = 8)), nrow = 2, scales = "free")+
  ylab("Percentage of Parent")+
  ggtitle("without dropped data")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/all_cells_no_drop.png", all_cells_no_drop, height = 8, width=8, bg="white")


# cell frequencies and malaria incidence ####

cells_incidence <- very_long_combo_combo %>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==1, Incidence_Type %in% c("inf_6_12", "symp_6_12", "inf_12_18", "symp_12_18"))%>%
  ggplot(., aes(x=Incidence_Value, y=cell_freq))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(aes(fill=cell_pop), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=cell_pop, group=Incidence_Value))+
  ggh4x::facet_nested_wrap(facets=c("Incidence_Type", "cell_pop"), labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow=4, scales = "free_y")+
  ggtitle("Malaria Incidence and Cell Frequencies (All Cells)")+
  ylab("Percentage of Parent")+
  xlab("Number of Infections")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/cells_incidence.png", cells_incidence, height = 8, width=20, bg="white", limitsize = FALSE)



  
cells_incidence_no_drop <- very_long_combo_combo%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  filter(timepoint==1, Incidence_Type %in% c("inf_6_12", "symp_6_12", "inf_12_18", "symp_12_18"))%>%
  ggplot(., aes(x=Incidence_Value, y=cell_freq))+
  #geom_violin(draw_quantiles = seq(0,1,0.25), scale = TRUE, )+
  geom_point(aes(fill=cell_pop), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=cell_pop, group=Incidence_Value))+
  ggh4x::facet_nested_wrap(facets=c("Incidence_Type", "cell_pop"), labeller = labeller(antigen = label_wrap_gen(width = 6)), nrow=4, scales = "free_y")+
  ggtitle("Malaria Incidence and Cell Frequencies (without dropped data)")+
  ylab("Percentage of Parent")+
  xlab("Number of Infections")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        
        axis.text.x = element_text(angle = 0, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/cells_incidence_no_drop.png", cells_incidence_no_drop, height = 8, width=20, bg="white", limitsize = FALSE)

# cells~placental malaria, pretern birth, maternal chemoprevention

histopathology_cells <- very_long_combo_combo %>%
  dplyr::select(-antigen, -conc) %>%
  filter(timepoint==3)%>%
  filter(!duplicated(id), !is.na(anyHPfinal))%>%
  # pivot_longer(cols = matches('tfh|abc|cd10'), names_to = "cell_pop", values_to = "cell_freq") %>%
  # mutate(cell_pop = recode(cell_pop, 
  #                          "abc_q1"="Activated B Cells",
  #                          "abc_q2"="Memory B Cells",
  #                          "abc_q3"="Naive B Cells",
  #                          "abc_q4"="Atypical B Cells"))%>%
  ggplot(., aes(x=anyHPfinalx, y=cell_freq))+
  geom_point(aes(fill=anyHPfinalx), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=anyHPfinalx))+
  facet_wrap(~cell_pop, scales = "free_y", nrow=2)+
  ggtitle("Histopathology")+
  xlab("\nHistopathology")+
  ylab("Frequency")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=rev(n_infection_cols))
  

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/histopathology_cells.png", histopathology_cells, height = 8, width=12, bg="white", limitsize = FALSE)



gestages_cont_cells <- very_long_combo_combo %>%
  filter(timepoint==3)%>%
  dplyr::select(-antigen, -conc)%>%
  group_by(id)%>%
  filter(!duplicated(cell_pop))%>%
  ggplot(., aes(x=gestage, y=cell_freq))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  geom_point(aes(fill=cell_pop), alpha=1, shape=21)+
  facet_wrap(~cell_pop, scales = "free_y", nrow=2)+
  ggtitle("Gestational Age")+
  xlab("\nGestational Age (Weeks)")+
  ylab("Frequency")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=pc1_cols)

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/gestages_cont_cells.png", gestages_cont_cells, height = 8, width=12, bg="white", limitsize = FALSE)




momrx_cells <- very_long_combo_combo %>%
  filter(timepoint==3)%>%
  dplyr::select(-antigen, -conc)%>%
  ggplot(., aes(x=MomFinalRx, y=cell_freq))+
  # geom_smooth(method="lm")+
  # ggpubr::stat_cor(method = "spearman", na.rm = TRUE)+
  geom_point(aes(fill=MomFinalRx), alpha=0.1, shape=21)+
  geom_boxplot(aes(fill=MomFinalRx))+
  facet_wrap(~cell_pop, scales = "free_y", nrow=2)+
  ggtitle("MomFinalRx")+
  xlab("\nMomFinalRx")+
  ylab("Frequency")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))+
  scale_fill_manual(values=rev(n_infection_cols))

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/momrx_cells.png", momrx_cells, height = 8, width=12, bg="white", limitsize = FALSE)

# gestational age, histopathology malaria incidence

gestages_incidence <- combo_data %>%
  filter(!is.na(gestage))%>%
  pivot_longer(cols = matches("symp|inf"), names_to = "Incidence_Type", values_to = "Incidence_Value")%>%
  # mutate(gestagef=factor(if_else(gestage<28, "<28 weeks", 
                                 # if_else(gestage>=28 & gestage<32, "28-32 weeks", 
                                 #         if_else(gestage>=32 & gestage<37, "32-37 weeks", 
                                 #                 if_else(gestage>=37, ">37 weeks", "whoops")))), levels=c("<28 weeks", "28-32 weeks", "32-37 weeks", ">37 weeks")))%>%
  ggplot(., aes(x=Incidence_Value, y=gestage, fill=factor(Incidence_Value)))+
  #geom_violin()+
  geom_boxplot(aes(fill=factor(Incidence_Value), group=Incidence_Value))+
  geom_point(aes(fill=factor(Incidence_Value)), alpha=0.1, shape=21)+
  facet_wrap(~Incidence_Type, labeller = labeller(antigen = label_wrap_gen(width = 6)), scales = "free")+
  ggtitle("Gestational Age")+
  xlab("Number of Infections")+
  ylab("Gestational Age (Weeks)")+
  theme_minimal()+
  scale_fill_manual(values = pc1_cols)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        strip.text = element_text(size=13.5))

ggsave("/Users/fbach/postdoc/stanford/clinical_data/BC1/figures_for_paper/gestages_incidence.png", gestages_incidence, height = 10, width=22, bg="white", limitsize = FALSE)

