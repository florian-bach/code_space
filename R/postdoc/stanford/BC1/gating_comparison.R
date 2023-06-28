library(dplyr)
library(tidyr)
library(ggplot2)
# old data ####


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



# new data ####

files <- list.files("~/Downloads/drive-download-20230620T190355Z-001", full.names = TRUE)

tfh_df <- data.frame()

for (i in files){
  
  df <- read.csv(i)
  tfh_df <- rbind(tfh_df, df)
  
  
}

tfh_df <- subset(tfh_df, tfh_df$X != "Mean" & tfh_df$X != "SD")
# colnames(tfh_df)[1:6]<- c("file_name", "Tfh_of_CD4", "Th1", "Th1_17", "Th17", "Th2")
colnames(tfh_df)[1:6]<- c("file_name", "tfh_of_cd4", "tfh_q1", "tfh_q2", "tfh_q3", "tfh_q4")

tfh_df <- tfh_df[,1:6]
tfh_df$id <- paste(1, substr(tfh_df$file_name, 15, 18), sep="")
tfh_df$gating <- "new"


subset_old <- combo_cells_data %>%
  filter(id %in% as.numeric(tfh_df$id)) %>%
  select(colnames(tfh_df)[2:7])%>%
  mutate("gating"="old")

# put it all together
combo_data <- rbind(subset_old, tfh_df[,-1])

long_combo_data <- combo_data %>%
  pivot_longer(cols = colnames(combo_data)[1:5], names_to = "cell_pop", values_to = "perc")%>%
  mutate("freq"=perc/100)%>%
  mutate("recode_cell_pop" = if_else(gating=="new",
                                     recode(cell_pop, "tfh_of_cd4"="% cTfh of CD4 T cells", "tfh_q1"="Th1", "tfh_q2" = "Th1 / Th17", "tfh_q3"="Th17","tfh_q4"="Th2"),
                                     recode(cell_pop, "tfh_of_cd4"="% cTfh of CD4 T cells", "tfh_q1"="Th17", "tfh_q2" ="Th1 / Th17", "tfh_q3"="Th1","tfh_q4"="Th2")
  ))


old_combo_data <- long_combo_data %>%
  filter(gating=="old")%>%
  mutate("old_freq"=perc/100)

new_combo_data <- long_combo_data %>%
  filter(gating=="new")%>%
  mutate("new_freq"=perc/100)

better_combo <- inner_join(old_combo_data,new_combo_data, by=c("id", "recode_cell_pop"))


pop_palette <- c("#81ebb6", "#f0d330", "#cc1dc9", "#70a7e6", "#d67c1c")

old_new_scatter <- ggplot(better_combo, aes(x=old_freq, y=new_freq, fill=recode_cell_pop))+
  geom_point(size=2, shape=21)+
  geom_blank(aes(x = new_freq, y = old_freq)) +
  scale_y_continuous(labels=scales::label_percent())+
  scale_x_continuous(labels=scales::label_percent())+
  geom_abline(slope = 1, intercept = 0, color="grey", linetype="dashed")+
  facet_wrap(~recode_cell_pop, scales = "free")+
  scale_fill_manual(values = pop_palette)+
  #facet_grid(~cell_pop, space="fixed", scales="free")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/old_new_scatter.png", old_new_scatter, height=5, width=8, bg="white")



old_new_mirror <- long_combo_data %>%
  mutate(inv_perc=if_else(gating=="old", perc*(-1), perc))%>%
  group_by(recode_cell_pop)%>%
  ggplot(aes(x=factor(id), y=inv_perc/100, fill=id))+
    geom_bar(stat="identity", position="dodge")+
    geom_blank(aes(x = factor(id), y = inv_perc/(-100))) +
    scale_y_continuous(labels=scales::label_percent())+
    ggtitle("spotcheck old (bottom) vs new (top) gating")+
    geom_hline(yintercept = 0)+
    xlab("")+
    ylab("frequency")+
    facet_wrap(~recode_cell_pop, scales = "free")+
    #scale_fill_manual(values = pop_palette)+
    #facet_grid(~cell_pop, space="fixed", scales="free")+
    theme_minimal()+
    theme(legend.position = "none",
        axis.text.x = element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/old_new_mirror.png", old_new_mirror, height=6, width=8, bg="white")


old_new_box <- ggplot(long_combo_data, aes(x=gating, y=freq, fill=recode_cell_pop))+
  geom_boxplot()+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(~recode_cell_pop, scales = "free")+
  scale_fill_manual(values = pop_palette)+
  ggtitle("spotcheck old vs new gating")+
  ylab("cell frequency")+
  
  #facet_grid(~cell_pop, space="fixed", scales="free")+
  theme_minimal()+
  theme(
        legend.title = element_blank()
    )

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/old_new_boxplot.png", old_new_box, height=6, width=8, bg="white")



old_new_vln <- ggplot(long_combo_data, aes(x=gating, y=freq, fill=recode_cell_pop))+
  geom_violin()+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(~recode_cell_pop, scales = "free")+
  scale_fill_manual(values = pop_palette)+
  ggtitle("spotcheck old vs new gating")+
  ylab("cell frequency")+
  
  #facet_grid(~cell_pop, space="fixed", scales="free")+
  theme_minimal()+
  theme(
    legend.title = element_blank()
  )

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/old_new_vln.png", old_new_vln, height=6, width=8, bg="white")



old_new_line <- ggplot(long_combo_data, aes(x=factor(gating, levels=c("old", "new")), y=freq, color=recode_cell_pop, group=id))+
  geom_point()+
  geom_line()+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(~recode_cell_pop, scales = "free", nrow=1)+
  scale_fill_manual(values = pop_palette)+
  ggtitle("spotcheck old vs new gating")+
  ylab("cell frequency")+
  
  #facet_grid(~cell_pop, space="fixed", scales="free")+
  theme_minimal()+
  theme(
    legend.position = "none",
    axis.title = element_blank()
  )

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/old_new_line.png", old_new_line, height=4, width=8, bg="white")


all_old <- combo_cells_data %>%
  pivot_longer(cols = matches('tfh_q[0-9]'), names_to = "cell_pop", values_to = "freq")%>%
  mutate("recode_cell_pop"=recode(cell_pop, "tfh_of_cd4"="% cTfh of CD4 T cells", "tfh_q1"="Th17", "tfh_q2" ="Th1 / Th17", "tfh_q3"="Th1","tfh_q4"="Th2"))%>%
  ggplot(aes(x="", y=freq/100, fill=recode_cell_pop))+
  geom_boxplot()+
  ggtitle("all old data")+
  ylab("cell frequency")+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(~recode_cell_pop, scales = "free", nrow=1)+
  scale_fill_manual(values = pop_palette)+
  #facet_grid(~cell_pop, space="fixed", scales="free")+
  theme_minimal()+
  theme(legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x= element_blank())

ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/all_old_boxplot.png", all_old, height=3, width=8, bg="white")

# correlate cell frequencies and antibody concentrations ####


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
tfh_df$id <- as.numeric(tfh_df$id)
twelve_cor_combo <- inner_join(wide_combo_12months, tfh_df, by="id")


long_twelve <- twelve_cor_combo%>%
  pivot_longer(cols = matches('tfh'), names_to = "cell_pop", values_to = "perc")%>%
  mutate("freq"=perc/100)%>%
  mutate("recode_cell_pop"=recode(cell_pop, "tfh_of_cd4"="% cTfh of CD4 T cells", "tfh_q1"="Th17", "tfh_q2" ="Th1 / Th17", "tfh_q3"="Th1","tfh_q4"="Th2"))%>%
  pivot_longer(cols = colnames(twelve_cor_combo)[4:24], names_to = "antigen", values_to = "ab_conc")


big_plot <- ggplot(long_twelve, aes(x=freq, y=ab_conc))+
  geom_point(aes(color=factor(id)))+
  scale_x_continuous(labels=scales::label_percent())+
  facet_wrap(recode_cell_pop~antigen, scales="free", ncol=21)+
  geom_smooth(method="lm", se = TRUE)+
  theme_minimal()+
  xlab("cell frequency")+
  ylab("antibody concentration")+
  theme(legend.position="none")


ggsave("~/postdoc/stanford/clinical_data/BC1/antibody_modelling/figures/new_gating_ab_corr.png", big_plot, width=40, height=12, bg="white", limitsize = FALSE)


cleaned_broomer <- long_twelve %>%
  group_by(recode_cell_pop, antigen)%>%
  do(broom::tidy(cor.test(.$freq, .$ab_conc, method="spearman")))%>%
  ungroup()%>%
  mutate("p_adj"=p.adjust(p.value))%>%
  ungroup()

sig_cleaned_broomer <- filter(cleaned_broomer, p.value<0.1)



