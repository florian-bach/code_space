library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


`%notin%` <- Negate(`%in%`)

time_cols <- list("baseline"="#E4DEBD",
                  "day0" = "#C03F3E",
                  "day7" = "#D87E1F",
                  "day14" = "#E6B85F")


metadata <- read.csv("~/Library/CloudStorage/Box-Box/Border Cohort Immunology (MUSICAL)/Data/musical_metadata_new_qpcr.csv")


metadata <- metadata%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %in% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"baseline",
                               timepoint_imm==-1 & id %notin% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"baseline",
                               timepoint_imm==-1 & infectiontype %in% c("A", "S") & id %in% c(134, 164, 176, 317, 331, 363, 379, 398, 402, 442, 577) ~"bad_baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))
  
#subset to only include musical folks

blood_counts <- readxl::read_xls("~/postdoc/stanford/clinical_data/MUSICAL/PRISMBC LYMPHOCYTE TRUCOUNTS.xls")

slim_blood_counts <- blood_counts%>%
  select(Id, `Draw date`, `CD3+  cells/ul blood`, `lymphocytes/ul blood`)%>%
  filter(!is.na(Id))%>%
  mutate(date=as.Date(`Draw date`), id=as.character(Id))%>%
  select(-`Draw date`, -Id)


clean_data_with_cell_counts <- metadata %>%
  mutate(date=as.Date(date), id=as.character(id))%>%
  left_join(., slim_blood_counts, by=c("id", "date"))%>%
  pivot_longer(cols=c( `CD3+  cells/ul blood`, `lymphocytes/ul blood`), names_to = "blood_count", values_to = "cbc_per_ul")%>%
  mutate(timepoint = factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))



cbc_boxplot <- clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  ggplot(., aes(x=infectiontype, y=cbc_per_ul, fill=factor(timepoint)))+
  geom_boxplot(outliers = F)+
  # ggpubr::stat_compare_means(paired = T, ref.group = ".all.")+
  facet_wrap(~blood_count)+
  theme_minimal()+
  scale_fill_manual(values=time_cols)

ggsave("~/postdoc/stanford/clinical_data/MUSICAL/cbc_boxplot.png", cbc_boxplot, height = 6, width = 8, dpi=444, bg="white")  

lymph_count_for_paper <-  clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"),
         timepoint %in% c("baseline", "day0", "day7", "day14"),
         blood_count!="lymphocytes/ul blood")%>%
  ggplot(., aes(x=timepoint,y=cbc_per_ul, fill=timepoint))+
  geom_boxplot(outliers = F)+
  # ggpubr::stat_compare_means(paired = T, ref.group = ".all.")+
  facet_wrap(~infectiontype, scales='free_x')+
  theme_minimal()+
  ylab("T cells / μL")+
  scale_fill_manual(values=time_cols)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tcounts.pdf", lymph_count_for_paper, height = 3, width = 6, bg="white", device = cairo_pdf)  



cbc_boxplot2 <- clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(blood_count, cbc_per_ul, id, timepoint, infectiontype)%>%
  filter(!is.na(cbc_per_ul))%>%
  pivot_wider(names_from = blood_count, values_from = cbc_per_ul, id_cols = c("id", "timepoint", "infectiontype"))%>%
  mutate("non T lymphocytes"=`lymphocytes/ul blood`-`CD3+  cells/ul blood`)%>%
  pivot_longer(cols=c(`lymphocytes/ul blood`, `CD3+  cells/ul blood`, `non T lymphocytes`), names_to = "blood_count", values_to = "cbc_per_ul" )%>%
  ggplot(., aes(x=timepoint, y=cbc_per_ul, fill=timepoint))+
  #geom_point()+
  geom_line(aes(group = id), alpha=0.5)+
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~blood_count+infectiontype, ncol=2, scales="free")+
  scale_fill_manual(values=viridis::magma(n=5))+
  da_boxplot_theme


clean_data_with_cell_counts %>%
  filter(infectiontype %in% c("S", "A"), timepoint %in% c("baseline", "day0", "day7", "day14"))%>%
  distinct(blood_count, cbc_per_ul, id, timepoint, infectiontype, parasitedensity)%>%
  pivot_wider(names_from = blood_count, values_from = cbc_per_ul)%>%
  mutate("non T lymphocytes"=`lymphocytes/ul blood`-`CD3+  cells/ul blood`)%>%
  pivot_longer(cols=c(`lymphocytes/ul blood`, `CD3+  cells/ul blood`, `non T lymphocytes`), names_to = "blood_count", values_to = "cbc_per_ul" )%>%
  ggplot(., aes(x=parasitedensity+0.1, y=cbc_per_ul+0.1))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  geom_smooth(method="lm")+
  theme_minimal()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("parasites / ul")+
  ylab("immune cells / ul")+
  facet_wrap(~blood_count+infectiontype, ncol=2, scales="free")+
  scale_fill_manual(values=viridis::magma(n=5))+
  da_boxplot_theme


# add flow cytometry data ####

# jason_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/mus1_mus6_flow_data.csv")
jason_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/df_jason_analysis.csv")

jason_data <- jason_data %>%
  mutate(id=as.numeric(cohortid),
         "timepoint"=case_when(
           timepoint==-1~"baseline",
           timepoint==0~"day0",
           timepoint==7~"day7",
           timepoint==14~"day14",
           timepoint==28~"day28"))%>%
  filter(!is.na(cohortid))%>%
  select(-c(X.2, index, Unnamed..0, X))

# write.csv(jason_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/mus1_mus8_flow_data_edit.csv")

combo_data <- clean_data_with_cell_counts%>%
  mutate(id=as.numeric(id))%>%
  filter(!is.na(cbc_per_ul))%>%
  right_join(., jason_data, by = c("id", "timepoint", "infectiontype"))%>%
  mutate("timepoint"=case_when(timepoint_imm==-2 & id %notin% c("176", "363", "577") ~"bad_baseline",
                               timepoint_imm==-2 & id %in% c("176", "363", "577") ~"baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28"))


cell_count_data <- combo_data%>%
  filter(blood_count=="CD3+  cells/ul blood" | is.na(blood_count))%>%
  mutate("absolute_Tr1"=cbc_per_ul*CD4_T_Cell_Frequency/100* Memory_CD4_T_Cell_Frequency/100 *Tr1_Frequency/100)%>%
  mutate("absolute_CD4_memory"=cbc_per_ul*CD4_T_Cell_Frequency/100* Memory_CD4_T_Cell_Frequency/100)%>%
  pivot_longer(cols = ends_with("Frequency"), names_to = "gate", values_to = "freq")%>%
  mutate(count = ifelse(grepl("Tr1", gate), absolute_Tr1*freq/100, 
                        ifelse(gate %in% c("IFNg_Frequency","IL10_Frequency", "IL21_Frequency"), absolute_CD4_memory*freq/100, NA)))

write.csv(cell_count_data, "~/postdoc/stanford/manuscripts/jason_tr1_2/revised_baselines/t_cell_cytometry_count_data.csv", row.names = F)




tr1_freq_plot <- cell_count_data %>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline")%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=absolute_Tr1*1000, fill=timepoint))+
  geom_line(aes(group=id), alpha=0.2)+
  geom_boxplot(outliers = F)+
  # scale_y_log10()+
  ylab("Tr1 cells / mL")+
  facet_wrap(~infectiontype+stim, scales = "free_x")+
  scale_fill_manual(values=time_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        legend.position = "none")
  
ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/tr1_freq_plot.png", tr1_freq_plot, width = 5.3333, height=3, dpi=444, bg="white")
ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_freq_plot.png", tr1_freq_plot, width =5.3333333, height=3, dpi=444, bg="white")


tr1_function_plot <- cell_count_data%>%
  filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", stim!="unstim", cell_pop!="Tr1_count")%>%
  mutate(cell_pop=case_when(cell_pop=="Tr1_IFNg_count"~ "IFNy",
                            cell_pop=="Tr1_IL10_count"~"IL10"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=cell_count*1000, fill=timepoint))+
  geom_point()+
  geom_line(aes(group=id))+
  geom_boxplot(outliers = F)+
  theme_minimal()+
  # scale_y_log10()+
  ylab("cyotkine+ Tr1 cells / mL")+
  facet_wrap(~infectiontype+cell_pop+stim, scales = "free", ncol=4)+
  scale_fill_manual(values=time_cols)+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size=7))

ggsave("~/postdoc/stanford/manuscripts/jason_tr1_2/tr1_cytokine.png", tr1_function_plot, width =8, height=5.333333, dpi=444, bg="white")

# sandbox####
# 
# count_nulisa_combo <- slim_nulisa_data%>%
#   mutate(id=as.integer(id))%>%
#   full_join(., cell_count_data, by = c("id", "timepoint", "infectiontype"))
# 
# 
# 
# grand_count_cor <- count_nulisa_combo%>%
#   filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in% c("A", "S"), !is.na(targetName))%>%
#   distinct(id, timepoint, infectiontype, absolute_Tr1, targetName, concentration)%>%
#   group_by(targetName)%>%
#   nest()%>%
#   mutate(correlation=map(data, ~cor.test(.$concentration, .$absolute_Tr1, method = "spearman")))%>%
#   mutate(p=map_dbl(correlation, ~.$p.value),
#          rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
#   ungroup()%>%
#   mutate(padj=p.adjust(p))
# 
# sig_count_cor <- grand_count_cor%>%
#   filter(p<0.05)
# 
# 
# grand_cor_plot <- count_nulisa_combo%>%
#   filter(gate=="Tr1_Frequency", stim=="unstim", infectiontype%in% c("A", "S"), targetName %in% sig_count_cor$targetName)%>%
#   distinct(id, timepoint, infectiontype, targetName, absolute_Tr1, concentration)%>%
#   ggplot(., aes(x=absolute_Tr1, y=concentration))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   ggpubr::stat_cor(method="spearman")+
#   facet_wrap(~targetName, scales = "free")+
#   theme_minimal()

# 
# long_combo%>%
#   filter(gate=="Tr1_IFNg_Frequency")%>%
#   distinct(id, infectiontype, timepoint, stim, freq)%>%
#   filter(infectiontype%in%c("A", "S"), timepoint!="bad_baseline")%>%
#   ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=freq, fill=timepoint))+
#   ggtitle("IFNg+ Tr1 Percentage")+
#   geom_boxplot()+
#   theme_minimal()+
#   facet_wrap(~infectiontype+stim, ncol=3)+
#   scale_fill_manual(values=viridis::magma(n=5))+
#   da_boxplot_theme
# 
# cell_count_data%>%
#   filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", cell_pop=="Tr1_IFNg_count")%>%
#   ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=cell_count*1000, fill=timepoint))+
#   geom_point()+
#   geom_line(aes(group=id))+
#   geom_boxplot(outliers = F)+
#   theme_minimal()+
#   ylab("cyotkine+ Tr1 cells / mL")+
#   ggtitle("IFNg+ Tr1 Count")+
#   facet_wrap(~infectiontype+stim, scales = "free", ncol=3)+
#   scale_fill_manual(values=viridis::magma(n=5))+
#   theme(axis.title.x = element_blank())
# 
# 
# 
# 
# 
# long_combo%>%
#   filter(gate=="Tr1_IL10_Frequency")%>%
#   distinct(id, infectiontype, timepoint, stim, freq)%>%
#   filter(infectiontype%in%c("A", "S"), timepoint!="bad_baseline")%>%
#   ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=freq, fill=timepoint))+
#   ggtitle("IL10+ Tr1 Percentage")+
#   geom_boxplot()+
#   theme_minimal()+
#   facet_wrap(~infectiontype+stim, ncol=3)+
#   scale_fill_manual(values=viridis::magma(n=5))+
#   da_boxplot_theme
# 
# cell_count_data%>%
#   filter(infectiontype %in% c("S", "A"), timepoint!="bad_baseline", cell_pop=="Tr1_IL10_count")%>%
#   ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), y=cell_count*1000, fill=timepoint))+
#   geom_point()+
#   geom_line(aes(group=id))+
#   geom_boxplot(outliers = F)+
#   theme_minimal()+
#   ylab("cyotkine+ Tr1 cells / mL")+
#   ggtitle("IL10+ Tr1 Count")+
#   facet_wrap(~infectiontype+stim, scales = "free", ncol=3)+
#   scale_fill_manual(values=viridis::magma(n=5))+
#   theme(axis.title.x = element_blank())
# 

# 
# na_counts <- clean_data_with_cell_counts%>%
#   filter(is.na(cbc_per_ul))
# 
# n_distinct(na_counts$sample_id)
# 
# na_counts_for_kenneth <- na_counts%>%
#   distinct(id, date, timepoint_imm)
# 
# write.csv(na_counts_for_kenneth,"~/Downloads/visits_without_cbc.csv", row.names = F)
# 
# good_counts <- clean_data_with_cell_counts%>%
#   filter(!is.na(cbc_per_ul))
# 
# n_distinct(good_counts$sample_id)