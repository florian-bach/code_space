library(WGCNA)
library(flashClust)
library(tidyr)
library(dplyr)
library(ggplot2)


protein_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

protein_data <- protein_data%>%
  filter(!is.na(new_qpcr))%>%
  mutate(day14_para=if_else(timepoint=="day14"&infectiontype=="A"& new_qpcr > 10, "parasitemic_day14", "no_parasites_day14"))%>%
  group_by(id)%>%
  mutate(class2= if_else(any(day14_para=="parasitemic_day14"), "more than 10 parasites / µL on day 14", "10 or fewer parasites / µL on day 14"))


baseline_qpcrs_plot <- protein_data%>%
  filter(timepoint=="baseline", infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, new_qpcr)%>%
  ggplot(., aes(x=class2, y=new_qpcr+0.001, fill=class2))+
  scale_y_log10(labels=scales::label_log())+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  ggpubr::stat_compare_means()+
  ylab("qPCR parasites / µL at baseline")+
  scale_fill_manual(values=c("darkgrey", "#636363"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/baseline_qpcrs_plot.png", baseline_qpcrs_plot, width = 5, height=4, dpi=444)

day0_qpcrs_plot <- protein_data%>%
  filter(timepoint=="day0", infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, new_qpcr)%>%
  ggplot(., aes(x=class2, y=new_qpcr+0.001, fill=class2))+
  scale_y_log10(labels=scales::label_log())+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  ggpubr::stat_compare_means()+
  ylab("qPCR parasites / μL at day0")+
  scale_fill_manual(values=c("darkgrey", "#636363"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/day0_qpcrs_plot.png", day0_qpcrs_plot, width = 5, height=4, dpi=444)



age_plot <- protein_data%>%
  filter(timepoint=="baseline", infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, ageyrs)%>%
  ggplot(., aes(x=class2, y=ageyrs, fill=class2))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  ggpubr::stat_compare_means()+
  ylab("age in years")+
  scale_fill_manual(values=c("darkgrey", "#636363"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/age_plot.png", age_plot, width = 5, height=4, dpi=444)


gender_plot <- protein_data%>%
  filter(timepoint=="baseline", infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, gender_categorical)%>%
  group_by(class2)%>%
  summarise("percentage_male"=sum(gender_categorical=="Male")/n())%>%
  ggplot(., aes(x=class2, y=percentage_male, fill=class2))+
  geom_bar(stat="identity")+
  ggpubr::stat_compare_means()+
  ylab("percentage male")+
  scale_fill_manual(values=c("darkgrey", "#636363"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  scale_y_continuous(limits=c(0,1), labels = scales::label_percent())+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/gender_plot.png", gender_plot, width = 5, height=4, dpi=444)


# infection duration ####

durations <- readxl::read_xls("~/postdoc/stanford/rna_seq/MUSICAL/Duration_for_Jason_Rev.xls")
colnames(durations)[c(1,2)] <- c("id", "date")

protein_data_with_duration <- durations%>%
  select(id, date, duration, durationcat)%>%
  mutate(date=as.character(date))%>%
  right_join(., protein_data, by=c("id", "date"))


duration_plot1 <- protein_data_with_duration%>%
  filter(infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, duration)%>%
  # group_by(class2)%>%
  # summarise("percentage_male"=sum(gender_categorical=="Male")/n())%>%
  ggplot(., aes(x=class2, y=duration, fill=class2))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.2))+
  ggpubr::stat_compare_means()+
  ylab("duration of infection in days")+
  scale_fill_manual(values=c("darkgrey", "#636363"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  # scale_y_continuous(limits=c(0,1), labels = scales::label_percent())+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/duration_plot1.png", duration_plot1, width = 5, height=4, dpi=444)


duration_plot2 <- protein_data_with_duration%>%
  filter(infectiontype=="A", !is.na(class2))%>%
  distinct(id, class2, duration, durationcat)%>%
  mutate(durationcat=ifelse(is.na(durationcat)&!is.na(duration)&duration<=14, "<14", durationcat))%>%
  group_by(durationcat)%>%
  summarise("percentage_control"=sum(class2=="10 or fewer parasites / µL on day 14")/n())%>%
  mutate(durationcat=factor(durationcat, levels=c("<14", ">14-<50",">=50-<90",">=90")))%>%
  ggplot(., aes(x=durationcat, y=percentage_control, fill=durationcat))+
  geom_bar(stat="identity")+
  ylab("percentage of controllers in duration stratum")+
  scale_fill_manual(values=colorspace::sequential_hcl(palette = "blues", n=5))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  scale_y_continuous(limits=c(0,1), labels = scales::label_percent())+
  theme_minimal(base_size = 12)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/stanford/rna_seq/MUSICAL/figures/duration_plot2.png", duration_plot2, width = 5, height=4, dpi=444)

