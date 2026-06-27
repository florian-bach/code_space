parasitemias <- read.csv("~/Downloads/boyle_chmi_parasitemias.csv")

clean_parasitemias <- parasitemias[,1:8]%>%
  mutate(timepoint=ifelse(grepl("^Pre", infectiontype), "baseline", timepoint))%>%
  mutate(infectiontype=gsub("Pre-", "", infectiontype))%>%
  mutate(infectiontype=gsub(" ", "", infectiontype, fixed = T))%>%
  mutate(timepoint=ifelse(timepoint=="baseline"&infectiontype=="IMP", "pre-treatment", timepoint))%>%
  mutate(parasites_per_500ul=gsub(",", "", parasites_per_500ul))%>%
  mutate(parasites_per_ml=ifelse(parasites_per_500ul=="ND", NA, as.numeric(parasites_per_500ul)*2))


pre_treatment_lines <- clean_parasitemias%>%
  filter(infectiontype %in% c("1st", "2nd"))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", paste(seq(4, 40, 0.5)))))%>%
  filter(timepoint %in% c("baseline", seq(4, 9, 0.5)))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, color=treatmentarm))+
  geom_line(aes(group=volunteer))+
  geom_point(aes(shape=treatmentarm))+
  scale_color_manual(values=c("blue", "red"))+
  scale_y_log10()+
  facet_wrap(~infectiontype)+
  theme_minimal(base_size = 20)

ggsave("~/postdoc/stanford/misc/rux_paper_figures/pre_treatment_lines.png", pre_treatment_lines, width=16, height=8, dpi=444, bg="white")



clean_parasitemias%>%
  filter(infectiontype %in% c("1st", "2nd"))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", paste(seq(4, 40, 0.5)))))%>%
  filter(timepoint %in% c("baseline", seq(4, 9, 0.5)))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, color=infectiontype))+
  geom_line(aes(group=volunteer))+
  geom_point(aes(shape=infectiontype))+
  scale_color_manual(values=c("darkgreen", "orange"))+
  scale_y_log10()+
  theme_minimal(base_size = 20)

clean_parasitemias%>%
  filter(infectiontype %in% c("2nd"))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", paste(seq(4, 40, 0.5)))))%>%
  filter(timepoint %in% c("baseline", seq(4, 10, 0.5)))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, fill=treatmentarm))+
  geom_boxplot()+
  geom_point(aes(color=treatmentarm), position = position_dodge(width=0.75))+
  scale_color_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000))+
  theme_minimal(base_size = 20)


clean_parasitemias%>%
  filter(infectiontype %in% c("1st", "2nd"))%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", paste(seq(4, 40, 0.5)))))%>%
  filter(timepoint %in% c("baseline", seq(4, 10, 0.5)))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, fill=infectiontype))+
  geom_boxplot()+
  geom_point(aes(color=infectiontype), position = position_dodge(width=0.75))+
  scale_color_manual(values=c("darkgreen", "orange"))+
  scale_fill_manual(values=c("darkgreen", "orange"))+
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000))+
  theme_minimal(base_size = 20)



# IMP ####
post_treatment_boxplots <- clean_parasitemias%>%
  filter(infectiontype %in% c("IMP"))%>%
  mutate(timepoint=factor(timepoint, levels=c("pre-treatment", paste(seq(2, 96, 2)))))%>%
  filter(timepoint %in% c("pre-treatment", paste(c(2, 4, 8, 12, 14, 16, 20, 24, 28, 36, 42, 48, 54, 60, 66, 72, 84, 96))))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, fill=treatmentarm))+
  geom_boxplot()+
  geom_point(aes(color=treatmentarm), position = position_dodge(width=0.75))+
  scale_color_manual(values=c("blue", "red"))+
  scale_fill_manual(values=c("blue", "red"))+
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000))+
  theme_minimal(base_size = 20)+
  theme(legend.position = "bottom")

ggsave("~/postdoc/stanford/misc/rux_paper_figures/post_treatment_boxplots.png", post_treatment_boxplots, width=16, height=8, dpi=444, bg="white")


post_treatment_line <- clean_parasitemias%>%
  filter(infectiontype %in% c("IMP"))%>%
  mutate(timepoint=factor(timepoint, levels=c("pre-treatment", paste(seq(2, 96, 2)))))%>%
  filter(timepoint %in% c("pre-treatment", paste(c(2, 4, 8, 12, 14, 16, 20, 24, 28, 36, 42, 48, 54, 60, 66, 72, 84, 96))))%>%
  ggplot(., aes(x=timepoint, y=parasites_per_ml, color=treatmentarm))+
  geom_line(aes(group=volunteer))+
  geom_point(aes(shape=treatmentarm))+
  scale_color_manual(values=c("blue", "red"))+
  scale_y_log10()+
  theme_minimal(base_size = 20)+
  theme(legend.position = "bottom")

ggsave("~/postdoc/stanford/misc/rux_paper_figures/post_treatment_line.png", post_treatment_line, width=16, height=8, dpi=444, bg="white")

