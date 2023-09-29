library(dplyr)
library(tidyr)
library(ggplot2)

comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/MICDROP all visit database through June 30th 2023.dta")
mic_drop <- mic_drop %>%
  mutate(mstatus = case_match(mstatus,
                              0~"asymptomatic",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"severe"))%>%
  mutate("any_parsdens" = qPCRparsdens)%>%
  mutate("parasitaemia_method" = if_else(qPCRdich==1, "qPCR", if_else(BSdich==1, "smear", "dunno")))%>%
  mutate(any_parsdens = if_else(is.na(qPCRparsdens) & !is.na(pardens), pardens, any_parsdens))%>%
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))
  

data_labs <- lapply(mic_drop, function(x) attributes(x)$labels)
listy_data_labs <- lapply(names(data_labs), function(x)ifelse(is.null(data_labs[[x]]), "n/a", data_labs[x]))
names(listy_data_labs) <- names(data_labs)


# view parasitaemia over time as order of infections ####
smear_positive <- mic_drop %>%
  filter(BSdich==1)%>%
  mutate(id=factor(id))%>%
  group_by(id)%>%
  add_count(name="total_positive_smears") %>%
  arrange(compage) %>%
  mutate(n_smear = seq(1, max(total_positive_smears)))

n_infection_cols <- c("white", colorspace::sequential_hcl(n=max(smear_positive$n_smear), palette = "Lajolla"))

pcr_pos <- mic_drop %>%
  filter(qPCRdich==1)%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_positive_qpcrs") %>%
  arrange(compage) %>%
  mutate(n_smear = seq(1, max(total_positive_qpcrs)))



# smear positive events n_infection
micro_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=pardens, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])


# pcr_positive evens n_infection
pcr_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=qPCRparsdens, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])

dens_n_smear_combo <- cowplot::plot_grid(micro_dens, pcr_dens,  nrow = 1)





# smear positive events by age
micro_age <- ggplot(smear_positive, aes(x=AGE, y=pardens, fill=factor(AGE)))+
  # geom_point(alpha=0.2)+
  # geom_line(alpha=0.2, aes(group=id))+
  geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fatten=0, fun.min = median, fun.max = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Age (Weeks)")+
  geom_smooth()+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE)


# pcr positive events by age
pcr_age <- ggplot(smear_positive, aes(x=AGE, y=qPCRparsdens+0.001, fill=factor(AGE)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  # geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Age (Weeks)")+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE)

dens_age_combo <- cowplot::plot_grid(micro_age, pcr_age,  nrow = 1)

# comp cases ####

comp_cases <- mic_drop %>%
  filter(mstatus %in% c("complicated", "severe"))

kids_with_comp <- comp_cases %>%
  reframe("id"=id)


# bs_pos_visits_of_kids_with_comp <- mic_drop %>%
#   filter(id %in% kids_with_comp$id, BSdich==1) %>%
#   group_by(id) %>%
#   add_count(name="total_positive_smears") %>%
#   arrange(compage) %>%
#   mutate(n_smear = seq(1, max(total_positive_smears)))


any_pos_visits_of_kids_with_comp <- mic_drop %>%
  filter(id %in% kids_with_comp$id, parasitaemia_method %in% c("qPCR", "smear")) %>%
  group_by(id) %>%
  add_count(name="total_positive_measurement") %>%
  arrange(compage) %>%
  mutate(n_para = seq(1, max(total_positive_measurement)))

# 
# ggplot(any_pos_visits_of_kids_with_comp, aes(x=factor(n_para), y=as.numeric(pardens)))+
#   geom_point(aes(color=factor(mstatus)))+
#   geom_line(alpha=0.3, aes(group=id))+
#   ylab("Microscopic Parasite Density")+
#   xlab("Order of Positive Measurement")+
#   scale_y_log10()+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=n_infection_cols[-1])+
#   scale_color_manual(values=comp_pal)
# 
# 
# ggplot(any_pos_visits_of_kids_with_comp, aes(x=, y=as.numeric(pardens)))+
#   geom_point(aes(color=factor(mstatus)))+
#   geom_line(alpha=0.3, aes(group=id))+
#   # geom_boxplot(aes(group=factor(n_smear), fill=factor(n_smear)))+
#   # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#   #              geom = "crossbar", width = 0.8, color="white")+
#   ylab("Microscopic Parasite Density")+
#   xlab("Age (months)")+
#   scale_y_log10()+
#   scale_x_continuous(breaks=seq(0,15,by=3))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_manual(values=n_infection_cols[-1])+
#   scale_color_manual(values=comp_pal)


# make a data frame where we gather the visits where complicated / severe occurred, and the previous
# parasitaemic event

n_inf_df <- any_pos_visits_of_kids_with_comp %>%
  select(id, AGE, n_para, mstatus)%>%
  group_by(id)%>%
  filter(mstatus %in% c("complicated", "severe"))%>%
  reframe(id, n_para)

label_df <- data.frame()

for(i in 1:nrow(n_inf_df)){
  #visit where complicated / severe malaria was diagnosed
  comp_visit <- n_inf_df[i,]
  # previous visit where treatment was administered (uncomplicated malaria was diagnosed)
  previous_treatment_visits <- any_pos_visits_of_kids_with_comp%>%
    filter(id==comp_visit$id,
           mstatus == "uncomplicated",
           n_para < comp_visit$n_para)%>%
    slice_max(n_para)
  
  most_recent_treament_visit <- previous_treatment_visits$n_para
  # label_visits <- filter(any_pos_visits_of_kids_with_comp, id==subject & n_para %in% c(comp_visit, comp_visit-1))
  label_visits <- filter(any_pos_visits_of_kids_with_comp, id==comp_visit$id & n_para %in% c(comp_visit$n_para, most_recent_treament_visit))
  label_df <- rbind(label_df, label_visits)
}

label_df$age_in_days <- label_df$date-label_df$dob

indie_comp_n_infection <- ggplot(any_pos_visits_of_kids_with_comp, aes(x=factor(n_para), y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(shape=parasitaemia_method,
    color=factor(mstatus, levels=c("asymptomatic",
                                                "uncomplicated",
                                                "complicated",
                                                "severe"))))+
  geom_line(alpha=0.3, aes(group=id))+
  facet_wrap(~ id, scales = "free_x")+
  ylab("Microscopic Parasite Density")+
  xlab("order of parasitaemic event")+
  scale_y_log10()+
  #scale_x_continuous(breaks=seq(0,15,by=3))+
  theme_minimal()+
  scale_shape_manual(values=c(16,15))+
  scale_color_manual(values=comp_pal)+
  guides(color=guide_legend(title="disease"))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/individual_complicated_histories.png", indie_comp_bs, height=6, width=8, dpi=444, bg="white")




(indie_comp_age <- 
    ggplot(any_pos_visits_of_kids_with_comp, aes(x=AGE, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(shape=parasitaemia_method,
                 color=factor(mstatus, levels=c("asymptomatic",
                                                "uncomplicated",
                                                "complicated",
                                                "severe"))))+
  geom_line(alpha=0.3, aes(group=id))+
  ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  ylab("qPCR Parasite Density")+
  xlab("Age (weeks)")+
  scale_y_log10()+
  scale_x_continuous(breaks=seq(0,72,by=4),
                     # labels = seq(0,72,by=8)
                     )+
  theme_minimal()+
  scale_shape_manual(values=c(16,15))+
  scale_fill_manual(values=n_infection_cols[-1])+
  scale_color_manual(values=comp_pal)+
  guides(color=guide_legend(title="disease")))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/individual_complicated_histories.png", indie_comp_age, height=10, width=18, dpi=444, bg="white")


# determine treatment failure ####

# keep kids with 1 or more malaria episodes (they are only treated when diagnosed);
# look at all visits where any parasites where measured, calculate the time interval between
more_than_one <- mic_drop %>%
  group_by(id) %>%
  mutate("number_of_malaria_episodes" = sum(mstatus=="uncomplicated"))%>%
  ungroup()%>%
  select(id, date, AGE, dob, mstatus, any_parsdens, parasitaemia_method, number_of_malaria_episodes)%>%
  filter(number_of_malaria_episodes>=1)%>%
  filter(any_parsdens>=1 & !is.na(any_parsdens)) %>%
  group_by(id) %>%
  add_count(name="total_positive_measurement")%>%
  mutate(lagdate = lag(date), "time_to_previous" = date - lagdate, "time_to_next"=lead(time_to_previous))
 
failure_to_clear <- more_than_one %>%
  filter(mstatus=="uncomplicated", time_to_next<=13)

# failure_to_cure <- more_than_one %>%
#   filter(mstatus=="uncomplicated", time_to_previous<=13)



treatment_failure_plot <- more_than_one %>%
  filter(id %in% unique(failure_to_clear$id))%>%
  ggplot(., aes(x=AGE, y=as.numeric(any_parsdens)+0.001))+
    geom_point(aes(shape=parasitaemia_method,
                   color=factor(mstatus, levels=c("asymptomatic",
                                                  "uncomplicated",
                                                  "complicated",
                                                  "severe"))))+
    geom_line(alpha=0.3, aes(group=id))+
    # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
    facet_wrap(~ id)+
    ylab("qPCR Parasite Density")+
    xlab("Age (weeks)")+
    scale_y_log10()+
    scale_x_continuous(breaks=seq(0,72,by=4),
                       # labels = seq(0,72,by=8)
    )+
    theme_minimal()+
    # scale_shape_manual(values=c(16,15))+
    scale_fill_manual(values=n_infection_cols[-1])+
    scale_color_manual(values=comp_pal)+
    guides(color=guide_legend(title="disease"))+
    theme(axis.text.x = element_text(size=5))


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/treatment_failures.png", treatment_failure_plot, width = 16, height=8, bg="white")

# add the shortest lag between treatment and parasitaemic event
# 
# label_df <- data.frame()
# 
# for(i in 1:nrow(n_inf_df)){
#   #visit where complicated / severe malaria was diagnosed
#   comp_visit <- n_inf_df[i,]
#   # previous visit where treatment was administered (uncomplicated malaria was diagnosed)
#   previous_treatment_visits <- any_pos_visits_of_kids_with_comp%>%
#     filter(id==comp_visit$id,
#            mstatus == "uncomplicated",
#            n_para < comp_visit$n_para)%>%
#     slice_max(n_para)
#   
#   most_recent_treament_visit <- previous_treatment_visits$n_para
#   # label_visits <- filter(any_pos_visits_of_kids_with_comp, id==subject & n_para %in% c(comp_visit, comp_visit-1))
#   label_visits <- filter(any_pos_visits_of_kids_with_comp, id==comp_visit$id & n_para %in% c(comp_visit$n_para, most_recent_treament_visit))
#   label_df <- rbind(label_df, label_visits)
# }
# 
# label_df$age_in_days <- label_df$date-label_df$dob
# 






# complicated syndromes ####
long_comp_cases <- comp_cases %>%
  select(id, date, AGE, mstatus, any_parsdens, conv, sitstand, vomit, lethargy, smal, unbreast, cmal, mseiz, rdist, othcriteria, hosp, fever, hb, plt, wbc, neutro, lymph, mono, eosino,hbgrade, pltgrade, wbcgrade, neutrograde, lymphgrade, monograde, eosinograde)%>%
  pivot_longer(cols=c(conv, sitstand, vomit, lethargy, smal, unbreast, cmal, mseiz, rdist, othcriteria), names_to = "symptom", values_to = "symptom_value", values_transform = as.character)%>%
  pivot_longer(cols=c(hb, plt, wbc, neutro, lymph, mono, eosino), names_to = "cell", values_to = "count")%>%
  pivot_longer(cols=c(hbgrade, pltgrade, wbcgrade, neutrograde, lymphgrade, monograde, eosinograde),values_to = "haematology_grade")
