library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
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

# treatment failures ####

tf <- c(10680, 10731, 10839, 10748, 11058, 11069, 11414)



actual_treatment_failure_plot <- more_than_one %>%
  filter(id %in% tf)%>%
  # mutate(id_dob=paste(id, dob))
  ggplot(., aes(x=AGE, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(shape=parasitaemia_method,
                 color=factor(mstatus, levels=c("asymptomatic",
                                                "uncomplicated",
                                                "complicated",
                                                "severe"))))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id,nrow = 4)+
  ylab("qPCR Parasite Density")+
  xlab("Age (weeks)")+
  scale_y_log10()+
  scale_x_continuous(breaks=seq(0,72,by=2),
                     # labels = seq(0,72,by=8)
  )+
  theme_minimal()+
  # scale_shape_manual(values=c(16,15))+
  scale_fill_manual(values=n_infection_cols[-1])+
  scale_color_manual(values=comp_pal)+
  guides(color=guide_legend(title="disease"))+
  theme(axis.text.x = element_text(size=5))


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/actual_treatment_failures.png", actual_treatment_failure_plot, width = 8, height=8, bg="white")


kids_with_tf <- data.frame("id"=c(rep(10680, 2), 10731, 10748, 10839, 11058, 11069, rep(11414, 2)),
                           "age_in_wk_at_tf"=c(16, 36, 55, 62, 46, 39, 46, 28, 30))

kids_with_tf$dob <- more_than_one$dob[match(kids_with_tf$id, more_than_one$id)]

tfs <- more_than_one %>%
  filter(id %in% kids_with_tf$id &AGE %in% kids_with_tf$age_in_wk_at_tf & mstatus!="asymptomatic")%>%
  select(date)

#there's a stray one
tfs <- tfs[-4,]

kids_with_tf$dotf <- tfs$date

actual_treatment_failure_timeline <- more_than_one %>%
  filter(id %in% tf)%>%
  # mutate(id_dob=paste(id, dob))
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(shape=parasitaemia_method,
                 color=factor(mstatus, levels=c("asymptomatic",
                                                "uncomplicated",
                                                "complicated",
                                                "severe"))))+
  geom_line(alpha=0.3, aes(group=id))+
  geom_vline(xintercept = kids_with_tf$dotf, linetype="dotted")+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  ylab("qPCR Parasite Density")+
  scale_x_date(date_breaks="2 months")+
  scale_y_log10()+
  theme_minimal()+
  scale_fill_manual(values=n_infection_cols[-1])+
  scale_color_manual(values=comp_pal)+
  guides(color=guide_legend(title=""),
         shape=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box="vertical",
        legend.margin=margin())


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/actual_treatment_failure_timeline.png", actual_treatment_failure_timeline, width = 8, height=4, bg="white")
# complicated syndromes ####

long_comp_cases <- comp_cases %>%
  select(id, date, AGE, mstatus, any_parsdens, conv, sitstand, vomit, lethargy, smal, unbreast, cmal, mseiz, rdist, othcriteria, hosp, fever, hb, plt, wbc, neutro, lymph, mono, eosino,hbgrade, pltgrade, wbcgrade, neutrograde, lymphgrade, monograde, eosinograde)%>%
  pivot_longer(cols=c(conv, sitstand, vomit, lethargy, smal, unbreast, cmal, mseiz, rdist, othcriteria), names_to = "symptom", values_to = "symptom_value", values_transform = as.character)%>%
  pivot_longer(cols=c(hb, plt, wbc, neutro, lymph, mono, eosino), names_to = "cell", values_to = "count")%>%
  pivot_longer(cols=c(hbgrade, pltgrade, wbcgrade, neutrograde, lymphgrade, monograde, eosinograde),values_to = "haematology_grade")
# more plots to be organised ####


big_plot <- mic_drop %>%
  # filter(!is.na(parasitaemia_method), !is.na(mstatus), any_parsdens != 0)%>%
  filter(!is.na(parasitaemia_method), !is.na(mstatus))%>%
  # filter(mstatus=="asymptomatic")%>%
  ggplot(aes(x=AGE, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(shape=parasitaemia_method,
                 fill=factor(mstatus, levels=c("asymptomatic",
                                               "uncomplicated",
                                               "complicated",
                                               "severe"))))+
  geom_line(alpha=0.3, aes(group=id))+
  ylab("Microscopic Parasite Density")+
  xlab("order of parasitaemic event")+
  scale_y_log10(limits=c(0.001, NA))+
  facet_wrap(~id)+
  #scale_x_continuous(breaks=seq(0,15,by=3))+
  theme_minimal()+
  scale_shape_manual(values=c(21,22,24))+
  scale_fill_manual(values=comp_pal)+
  guides(fill=guide_legend(title="disease"))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/big_parsdens.png", big_plot, width = 40, height = 25, bg="white")





#plot the evolution of complicated malaria risk as a function of the order of infection

comp_model <- glm(comp_severe/uncomplicated~n_infection, weights = uncomplicated+comp_severe, family="binomial", data=clinical_summary_df)
comp_model_fun <- function(x){exp(comp_model$coefficients[1])*exp(comp_model$coefficients[2])^x}


# clinical_summary_df <- all_kids_with_clinical_episodes %>%
#   group_by(n_infection, mstatus)%>%
#   count()%>%
#   pivot_wider(names_from = mstatus, values_from = n)%>%
#   mutate(comp_severe=sum(complicated, severe, na.rm = TRUE))%>%
#   mutate(comp_severe=if_else(is.na(comp_severe), 0, comp_severe))

clinical_summary_df <- all_kids_with_clinical_episodes %>%
  mutate("age_in_months"=round(AGE/4.33))%>%
  group_by(age_in_months, mstatus)%>%
  count()%>%
  pivot_wider(names_from = mstatus, values_from = n)%>%
  mutate(comp_severe=sum(complicated, severe, na.rm = TRUE))%>%
  mutate(comp_severe=if_else(is.na(comp_severe), 0, comp_severe))


ggplot(clinical_summary_df, aes(x=age_in_months, y=comp_severe/uncomplicated))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data= model_visualiser(comp_model, "n_infection"), aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun=comp_model_fun, colour="black")+
  scale_y_continuous(limits=c(0,0.1), labels=scales::label_percent())+
  scale_x_continuous(limits=c(0,19), breaks = seq(1,10,1))+
  ylab("fraction of complicated cases")+
  theme(legend.position = "none")

# plot parasitaemia patterns by binning parasitaemic months ####
binned_parasite_densities <- mic_drop %>%
  # mutate("age_in_months"=round(AGE/4.33))%>%
  filter(!is.na(pardens), dob < "2022-04-30", AGE<=52, AGE>=8)%>%
  mutate("binned_age"=case_when(AGE <= 18 ~ "8-18 weeks",
                                AGE > 18 & AGE <= 28 ~ "19-29 weeks",
                                AGE > 28 & AGE <= 40 ~ "30-40 weeks",
                                AGE > 40 & AGE <= 52 ~ "41-52 weeks"))%>%
  group_by(id, binned_age)%>%
  # summarise("geo_parasites"=exp(mean(log(pardens+0.0001))))%>%
  summarise("geo_parasites"=exp(mean(log(pardens+1))))%>%
  pivot_wider(names_from = binned_age, values_from = geo_parasites)%>%
  select(id, "8-18 weeks", "19-29 weeks","30-40 weeks", "41-52 weeks")


binned_parasite_densities <- replace(binned_parasite_densities, is.na(binned_parasite_densities), 0)
log_binned_parasite_densities <- binned_parasite_densities
# log_binned_parasite_densities[,2:5] <- log10(binned_parasite_densities[,2:5]+0.00001)
log_binned_parasite_densities[,2:5] <- log10(binned_parasite_densities[,2:5]+1)



no_zero <- as.matrix(log_binned_parasite_densities[,2:5])
baseline_dist <- dist(no_zero, method = "euclidean")

baseline_hclust <- hclust(baseline_dist)

# cut_tree <- cutree(baseline_hclust, k=6) 
# colnames(no_zero) <- c("2-3 months", "4-5 months", "6-7 months", "8-9 months")

inferno <- colorspace::sequential_hcl("inferno", n=9)
viridis <- colorspace::sequential_hcl("viridis", n=9)
# col_inferno <- circlize::colorRamp2(c(-5, -1.5, -1, -0.5, 0, 1, 2, 3, 4), inferno)
# col_inferno2 <- circlize::colorRamp2(c(seq(-0.5, 6, by=7/9)), colorspace::sequential_hcl("viridis", n=9))

disease_phenotypes <- mic_drop %>%
  group_by(id)%>%
  summarise("disease"=if_else(
    all(mstatus == "asymptomatic"), "asymptomatic", if_else(
      any(mstatus == "severe"), "severe", if_else(
        any(mstatus == "complicated"), "complicated", if_else(
          any(mstatus=="uncomplicated"), "uncomplicated", "dont_know")
      )
    )
  )
  )

max_parasites <- mic_drop %>%
  group_by(id)%>%
  summarise("max_parasites"=log10(max(any_parsdens, na.rm = TRUE)+0.001))


col_inferno <- circlize::colorRamp2(seq(min(no_zero), max(no_zero), length.out=length(inferno)), inferno)
col_inferno2 <- circlize::colorRamp2(seq(min(max_parasites$max_parasites), max(max_parasites$max_parasites), length.out=length(viridis)), viridis)

disease_phenotypes <- disease_phenotypes[match(log_binned_parasite_densities$id, disease_phenotypes$id),]
max_parasites <- max_parasites[match(log_binned_parasite_densities$id, max_parasites$id),]

right_anno <-  rowAnnotation(#annotation_name_gp = gpar(fontsize=10),
  # annotation_name_rot = 45,
  gap = unit(1.5, "mm"),
  "disease phenotype"= disease_phenotypes$disease,
  "max parasitaemia"=max_parasites$max_parasites,
  show_legend = TRUE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm"), # width of the significance bar
  col=list("disease phenotype" = c("asymptomatic"="darkgrey",
                                   "uncomplicated"="#ffda03",
                                   "complicated"="darkorange",
                                   "severe"="darkred"),
           "max parasitaemia"=col_inferno2),
  annotation_legend_param = list("disease phenotype" = list(title = "disease phenotype",
                                                            at = c("asymptomatic", "uncomplicated", "complicated", "severe"),
                                                            title_gp=gpar(angle=45),
                                                            legend_gp = gpar(fill = c("darkgrey","#ffda03", "darkorange", "darkred")),
                                                            title_position = "topleft"),
                                 "max parasitaemia" = list(title = "highest log10\nparasitaemia")
  )
  
)


internal_split_plot <- Heatmap(matrix = no_zero,
        right_annotation=right_anno,
        cluster_rows = baseline_hclust,
        # cluster_rows=TRUE,
        cluster_columns=FALSE,
        show_row_dend = TRUE,
        cluster_row_slices = TRUE,
        show_heatmap_legend = TRUE,
        # name = "geometric mean log10 parasites in time window",
        #cluster_columns = FALSE,
        row_split=8,
        column_names_gp = gpar(fontsize = 6),
        column_names_centered=TRUE,
        heatmap_legend_param = list(title = "geometric mean log10\nparasitaemia in time window", title_position = "topleft", grid_width=unit(7, "mm")),
        row_names_gp = gpar(fontsize = 6),
        row_names_side = "left",
        col = col_inferno,
        column_names_rot = 0
)

png("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_06/figures/internal_split_plot_euclidean.png", width = 6, height=12, units = "in", res = 444)
draw(internal_split_plot)
dev.off()


# sandbox ####
# 
# asymp <- mic_drop %>%
#   filter(!is.na(parasitaemia_method), !is.na(mstatus), any_parsdens != 0)%>%
#   filter(mstatus=="asymptomatic")%>%
#   ggplot(aes(x=date, y=as.numeric(any_parsdens)+0.001))+
#   geom_point(aes(shape=parasitaemia_method,
#                  fill=factor(mstatus, levels=c("asymptomatic",
#                                                "uncomplicated",
#                                                "complicated",
#                                                "severe"))))+
#   # geom_line(alpha=0.3, aes(group=id))+
#   ylab("Microscopic Parasite Density")+
#   xlab("order of parasitaemic event")+
#   scale_y_log10(limits=c(0.001, NA))+
#   # facet_wrap(~mstatus)+
#   #scale_x_continuous(breaks=seq(0,15,by=3))+
#   theme_minimal()+
#   scale_shape_manual(values=c(21,22,24))+
#   scale_fill_manual(values=comp_pal)+
#   guides(fill=guide_legend(title="disease"))
# 
# uncomp <- mic_drop %>%
#   filter(!is.na(parasitaemia_method), !is.na(mstatus), any_parsdens != 0)%>%
#   filter(mstatus=="uncomplicated")%>%
#   ggplot(aes(x=date, y=as.numeric(any_parsdens)+0.001))+
#   geom_point(aes(shape=parasitaemia_method,
#                  fill=factor(mstatus, levels=c("asymptomatic",
#                                                "uncomplicated",
#                                                "complicated",
#                                                "severe"))))+
#   # geom_line(alpha=0.3, aes(group=id))+
#   ylab("Microscopic Parasite Density")+
#   xlab("order of parasitaemic event")+
#   scale_y_log10(limits=c(0.001, NA))+
#   # facet_wrap(~mstatus)+
#   #scale_x_continuous(breaks=seq(0,15,by=3))+
#   theme_minimal()+
#   scale_shape_manual(values=c(21,22,24))+
#   scale_fill_manual(values=comp_pal)+
#   guides(fill=guide_legend(title="disease"))
# 
# comp <- mic_drop %>%
#   filter(!is.na(parasitaemia_method), !is.na(mstatus), any_parsdens != 0)%>%
#   filter(mstatus=="complicated")%>%
#   ggplot(aes(x=date, y=as.numeric(any_parsdens)+0.001))+
#   geom_point(aes(shape=parasitaemia_method,
#                  fill=factor(mstatus, levels=c("asymptomatic",
#                                                "uncomplicated",
#                                                "complicated",
#                                                "severe"))))+
#   # geom_line(alpha=0.3, aes(group=id))+
#   ylab("Microscopic Parasite Density")+
#   xlab("order of parasitaemic event")+
#   scale_y_log10(limits=c(0.001, NA))+
#   # facet_wrap(~mstatus)+
#   #scale_x_continuous(breaks=seq(0,15,by=3))+
#   theme_minimal()+
#   scale_shape_manual(values=c(21,22,24))+
#   scale_fill_manual(values=comp_pal)+
#   guides(fill=guide_legend(title="disease"))
# 
# 
# 
# severe <- mic_drop %>%
#   filter(!is.na(parasitaemia_method), !is.na(mstatus), any_parsdens != 0)%>%
#   filter(mstatus=="severe")%>%
#   ggplot(aes(x=date, y=as.numeric(any_parsdens)+0.001))+
#   geom_point(aes(shape=parasitaemia_method,
#                  fill=factor(mstatus, levels=c("asymptomatic",
#                                                "uncomplicated",
#                                                "complicated",
#                                                "severe"))))+
#   # geom_line(alpha=0.3, aes(group=id))+
#   ylab("Microscopic Parasite Density")+
#   xlab("order of parasitaemic event")+
#   scale_y_log10(limits=c(0.001, NA))+
#   # facet_wrap(~mstatus)+
#   #scale_x_continuous(breaks=seq(0,15,by=3))+
#   theme_minimal()+
#   scale_shape_manual(values=c(21,22,24))+
#   scale_fill_manual(values=comp_pal)+
#   guides(fill=guide_legend(title="disease"))
# 
# 
# ggExtra::ggMarginal(asymp, type = "histogram", groupFill = TRUE)
# ggExtra::ggMarginal(uncomp, type = "histogram", groupFill = TRUE)
# ggExtra::ggMarginal(comp, type = "histogram", groupFill = TRUE)
# ggExtra::ggMarginal(severe, type = "histogram", groupFill = TRUE)


# 
# bds <- mic_drop %>%
#   group_by(id)%>%
#   summarise("dob"=unique(dob))
# 
# ggplot(bds, aes(x=dob))+
#   geom_histogram(binwidth = 7)+
#   scale_x_date(breaks="1 week")+
#   theme_minimal()+
  # theme(axis.text.x = element_text(angle=90))
# disease_phenotypes <- mic_drop %>%
#   group_by(id)%>%
#   summarise("disease"=if_else(
#     all(mstatus == "asymptomatic"), "asymptomatic", if_else(
#       any(mstatus == "severe"), "severe", if_else(
#         any(mstatus == "complicated"), "complicated", if_else(
#           any(mstatus=="uncomplicated"), "uncomplicated", "dont_know")
#       )
#     )
#   )
#   )
# 
# disease_phenotypes <- disease_phenotypes[match(log_binned_parasite_densities$id, disease_phenotypes$id),]
# 
# right_anno <-  rowAnnotation(#annotation_name_gp = gpar(fontsize=10),
#   # annotation_name_rot = 45,
#   gap = unit(1.5, "mm"),
#   "disease phenotype"= disease_phenotypes$disease,
#   show_legend = TRUE,
#   show_annotation_name = FALSE,
#   simple_anno_size = unit(3, "mm"), # width of the significance bar
#   col=list("disease phenotype" = c("asymptomatic"="darkgrey",
#                                    "uncomplicated"="#ffda03",
#                                    "complicated"="darkorange",
#                                    "severe"="darkred")),
#   annotation_legend_param = list("disease phenotype" = list(title ="disease phenotype",
#                                                             at = c("asymptomatic", "uncomplicated", "complicated", "severe"),
#                                                             title_gp=gpar(angle=45),
#                                                             legend_gp = gpar(fill = c("darkgrey","#ffda03", "darkorange", "darkred")),
#                                                             title_position = "topleft")
#   )
#   
# )
# 
# 
# Heatmap(matrix = no_zero,
#         right_annotation=right_anno,
#         
#         cluster_rows = baseline_hclust,
#         # cluster_rows=FALSE,
#         cluster_columns=FALSE,
#         show_row_dend = TRUE,
#         cluster_row_slices = FALSE,
#         show_heatmap_legend = TRUE,
#         # name = "geometric mean log10 parasites in time window",
#         #cluster_columns = FALSE,
#         row_split=8,
#         column_names_gp = gpar(fontsize = 6),
#         column_names_centered=TRUE,
#         heatmap_legend_param = list(title = "geometric mean log10\nparasitaemia in time window", title_position = "topleft", grid_width=unit(7, "mm")),
#         row_names_gp = gpar(fontsize = 6),
#         row_names_side = "left",
#         col = col_inferno,
#         column_names_rot = 0
# )
# 
# 
# bds <- mic_drop %>%
#   group_by(id)%>%
#   summarise("dob"=unique(dob))
# 
# ggplot(bds, aes(x=dob))+
#   geom_histogram(binwidth = 7)+
#   scale_x_date(breaks="1 week")+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle=90))
# 
# 
# 
# 
# 
# 
# complicated_data%>%
#   arrange(desc(disease))%>%
#   ggplot(aes(x=factor(n_infection), y=parsdens, fill=factor(n_infection)))+
#   geom_boxplot()+
#   # geom_point(color="darkred", aes(alpha=factor(disease)))+
#   geom_point(aes(color=factor(disease)))+
#   theme_minimal()+
#   ylab("Parasites /  μL")+
#   scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
#   scale_color_manual(values = c("darkred", "orange"))+
#   scale_alpha_manual(values = c("uncomplicated"=0.1, "complicated"=1))+
#   scale_y_log10()+
#   theme(legend.position = "none")
# xlab("Order of Infection")
# 
# complicated_data %>%
#   group_by(id)%>%
#   mutate("high_before_six"=if_else(any(parsdens>10000) & age>=0.5, "yes", "no"))%>%
#   arrange(desc(disease))%>%
#   ggplot(aes(x=factor(n_infection), y=parsdens, fill=factor(n_infection)))+
#   geom_boxplot()+
#   # geom_point(color="darkred", aes(alpha=factor(disease)))+
#   geom_point(aes(color=factor(disease)))+
#   facet_wrap(~high_before_six)+
#   theme_minimal()+
#   ylab("Parasites /  μL")+
#   scale_fill_manual(values = colorspace::sequential_hcl(10, palette = "Purple Yellow"))+
#   scale_color_manual(values = c("darkred", "orange"))+
#   scale_alpha_manual(values = c("uncomplicated"=0.1, "complicated"=1))+
#   scale_y_log10()+
#   theme(legend.position = "none")
# xlab("Order of Infection")
# 
# 
# complicated_data <- complicated_data %>%
#   group_by(id)%>%
#   mutate("high_before_six"=if_else(any(parsdens>10000) & age <= 0.5, "yes", "no"))
# 
# df <- complicated_data %>%
#   group_by(high_before_six, n_infection)%>%
#   count(disease)%>%
#   pivot_wider(values_from = n, names_from = disease)%>%
#   mutate(complicated=replace_na(complicated, 0))%>%
#   mutate(risk=complicated/(uncomplicated+complicated))
# 
# 
# 
# ggplot(df, aes(x=n_infection, y=risk, fill=high_before_six))+
#   geom_bar(stat="identity")+
#   geom_text(aes(label= paste0("frac(",complicated, ",", uncomplicated+complicated,")")),parse = TRUE, vjust= -0.2, size=3.5)+
#   facet_wrap(~high_before_six)+
#   scale_x_continuous(breaks = seq(1, 7), limits=c(0, 7))+
#   scale_y_continuous(labels=scales::label_percent(), limits = c(0,0.15))+
#   scale_fill_manual(values=c("magenta4", "darkred"))+
#   theme_minimal()
# 





