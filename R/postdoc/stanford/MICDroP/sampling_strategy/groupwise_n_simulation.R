library(dplyr)
library(ggplot2)
library(patchwork)


comp_pal <- c("asymptomatic"="lightgrey",
              "uncomplicated"="black",
              "complicated"="orange",
              "severe"="darkred")


mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_04/MICDROP expanded database through April 30th 2024.dta")


# n_infection ####
# filter here for les than one year so aonly infections in first year of life are counted?
kids_with_malaria <- mic_drop %>%
  filter(!is.na(mstatus) & mstatus !=0, ageyrs<=1.02)%>%
  # mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_malaria_episodes")%>%
  filter(!duplicated("id"))%>%
  select(id, total_malaria_episodes)

# DPSP is a trial with a 1:1:1 split, DP, SP, DPSP
# MICDROP is a trial with a 1:1:1 split DP1Y, DP2Y, Placebo
kids_still_in_study <- mic_drop %>%
  group_by(id)%>%
  mutate("enrolled"=if_else(any(!is.na(pdeath) | !is.na(withdrawaldate)), "not_enrolled", "still_enrolled"))%>%
  filter(enrolled=="still_enrolled")%>%
  select("id")%>%
  filter(!duplicated(id))
  

fake_randomisation <- data.frame("id"=unique(kids_still_in_study$id))
set.seed(1234)
fake_randomisation <- fake_randomisation %>%
  mutate("fake_maternal_treatment"=sample(x = c("DP", "SP", "DPSP"),
                                          size=nrow(fake_randomisation),
                                          prob = c(1/3, 1/3, 1/3),
                                          replace=TRUE),
         "fake_infant_treatment"=sample(x = c("DP1Y", "DP2Y", "Placebo"),
                                        size=nrow(fake_randomisation),
                                        prob = c(1/3, 1/3, 1/3),
                                        replace=TRUE),
         "fake_infant_dp"=ifelse(fake_infant_treatment %in% c("DP1Y", "DP2Y"), "DP", "Placebo"),
         "total_n_infection"=kids_with_malaria$total_malaria_episodes[match(.$id, kids_with_malaria$id)],
         total_n_infection=ifelse(is.na(total_n_infection), 0, total_n_infection))



old_simulation <- data.frame(matrix(nrow = 0, ncol=5))
# colnames(old_simulation) <- colnames(fake_randomisation)


for(i in 1:1000){
  for(j in seq(30, 80, by=10)){
  
    
    fake_randomisation %>%
      mutate("fake_maternal_treatment"=sample(x = c("DP", "SP", "DPSP"),
                                              size=nrow(fake_randomisation),
                                              prob = c(1/3, 1/3, 1/3),
                                              replace=TRUE),
             "fake_infant_treatment"=sample(x = c("DP1Y", "DP2Y", "Placebo"),
                                            size=nrow(fake_randomisation),
                                            prob = c(1/3, 1/3, 1/3),
                                            replace=TRUE),
             "fake_infant_dp"=ifelse(fake_infant_treatment %in% c("DP1Y", "DP2Y"), "DP", "Placebo"),
             "total_n_infection"=kids_with_malaria$total_malaria_episodes[match(.$id, kids_with_malaria$id)],
             total_n_infection=ifelse(is.na(total_n_infection), 0, total_n_infection))%>%
      filter(fake_maternal_treatment != "DP")%>%
      group_by(fake_maternal_treatment, fake_infant_dp)%>%
      slice_sample(n=j)%>%
      mutate("n"=paste("n = ", j, sep=""),
             "n_sim"=i) -> new_simulation
  
    old_simulation <- rbind(old_simulation, new_simulation)
}}

# old_simulation %>%
#   group_by(n, n_sim)%>%
  
simulation_summary <- old_simulation%>%
  group_by(n, n_sim)%>%
  count(name="infections_in_first_year", total_n_infection)%>%
  group_by(n, total_n_infection)%>%
  summarise("mean"=mean(infections_in_first_year), "stdv"=sd(infections_in_first_year), "se"=stdv/sqrt(nrow(.)))%>%
  mutate("upper_ci"=mean+(1.96*stdv),
         "lower_ci"=mean-(1.96*stdv))



simulation_plot <- ggplot(simulation_summary, aes(x=total_n_infection, y=mean, color=n))+
  geom_point()+
  # geom_line()+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.5, linewidth = 0) +
  scale_x_continuous(breaks = seq(0, max(kids_with_malaria$total_malaria_episodes)))+
  scale_y_continuous(limits = c(0,300))+
  scale_color_viridis_d(option = "D")+
  geom_text(aes(label=round(mean), vjust = -2))+
  facet_grid(n~1)+
  theme_minimal()+
  theme(legend.position = "none")


ggsave("~/postdoc/stanford/clinical_data/MICDROP/sampling_strategy/n_infection_simulation.png", width=6, height = 15, dpi = 444, bg="white")








# mic_drop%>%
#   filter(id %in% kids_with_malaria$id, !is.na(mstatus))%>%
#   mutate(mstatus = case_match(mstatus,
#                               0~"asymptomatic",
#                               1~"uncomplicated",
#                               2~"complicated",
#                               3~"uncomplicated",
#                               4~"uncomplicated"))%>%
#   ggplot(., aes(x=AGE, y=pardens+0.001))+
#   geom_point(aes(color=factor(mstatus)))+
#   geom_line()+
#   scale_x_continuous(limits = c(0, 54))+
#   scale_y_log10()+
#   facet_wrap(~id)+
#   scale_color_manual(values=comp_pal)
#   theme_minimal()
