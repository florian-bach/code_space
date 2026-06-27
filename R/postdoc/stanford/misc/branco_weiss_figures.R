
(uncomplicated_comp_plot_for_bw <- all_malaria_summary%>%
  ggplot(., aes(x=n_any_malaria, y=risk))+
  geom_point(color="darkred")+
  geom_ribbon(data=malaria_summary_model_prediction$se, aes(x=n_any_malaria, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = malaria_summary_model_prediction$fun, colour="black")+
  # geom_text(aes(y=0.17, label= paste0("frac(",comp_1, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits = c(1,18))+
  scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  # ggtitle("Placebo")+
  # ggtitle(paste(sum(all_malaria_summary$comp_0+all_malaria_summary$comp_1), " malaria episodes (", sum(all_malaria_summary$comp_1), " complicated)", sep = ""))+
  xlab("Order of Malaria Episode")+
  ylab("Risk of Complicated Malaria During Episode"))+
  theme_minimal(base_size = 13)+
  theme(axis.text.x = element_text(size=11))

ggsave("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Florian_Grant_Applications/2026/Branco_Weiss/figures/uncomplicated_comp_plot_for_bw.png", width=4.5, height=4, dpi=444, bg="white")






mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/Specimens/Jun25/MICDSpecimenBoxJun25_withclinical.dta")

mic_drop_data <- mic_drop %>%
  filter(mstatus != 0, !is.na(mstatus), incidentmalaria==1 | is.na(incidentmalaria) & mstatus==2 | incidentmalaria==0 & mstatus==2)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(date=as.character(date))%>%
  select(id, date, n_infection)

infections <- clean_data%>%
  filter(mstatus==1)%>%
  distinct(id, date, ageinwks)


n_infections <- left_join(infections, mic_drop_data, by=c("id", "date"))


analytes <- unique(clean_data$targetName)



nulisa_to_show <- clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  # mutate(n_infection=ifelse(n_infection>=4, "4+", n_infection))%>%
  filter(mstatus==1, !is.na(n_infection), targetName %in% c("CTLA4", "CD274", "CD80", "IL22"))%>%
  mutate(targetNamef=factor(targetName, levels=c("CTLA4", "CD274", "CD80", "IL22")))%>%
  ggplot(aes(x=factor(n_infection), y=conc, fill=factor(n_infection)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetNamef, scales = "free", nrow=2)+
  scale_fill_manual(values=viridis::rocket(n = 5))+
  xlab("Order of Malaria Episode")+
  ylab("NULISA Units in  Blood")+
  theme_minimal(base_size = 13)+
  theme(legend.position = "none")

ggsave("~/Library/CloudStorage/Box-Box/Jagannathan_Lab_Folder/PROJECTS/Florian_Grant_Applications/2026/Branco_Weiss/figures/nulisa_during_malaria.png", nulisa_to_show, height=3.2, width=3.6, dpi = 444, bg="white")

for(i in 1:5){
  
  analytes_to_show <- 0:100[((i-1)*9)+1:i*9]

  print(analytes_to_show)
}

(nulisa9 <- clean_data %>%
  left_join(., n_infections, by=c("id", "date"))%>%
  mutate(n_infection=ifelse(n_infection>=4, "4+", n_infection))%>%
  filter(mstatus==1)%>%
  filter(targetName %in% analytes_to_show%>%
  ggplot(aes(x=factor(n_infection), y=conc, fill=factor(mstatus)))+
  geom_boxplot(outliers = FALSE)+
  facet_wrap(~targetName, scales = "free", nrow=2)+
  scale_fill_manual(values=viridis::rocket(n = 5))+
  xlab("Order of infection")+
  ylab("Concentration in Blood")+
  theme_minimal(base_size = 16)+
  theme(legend.position = "none"))
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/n_infection/", seq( ((i*9)-8), i*9)[1], "_", seq( ((i*9)-8), i*9)[9], ".png", sep=""), nulisa9, height=8, width=8, dpi = 444, bg="white")
}


