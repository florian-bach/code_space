
mic_drop <-  haven::read_dta("~/Library/CloudStorage/Box-Box/MIC_DroP IPTc Study/Data/MICDroP Data/MICDROP expanded database through November 30th 2024.dta")

mic_drop_hbs <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP SickleTr final.dta")

mic_drop_data <- mic_drop %>%
  filter(!is.na(qPCRparsdens), !is.na(mstatus))%>%
  group_by(id) %>%
  add_count(name="total_n_visits") %>%
  mutate(n_visit = seq(1, max(total_n_visits)))%>%
  mutate("total_n_para"=sum(qPCRparsdens!=0),
         "total_n_malaria"=sum(qPCRparsdens!=0),
         "n_para"=if_else(pardens!=0, cumsum(qPCRparsdens!=0), NA),
         "n_malaria"=if_else(mstatus!=0, cumsum(qPCRparsdens!=0), NA))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  ungroup()%>%
  mutate(hbs=mic_drop_hbs$HbS[match(as.numeric(id), mic_drop_hbs$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))


mic_drop_for_merge <- mic_drop_data %>%
  dplyr::select(id, date, dob, hbs, mstatus, qPCRparsdens, total_n_visits, total_n_para, total_n_malaria, n_para, n_malaria, temp)%>%
  mutate("study"="micdrop")%>%
  mutate("age"=date-dob)

combo_data <- mic_drop_for_merge


three_month_labels <- paste0(seq(0, 21, by=3), " to ", seq(3, 24, by=3), " months")

all_malaria <- combo_data %>%
  ungroup()%>%
  filter(total_n_malaria>=1)%>%
  group_by(id)%>%
  mutate("age_at_second"=ifelse(any(n_malaria==2), age[n_malaria==2], NULL))%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(age_at_second_bins=cut(as.numeric(age_at_second), breaks = seq(0, 730, by=90), labels = three_month_labels))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", if_else(mstatus=="no malaria", "asymptomatic", "uncomplicated")),
         complicated=if_else(mstatus=="complicated", 1, 0), 
         symptoms=if_else(disease != "asymptomatic", 1, 0))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         age_months=as.numeric(factor(agebins)),
         log_pardens=log10(qPCRparsdens+0.1))


n_para_comp <- all_malaria %>%
  group_by(disease, n_para)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections,
         asymp_prob=asymptomatic/total_infections)%>%
  filter(n_para>0 & n_para <12)


comp_model <- glm(risk~n_para+I(n_para^2), family = "binomial", weights = total_infections, data = n_para_comp)

comp_model2 <- glm(complicated~n_para+age, family = "binomial", data = all_malaria)

comp_model3 <- lme4::glmer(complicated~n_para+I(n_para^2)+age_at_second+(1|id), family = "binomial", data = all_malaria)

## risk~n_infection ####

#save model as function for plotting
comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_para = seq(from = 1, to = 19, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit


n_para_comp_plot <- ggplot(n_para_comp[1:19,], aes(x=n_para, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_para, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.10, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,19))+
  scale_y_continuous(limits = c(0,0.12), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")+
  ggtitle("827 children, aged 8 weeks - 120 weeks
3541 parasitemic episodes")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/micdrop_only_qpcr_n_para_comp_plot.png", n_para_comp_plot, height = 4.5, width=7.5, dpi=444, bg="white")
