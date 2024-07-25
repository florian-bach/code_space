promote <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL.dta")
bc3 <- haven::read_dta("~/postdoc/stanford/clinical_data/PROMOTE/BC-3 childs all visit database FINAL_withIncidentNMF.dta")

promote_data <- promote %>%
  filter(!is.na(parsdens), !is.na(mstatus), parsdens!=0)%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  mutate(n_infection = seq(1, max(total_n_infection)),
         mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate("hbs"=bc3$hbs[match(id, bc3$id)],
         hbs=case_match(hbs,
                        1~"HbAA",
                        2~"HbAS",
                        3~"HbSS"))%>%
  mutate("age"=date-dob)%>%
  mutate(pardens=parsdens)%>%
  dplyr::select(-parsdens)%>%
  filter(!is.na(hbs))

         

all_malaria <- promote_data %>%
  # filter(mstatus != "no malaria", !is.na(mstatus))%>%
  group_by(id)%>%
  mutate("age_at_first"=age[n_infection==1])%>%
  ungroup()%>%
  mutate("treatment_failure"=if_else(mstatus%in% c("quinine for AL failure", "Q/AS failure"), 1, 0))%>%
  mutate(agebins=cut(as.numeric(age), breaks = seq(0, max(as.numeric(age)), by=30)))%>%
  mutate(age_at_first_bins=cut(as.numeric(age_at_first), breaks = seq(0, 360, by=90), labels = three_month_labels[1:4]))%>%
  mutate(disease=if_else(mstatus=="complicated", "complicated", if_else(mstatus=="no malaria", "asymptomatic", "uncomplicated")),
         complicated=if_else(mstatus=="complicated", 1, 0),
         clinical=if_else(mstatus!="no malaria", 1, 0))%>%
  mutate(id=factor(id),
         age=as.numeric(age),
         log_pardens=log10(pardens+0.1))



comp_model <- glm(complicated~n_infection+I(n_infection^2)+hbs, family = "binomial", data = all_malaria)

aa_only <- filter(all_malaria, hbs=="HbAA")
comp_model_aa <- glm(complicated~n_infection+I(n_infection^2), family = "binomial", data = aa_only)


combo_comp <- all_malaria %>%
  group_by(disease, n_infection, hbs)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections)


hbs_comp <- combo_comp%>%
  filter(hbs!="HbSS")%>%
  ggplot(aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  # geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = comp_model_fun, colour="black")+
  geom_text(aes(y=0.15, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,15))+
  scale_y_continuous(limits = c(0,0.17), labels = scales::label_percent())+
  geom_smooth(method="glm", formula = y~x+I(x^2))+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")+
  facet_wrap(~hbs)+
  ggtitle("374 HbAA, 72 HbAS children, aged 8 weeks - 2 years\n1444 and 246 parasitemic episodes, respectively")

ggsave("~/postdoc/stanford/clinical_data/hbs_comp.png", hbs_comp, width=8, height=4, bg="white", dpi=444)





#save model as function for plotting
comp_model_aa_fun <- function(x){
  exp(comp_model_aa$coefficients[1])*
    exp(comp_model_aa$coefficients[2])^x*
    exp(comp_model_aa$coefficients[3])^x^2}

#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 13, length.out = 100))
err <- predict(comp_model_aa, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit





aa_only_comp <- aa_only %>%
  group_by(disease, n_infection)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated),
         total_infections=complicated+uncomplicated+asymptomatic,
         risk=complicated/total_infections)

ggplot(aa_only_comp, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_aa_fun, colour="black")+
  geom_text(aes(y=0.15, label= paste0("frac(",complicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  scale_x_continuous(breaks = 1:50, limits=c(1,13))+
  scale_y_continuous(limits = c(0,0.17), labels = scales::label_percent())+
  xlab("Order of Infection")+
  ylab("Risk of Complicated")+
  # facet_wrap(~hbs)+
  ggtitle("374 HbAA children, aged 8 weeks - 2 years\n 3430 parasitemic episodes")


# clinical immunity ####


combo_clin <- all_malaria %>%
  group_by(disease, n_infection, hbs)%>%
  summarise("n"=n())%>%
  pivot_wider(names_from = disease, values_from = n)%>%
  mutate(complicated=if_else(is.na(complicated), 0, complicated))%>%
  mutate(uncomplicated=if_else(is.na(uncomplicated), 0, uncomplicated))%>%
  mutate(asymptomatic=if_else(is.na(asymptomatic), 0, asymptomatic))%>%
  ungroup()%>%
  group_by(hbs, n_infection)%>%
  mutate(#complicated=if_else(is.na(complicated), 0, complicated),
    total_infections=sum(complicated,uncomplicated,asymptomatic,na.rm = TRUE),
    risk=sum(complicated,uncomplicated, na.rm = TRUE)/total_infections
  )


clin_model <- glm(risk~n_infection+hbs, family = "binomial", weights = total_infections, data = combo_clin[1:10,])

no_hbss <-combo_clin%>%
  filter(hbs != "HbSS")

clin_model <- glm(clinical~n_infection+hbs, family = "binomial", data = all_malaria)

clin_model_hbaa <- function(x){
  exp(clin_model$coefficients[1])*
    # exp(clin_model$coefficients[3])*
    exp(clin_model$coefficients[2])^x
    }



clinical_malaria_risk_hbs <-  ggplot(no_hbss, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_smooth(method="glm", method.args=list(family="binomial", weights=no_hbss$total_infections))+
  facet_wrap(~hbs)+
  geom_function(fun = clin_model_hbaa, colour="black")+
  scale_x_continuous(breaks = 1:50, limits=c(1,10))+
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, NA))+
  geom_text(aes(y=1.1, label= paste0("frac(",complicated+uncomplicated, ",", total_infections,")")),parse = TRUE, size=2.5)+
  # ggtitle("903 children, aged 8 weeks - 2 years\n3438 parsitemic episodes")+
  xlab("Order of Infection")+
  ylab("Risk of Malaria When Parasitemic")
  
ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/clinical_malaria_risk_hbs.png", clinical_malaria_risk_hbs, height = 4.5, width=7.5, dpi=444, bg="white")




hbs_n_infection_para <- all_malaria %>%
  filter(n_infection<=10, hbs!= "HbSS")%>%
  # mutate(quarter=cut(as.numeric(factor(agebins)),breaks = seq(0, 24, by=3)))%>%
  ggplot(aes(x=n_infection, y=pardens+0.1))+
  # geom_point(aes(color=age_at_first_bins))+
  geom_boxplot(aes(fill=factor(n_infection)), outlier.shape = NA)+
  scale_y_log10()+
  scale_x_continuous(breaks=seq(1, 10))+
  # geom_smooth(method="lm", formula = y~x+I(x^2))+
  facet_wrap(~hbs)+
  scale_fill_manual(values=colorspace::sequential_hcl(n=25, palette = "Lajolla"))+
  xlab("order of infection")+
  ylab("parasites / mL")+
  theme(#axis.text.x = element_text(angle=90),
    legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/complicated_malaria/hbs_n_infection_para.png", hbs_n_infection_para, height = 4, width=8, dpi=444, bg="white")


