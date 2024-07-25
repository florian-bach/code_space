impact <- haven::read_dta("~/postdoc/stanford/clinical_data/IMPACT/IMPACT all visits database through January 31st 2024.dta")


complicated_promote <- impact %>%
  filter(mstatus %in% c(1,2))%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_n_infection") %>%
  arrange(AGE) %>%
  mutate(n_infection = seq(1, max(total_n_infection)))%>%
  mutate(disease=case_when(mstatus==2 ~ "complicated",
                           mstatus==1 ~ "uncomplicated"))%>%
  group_by(id) %>%
  #there are 2 complicated episodes that are non-incident, and 11 total. we keep the later one for all, allowing resolution of complicated episodes.
  # mutate(when_not_incident=ifelse(any(incidentmalaria=="not_incident"), which(incidentmalaria=="not_incident"), 200))%>%
  # mutate(n_infection=ifelse(n_infection>=when_not_incident, n_infection-1, n_infection))%>%
  group_by(id, n_infection)%>%
  slice_max(order_by = date, n=1, with_ties = FALSE)


promote_total_infections <- complicated_promote %>%
  group_by(n_infection) %>%
  count(name = "total_infections")

promote_complicated_df <- as.data.frame(table(complicated_promote$n_infection, complicated_promote$disease))
colnames(promote_complicated_df) <- c("n_infection", "disease", "complicated_episodes")
promote_complicated_df <- subset(promote_complicated_df, promote_complicated_df$disease=="complicated")
promote_complicated_df$total_infections <- promote_total_infections$total_infections
promote_complicated_df$risk <- promote_complicated_df$complicated_episodes/promote_complicated_df$total_infections
promote_complicated_df$n_infection <- as.numeric(promote_complicated_df$n_infection)

# binomial regression
comp_model <- glm(risk~n_infection+I(n_infection^2), family = "binomial", weights = total_infections, data = promote_complicated_df)
# weighted linear regression with quadratic term

#save model as function for plotting
comp_model_fun2 <- function(x){
  exp(comp_model2$coefficients[1])*
    exp(comp_model2$coefficients[2])^x*
    exp(comp_model2$coefficients[3])^x^2}



#calculate SE for plotting
prd <- data.frame(n_infection = seq(from = 1, to = 10, length.out = 100))
err <- predict(comp_model, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit



ggplot(complicated_df, aes(x=n_infection, y=risk))+
  geom_point(color="darkred")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  ylab("Risk of Complicated Episodes")+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(limits = c(0,0.6))+
  xlab("Order of Infection")
