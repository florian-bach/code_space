model_visualiser <- function(model=NULL, xvar=NULL){
  model_data = model$data
  xvar_range = range(model_data[[xvar]])
  pre_df = data.frame(seq(from = xvar_range[1], to = xvar_range[2], length.out = 100))
  colnames(pre_df)=xvar
  error = predict(model, newdata = pre_df, se.fit = TRUE)
  pre_df$lci = error$fit - 1.96 * error$se.fit
  pre_df$fit = error$fit
  pre_df$uci = error$fit + 1.96 * error$se.fit
  return(pre_df)
}



goncalves_data <- data.frame("N_Infection"=seq(0,14),
                             "Severe"=c(0, 32, 18, 20, 9, 9, 2, 2, 8, NA, 1, NA, NA, NA, 1),
                             "Infections"=c(0, 715, 540, 424, 330, 267, 220, 168, 131, NA, 90, NA, NA, NA, 25),
                             "Susceptible"=c(102, 70, 52, 32, 23, 14, 12, 10, 2, 2, 1, 1, 1, 1, 0),
                             "Severe_Moderate"=c(0, 76, 41, 31, 15, 12, 5, 8, 5, 1, 1, 3, NA, NA, 1),
                             "Supp_Infections"=c(0, 715, 501, 369, 276, 213, 169, 120, 86, 62, 55, 47, NA, NA, 25))

goncalves_data$Severe_Moderate <- ifelse(is.na(goncalves_data$Severe_Moderate), 0, goncalves_data$Severe_Moderate)

goncalves_data$Severe <- ifelse(is.na(goncalves_data$Severe), 0, goncalves_data$Severe)

# goncalves_plot_data <- subset(goncalves_data, goncalves_data$N_Infection>0&goncalves_data$N_Infection<12)

supp_infection_decay <- lm(log(Supp_Infections)~N_Infection, data=goncalves_data[2:15,], weights = Supp_Infections)

supp_prd <- data.frame(N_Infection = seq(from = 1, to = 14), length.out=1000)
supp_err <- predict(supp_infection_decay, newdata = supp_prd, se.fit = TRUE)
supp_prd$modelled_cases <- exp(supp_err$fit)
supp_prd$round_model <- round(supp_prd$modelled_cases)

goncalves_data$Modelled_Supp_Cases <- ifelse(is.na(goncalves_data$Supp_Infections), supp_prd$round_model, goncalves_data$Supp_Infections)

goncalves_data$Supp_Imputed <- ifelse(is.na(goncalves_data$Supp_Infections), TRUE, FALSE)





main_infection_decay <- lm(log(Infections)~N_Infection, data=goncalves_data[2:15,], weights = Infections)

main_prd <- data.frame(N_Infection = seq(from = 1, to = 14), length.out=100)
main_err <- predict(main_infection_decay, newdata = main_prd, se.fit = TRUE)
main_prd$modelled_cases <- exp(main_err$fit)
main_prd$round_model <- round(main_prd$modelled_cases)

goncalves_data$Modelled_Main_Cases <- ifelse(is.na(goncalves_data$Infections), main_prd$round_model, goncalves_data$Infections)

goncalves_data$Main_Imputed <- ifelse(is.na(goncalves_data$Infections), TRUE, FALSE)






moderate_glm <- glm(Severe_Moderate/Modelled_Supp_Cases~N_Infection, weights = Modelled_Supp_Cases, family="binomial", data=goncalves_data)
moderate_exp_fun <- function(x){exp(moderate_glm$coefficients[1])*exp(moderate_glm$coefficients[2])^x}




comp_glm_model <- ggplot(goncalves_data, aes(x = N_Infection, y=Severe_Moderate/Supp_Infections))+
  geom_point(colour="darkred")+
  geom_ribbon(data= model_visualiser(moderate_glm, "N_Infection"), aes(x=N_Infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = moderate_exp_fun, colour="black")+
  geom_text(aes(y=0.14, label= paste0("frac(",Severe_Moderate, ",", Modelled_Supp_Cases,")")),parse = TRUE, size=2.5)+
  ylab("Risk of Complicated Malaria")+
  xlab("Order of Infection")+
  scale_y_continuous(label=scales::percent, limits=c(0, 0.15))+
  scale_x_continuous(limits=range(1, nrow(goncalves_data)), breaks=seq(0, nrow(goncalves_data)))+
  theme_minimal()+
  theme(plot.title = element_text(size=15, hjust = 0.5),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13.5)
  )

ggsave("~/postdoc/edinburgh/goncalves_modelling/complicated_malaria.pdf", comp_glm_model, height = 4, width=5, dpi=444, bg="white")




severe_glm <- glm(Severe/Modelled_Main_Cases~N_Infection, weights = Modelled_Main_Cases, family="binomial", data=goncalves_data)
severe_exp_fun <- function(x){exp(severe_glm$coefficients[1])*exp(severe_glm$coefficients[2])^x}


(
severe_glm_model <- ggplot(goncalves_data, aes(x = N_Infection, y=Severe/Modelled_Main_Cases))+
  geom_point(colour="darkred")+
  geom_ribbon(data= model_visualiser(severe_glm, "N_Infection"), aes(x=N_Infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = severe_exp_fun, colour="black")+
  geom_text(aes(y=0.07, label= paste0("frac(",Severe, ",", Modelled_Main_Cases,")")),parse = TRUE, size=2.5)+
  ylab("Risk of Complicated Malaria")+
  xlab("Order of Infection")+
  scale_y_continuous(label=scales::percent, limits=c(0, 0.075))+
  scale_x_continuous(limits=range(1, nrow(goncalves_plot_data)), breaks=seq(0, nrow(goncalves_plot_data)))+
  theme_minimal()+
  theme(plot.title = element_text(size=15, hjust = 0.5),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13.5)
  )
)
ggsave("~/postdoc/edinburgh/goncalves_modelling/severe_malaria.png", severe_glm_model, height = 4, width=5, dpi=444, bg="white")






supp_case_decay <- ggplot()+
  geom_point(data=goncalves_data, colour="darkblue", aes(x=N_Infection, y=Modelled_Supp_Cases, shape=Supp_Imputed))+
  geom_line(data = supp_prd, aes(x = N_Infection, y=round_model), colour="blue", linetype="dashed")+
  ylab("\n\nNumber of Cases\n")+
  scale_x_continuous(breaks=seq(1, nrow(goncalves_plot_data)-1),
                     limits=range(1,nrow(goncalves_plot_data)-1),
                     expand = expansion(add = c(0.2,0.2)))+
  ggtitle("Number of Cases with Imputed Values")+
  theme_minimal()+
  guides(shape=guide_legend(title="Imputed"))+
  theme(plot.title = element_text(size=15, hjust = 0.5),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13.5)
        #plot.margin = unit(c(2,2,2,2), "lines")
  )
