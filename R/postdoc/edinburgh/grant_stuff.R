library(pwr)
library(ggplot2)

# parasitaemia power plot ####
vac63_pmr <- data.frame("N_infection"=c("First", "Second", "Third"),
                        "mean"=c(9.87, 9.54, 8.26),
                        "SD"=c(1.927735036, 2.8272888332, 2.5396508946))


power_values <- c(seq(0.5, 0.9, by=0.1))
effect_size_values <- seq(0.6, 0.9, by=0.01)


results_df <- data.frame(matrix(ncol=4))

colnames(results_df) <- c("sample_size", "effect_size", "power", "percentage_shift")

for(power_value in power_values){
  
  for(effect_value in effect_size_values){
  
new_test <- pwr.t.test(n=NULL, 
           d = abs(9.87-(9.87*effect_value))/vac63_pmr$SD[1],
           sig.level = 0.05, 
           power = power_value,
           type=c("paired"),
           alternative = "two.sided")

new_row <- c(new_test$n, new_test$d, new_test$power, 1-effect_value)
results_df <- rbind(results_df, new_row)

}
}


results_df <- results_df[-1,]
# 
# powers_to_show <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# 
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# 

#results_df <- results_df[results_df$sample_size<=30,]

(power_plot <- ggplot(results_df, aes(x=percentage_shift, y=sample_size, colour=factor(power), group=factor(power)))+
  geom_line()+
  theme_minimal()+
  scale_y_continuous(breaks=c(seq(0,30, by=5)))+
  scale_x_continuous(labels = scales::label_percent(accuracy = 1))+
#  scale_color_manual(values=hcl.colors(7, palette = "Sunset"))+
  ylab("Sample Size")+
  xlab("Percentage Shift in Parasitaemia")+
  guides(colour=guide_legend(title="Power", reverse = TRUE))+
  coord_cartesian(ylim=c(0,30))
  )
ggsave("~/postdoc/grant_stuff/wiebke_power_plot.pdf", power_plot, height=4, width=4)
  

# goncalves figure ####


goncalves_data <- data.frame("N_Infection"=seq(0,14),
                             "Severe"=c(0, 32, 18, 20, 9, 9, 2, 2, 8, NA, 1, NA, NA, NA, 1),
                             "Infections"=c(0, 715, 540, 424, 330, 267, 220, 168, 131, NA, 90, NA, NA, NA, 25),
                             "Susceptible"=c(102, 70, 52, 32, 23, 14, 12, 10, 2, 2, 1, 1, 1, 1, 0),
                             "Severe_Moderate"=c(0, 76, 41, 31, 15, 12, 5, 8, 5, 1, 1, 3, NA, NA, 1),
                             "Supp_Infections"=c(0, 715, 501, 369, 276, 213, 169, 120, 86, 62, 55, 47, NA, NA, 25))


ggplot()+
  geom_point(data=goncalves_data[2:15,], aes(x=N_Infection, y=Supp_Infections), colour="blue")+
  geom_line(data=prd, aes(x=N_Infection, y=modelled_cases), colour="purple", linetype="dashed")+
  #scale_y_log10()+
  theme_minimal()+
  scale_x_continuous(breaks = seq(1,14))
  
  


goncalves_data <- goncalves_data%>%
  replace(is.na(.), 0)

goncalves_plot_data <- subset(goncalves_data, goncalves_data$N_Infection>0&goncalves_data$N_Infection<12)

fit <- glm(Severe_Moderate/Supp_Infections~N_Infection, weights = Supp_Infections, family="binomial", data=goncalves_plot_data)
moderate_exp_fun <- function(x){exp(fit$coefficients[1])*exp(fit$coefficients[2])^x}

prd <- data.frame(N_Infection = seq(from = 1, to = 11, length.out = 100))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit



pi<-summary(exponential_decay_estimate)@coef[1]
pi_se <- summary(exponential_decay_estimate)@coef[3]
alpha <- summary(exponential_decay_estimate)@coef[2]
alpha_se <- summary(exponential_decay_estimate)@coef[4]

prd2 <-  data.frame(N_Infection = seq(from = 1, to = 11, length.out = 100))
prd2$fit <- exp_model(prd2$N_Infection)
prd2$uci <- (pi+1.96*pi_se)*exp((-alpha+1.96*alpha_se)*prd2$N_Infection)
prd2$lci <- (pi-1.96*pi_se)*exp((-alpha-1.96*alpha_se)*prd2$N_Infection)


glm_model <- ggplot(goncalves_plot_data, aes(x = N_Infection, y=Severe_Moderate/Supp_Infections))+
  ylab("\n\nRisk of Severe / Moderate Disease\n")+
  scale_y_continuous(label=scales::percent, limits=c(0, 0.15))+
  scale_x_continuous(limits=range(0, nrow(goncalves_plot_data)), breaks=seq(0, nrow(goncalves_plot_data)), expand = expansion(add = c(0.2,0)))+
  geom_point(colour="blue")+
  ggtitle("Severe / Moderate Disease")+
  theme_minimal()+
  geom_ribbon(data=prd, aes(x=N_Infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  # geom_ribbon(data=prd2, aes(x=N_Infection, ymin = lci, ymax = uci),
  #             alpha = 0.2, inherit.aes = FALSE)+
  # geom_function(fun = exp_model, color="blue")+
  geom_function(fun = moderate_exp_fun, colour="purple")+
  theme(plot.title = element_text(size=15, hjust = 0.5),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13.5)
        #plot.margin = unit(c(2,2,2,2), "lines")
  )
  


ggsave("~/postdoc/grant_stuff/goncalves_figure.pdf", glm_model, height = 4, width=4*1.8)



# fraction of cd4 t cell memory plot, go to ~/code_space/R/vac69b/vac69b_activation.R ####

vac69ab_memory_fraction_activated <- memory_acti_data %>%
       group_by(volunteer, timepoint, n_infection)%>%
       summarise("percentage_of_cd4_memory"=sum(perc_of_memory))

write.table(vac69ab_memory_fraction_activated, "~/postdoc/grant_stuff/vac69ab_memory_fraction_activated.csv", sep=",", row.names = FALSE)


# neut power plot ####

neut_mean <- 2.110531047
neut_sd <- 0.3523358831

power_values <- c(seq(0.5, 0.9, by=0.1))
sample_size_values <- seq(8, 12, by=1)


neut_results_df <- data.frame(matrix(ncol=3))

colnames(neut_results_df) <- c("sample_size", "effect_size", "power")

for(power_value in power_values){
  
  for(n_value in sample_size_values){
    
    new_test <- pwr.t.test(n=n_value, 
                           d = NULL,
                           sig.level = 0.05, 
                           power = power_value,
                           type=c("two.sample"),
                           alternative = "two.sided")
    
    new_row <- c(new_test$n, new_test$d, new_test$power)
    neut_results_df <- rbind(neut_results_df, new_row)
    
  }
}

neut_results_df$percentage_shift <- neut_results_df$effect_size*neut_sd/neut_mean*100
#abs(neut_mean-neut_mean*effect_value)/neut_sd

neut_results_df <- neut_results_df[-1,]

write.table(neut_results_df, "~/postdoc/grant_stuff/neutralising_power_calc.csv", sep=",", row.names=FALSE)

# powers_to_show <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# 
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# 

#results_df <- results_df[results_df$sample_size<=30,]

(neut_power_plot <- ggplot(neut_results_df, aes(x=percentage_shift/100, y=sample_size, colour=factor(power), group=factor(power)))+
    geom_line()+
    theme_minimal()+
    scale_x_continuous(labels = scales::label_percent(accuracy = 1))+
    #scale_color_brewer(palette = pal_unikn_dark)+
    scale_color_manual(values=hcl.colors(7, palette = "Sunset"))+
    ylab("Sample Size")+
    xlab("Percentage Shift in [Neutralising Antibodies]")+
    guides(colour=guide_legend(title="Power", reverse = TRUE))+
    coord_cartesian(xlim = c(0.1, 0.3))
)
ggsave("~/postdoc/grant_stuff/neut_wiebke_power_plot.pdf", neut_power_plot, height=4, width=4)




# total ab power plot ####

ab_mean <- 2.47237292
ab_sd <- 0.35724909

power_values <- c(seq(0.5, 0.9, by=0.1))
sample_size_values <- seq(8, 12, by=1)

ab_results_df <- data.frame(matrix(ncol=3))

colnames(ab_results_df) <- c("sample_size", "effect_size", "power")

for(power_value in power_values){
  
  for(n_value in sample_size_values){
    
    new_test <- pwr.t.test(n=n_value, 
                           d = NULL,
                           sig.level = 0.05, 
                           power = power_value,
                           type=c("two.sample"),
                           alternative = "two.sided")
    
    new_row <- c(new_test$n, new_test$d, new_test$power)
    ab_results_df <- rbind(ab_results_df, new_row)
    
  }
}

ab_results_df$percentage_shift <- ab_results_df$effect_size*ab_sd/ab_mean*100


ab_results_df <- ab_results_df[-1,]

write.table(ab_results_df, "~/postdoc/grant_stuff/total_ab_power_calc.csv", sep=",", row.names=FALSE)




(total_ab_power_plot <- ggplot(ab_results_df, aes(x=percentage_shift/100, y=sample_size, colour=factor(power), group=factor(power)))+
    geom_line()+
    theme_minimal()+
    scale_x_continuous(labels = scales::label_percent(accuracy = 1))+
    #scale_color_brewer(palette = pal_unikn_dark)+
    scale_color_manual(values=hcl.colors(7, palette = "Sunset"))+
    ylab("Sample Size")+
    xlab("Percentage Shift in [Total IgG]")+
    guides(colour=guide_legend(title="Power", reverse = TRUE))+
    coord_cartesian(xlim = c(0.1, 0.3))
)
ggsave("~/postdoc/grant_stuff/total_igg_wiebke_power_plot.pdf", total_ab_power_plot, height=4, width=4)


# vivax pmr power plot ####



vivax_pmr <- readxl::read_excel("~/postdoc/grant_stuff/vivax_PMR.xlsx")

vivax_mean <- mean(vivax_pmr$vivax_PMR)
vivax_sd <- sd(vivax_pmr$vivax_PMR)




power_values <- c(seq(0.5, 0.9, by=0.1))
effect_size_values <- seq(0.6, 0.9, by=0.01)


vivax_results_df <- data.frame(matrix(ncol=4))

colnames(vivax_results_df) <- c("sample_size", "effect_size", "power", "percentage_shift")

for(power_value in power_values){
  
  for(effect_value in effect_size_values){
    
    new_test <- pwr.t.test(n=NULL, 
                           d = abs(vivax_mean-(vivax_mean*effect_value))/vivax_sd,
                           sig.level = 0.05, 
                           power = power_value,
                           type=c("paired"),
                           alternative = "two.sided")
    
    new_row <- c(new_test$n, new_test$d, new_test$power, 1-effect_value)
    vivax_results_df <- rbind(vivax_results_df, new_row)
    
  }
}


vivax_results_df <- vivax_results_df[-1,]
# 
# powers_to_show <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
# 
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# results_df_one <- dplyr::filter(results_df, results_df$power %in% powers_to_show)
# 

#results_df <- results_df[results_df$sample_size<=30,]

(vivax_power_plot <- ggplot(vivax_results_df, aes(x=percentage_shift, y=sample_size, colour=factor(power), group=factor(power)))+
    geom_line()+
    theme_minimal()+
    scale_y_continuous(breaks=c(seq(0,30, by=5)))+
    scale_x_continuous(labels = scales::label_percent(accuracy = 1))+
    #scale_color_brewer(palette = pal_unikn_dark)+
    scale_color_manual(values=hcl.colors(7, palette = "Sunset"))+
    ylab("Sample Size")+
    xlab("Percentage Shift in Parasitaemia")+
    guides(colour=guide_legend(title="Power", reverse = TRUE))+
    coord_cartesian(ylim=c(0,30))
)

ggsave("~/postdoc/grant_stuff/wiebke_vivax_power_plot.pdf", vivax_power_plot, height=4, width=4)
