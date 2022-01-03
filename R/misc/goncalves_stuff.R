library(ggplot2)

# figures ####

goncalves_data <- data.frame("Infection"=seq(0,14),
                             "Severe"=c(0, 32, 18, 20, 9, 9, 2, 2, 8, NA, 1, NA, NA, NA, 1),
                             "Infections"=c(0, 715, 540, 424, 330, 267, 220, 168, 131, NA, 90, NA, NA, NA, 25),
                             "Susceptible"=c(102, 70, 52, 32, 23, 14, 12, 10, 2, 2, 1, 1, 1, 1, 0))


decay <- function(x){0.7^x}
decay3a <- function(x){0.78^(x)}
# decay4a <- function(x){0.6^(x)}
# decay5a <- function(x){0.5^(x)}

(kaplan_meyer_esque <- ggplot(goncalves_data, aes(x=Infection, y=Susceptible/102))+
  ylab("Percentage Susceptible to Severe Disease")+
  xlab("Order of Infection")+
  scale_y_continuous(label=scales::percent, expand = expansion(add = c(0.02,0)))+
  #geom_segment(aes(x=0, xend=2, y=0.5, yend=0.5), linetype="dashed", size=0.3)+
  #geom_segment(aes(x=2, xend=2, y=0, yend=0.5), linetype="dashed", size=0.3)+
  geom_function(fun=decay, inherit.aes = TRUE, colour="purple")+
  geom_function(fun=decay3a, inherit.aes = TRUE, colour="red")+
  #geom_function(fun=decay4a, inherit.aes = TRUE, colour="green")+
  #geom_function(fun=decay5a, inherit.aes = TRUE, colour="pink")+
  geom_text(aes(x=10, y=0.6), label=expression(paste("f(x)"~"="~0.78^x)), colour="red", size=6)+
  geom_text(aes(x=10, y=0.5), label=expression(paste("f(x)"~"="~0.7^x)), colour="purple", size=6)+
  scale_x_continuous(breaks = seq(0,14, by=1),
    #labels = c(seq(0,10,by=1), 12, 14),
    expand = expansion(add = c(0.2,0)))+
  #geom_smooth(color="black", method = "loess", fullrange=TRUE)+
  geom_point(color="black")+
  theme_minimal()+
    theme()
)

ggsave("~/postdoc/goncalves_modelling/png/kaplan_meyer_esque.png", kaplan_meyer_esque, height=4, width=3.5, bg="white", dpi=444)
ggsave("~/postdoc/goncalves_modelling/pdf/kaplan_meyer_esque.pdf",  kaplan_meyer_esque, height=4, width=3.5, bg="white")



decay2 <- function(x){715*0.7^(x-1)}
decay3 <- function(x){715*0.78^(x-1)}
decay4 <- function(x){715*0.6^(x-1)}
decay5 <- function(x){715*0.5^(x-1)}

(infection_decay <- ggplot(goncalves_data, aes(x = Infection, y=Infections))+
  ylab("Number of Infections")+
  geom_function(fun=decay2, inherit.aes = TRUE, colour="purple")+
    geom_function(fun=decay3, inherit.aes = TRUE, colour="red")+
   #geom_function(fun=decay4, inherit.aes = TRUE, colour="green")+
   #geom_function(fun=decay5, inherit.aes = TRUE, colour="pink")+
   #geom_text(aes(x=10, y=400), label=expression(paste("f(x)"~"="~715~"*"~0.78^(x-1))), colour="red", size=4)+
  geom_text(aes(x=10, y=350), label=expression(paste("f(x)"~"="~715~"*"~0.7^(x-1))), colour="purple", size=4)+
  geom_text(aes(x=10, y=400), label=expression(paste("f(x)"~"="~715~"*"~0.78^(x-1))), colour="red", size=4)+
    
  #scale_y_continuous(label=scales::percent, limits=c(0, 0.06))+
  scale_x_continuous(limits=c(1,14), breaks = seq(1,14, by=1), expand = expansion(add = c(0.2,0)))+
  geom_point()+
  theme_minimal())
  

ggsave("~/postdoc/goncalves_modelling/pdf/infection_decay.pdf", infection_decay, height=4, width=3.5, bg="white")
ggsave("~/postdoc/goncalves_modelling/png/infection_decay.png", infection_decay, height=4, width=3.5, bg="white", dpi=444)




ggplot(goncalves_data, aes(x = Infection, y=Severe/Infections))+
  ylab("Risk of Severe Disease")+
  scale_y_continuous(label=scales::percent, limits=c(0, 0.06))+
  scale_x_continuous(breaks = seq(0,14, by=1), expand = expansion(add = c(0.2,0)))+
  geom_point()+
  #geom_smooth()+
  theme_minimal()

# glm modelling ####

#goncalves_model_data <- subset(goncalves_data, Infection<9&Infection>0)
goncalves_model_data <- subset(goncalves_data, Infection<6&Infection>0)
goncalves_model_data$Severe_Div_Infection <- goncalves_model_data$Severe/goncalves_model_data$Infections
goncalves_model_data$Severe_Div_Susceptible <- goncalves_model_data$Severe/(goncalves_model_data$Susceptible+goncalves_model_data$Severe)

# negative binomial; P(severe | infection); first five infections, order of infection not significant; --> same risk
glm1 <- glm(Severe_Div_Infection~Infection, weights = Infections, data=goncalves_model_data, family = "binomial")

# negative binomial; first five infections, order of infection not significant; --> same risk
glm2 <- glm(Severe_Div_Susceptible~Infection, weights = Susceptible+Severe, data=goncalves_model_data, family = "binomial")


simulation_hundred <- data.frame("First"=NA, "Second"=NA, "Third"=NA, "Fourth"=NA, "Fifth"=NA)

for(i in 1:5000){
  new_entry <- rbinom(n=goncalves_model_data$Severe, size=goncalves_model_data$Infections, prob = median(goncalves_model_data$Severe_Div_Infection, na.rm = TRUE))
  #names(new_entry) <- c("First", "Second", "Third", "Fourth", "Fifth")
  simulation_hundred[i,]<-new_entry
  
}
  
long_simulation_hundred <- tidyr::pivot_longer(simulation_hundred, cols=c(First, Second, Third, Fourth, Fifth), names_to = "n_infection", values_to = "severe_cases")

bad_numbers <- c("First", "Second", "Third", "Fourth", "Fifth")
good_numbers <- as.character(seq(1:5))

long_simulation_hundred$n_infection_num <- stringr::str_replace(long_simulation_hundred$n_infection, bad_numbers, good_numbers)
# 

summary_df <- long_simulation_hundred %>%
  group_by(n_infection)%>%
  summarise("mean"=mean(severe_cases), "sd"=sd(severe_cases))

summary_df$mean_plus2 <- summary_df$mean+(summary_df$sd)*2
summary_df$mean_minus2 <- summary_df$mean-(summary_df$sd)*2

summary_df$mean_plus <- summary_df$mean+summary_df$sd
summary_df$mean_minus <- summary_df$mean-summary_df$sd

summary_df$n_infection_num <- c(5, 1, 4, 2, 3)

(simulation_hundred_plot <- ggplot(long_simulation_hundred, aes(x=as.numeric(n_infection_num), y=severe_cases))+
  #geom_point(aes(x=as.numeric(n_infection_num), y=severe_cases), position = position_jitter(width=0.1, height = 0), inherit.aes = FALSE)+
  geom_violin(aes(x=n_infection_num, y=severe_cases), inherit.aes = FALSE)+
  geom_smooth(method = "loess", se=FALSE)+
  geom_point(data=goncalves_model_data, aes(x = Infection, y=Severe), color="red", inherit.aes = FALSE)+
  #geom_point(data=summary_df, aes(x = Infection, y=Severe), color="red", inherit.aes = FALSE)+
  xlab("Order of Infection")+
  geom_ribbon(aes(x=n_infection_num, ymin = mean_minus2, ymax = mean_plus2), alpha = 0.1, data = summary_df, inherit.aes = FALSE)+
  geom_ribbon(aes(x=n_infection_num, ymin = mean_minus, ymax = mean_plus), alpha = 0.2, data = summary_df, inherit.aes = FALSE)+
  ylab("Number of Severe Cases")+
  theme_minimal())

ggsave("~/postdoc/goncalves_modelling/png/simulation_violin.png", simulation_hundred_plot, height=4, width=5, bg="white", dpi=444)
ggsave("~/postdoc/goncalves_modelling/pdf/simulation_violin.pdf", simulation_hundred_plot, height=4, width=5, bg="white")




median_normality <- long_simulation_hundred %>%
   group_by(n_infection)%>%
   slice_sample(n = 100)%>%
   ggpubr::ggqqplot(.$severe_cases)
   
