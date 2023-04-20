library(ggplot2)


# sequential acquisition of immunity ####

severe <- data.frame("severe"=dnorm(x=1:1000, mean = 2, sd = 1))

clinical <- data.frame("clinical"=rnorm(1000, mean=6, sd=2))

asymptomatic <- data.frame("asymptomatic"=rnorm(1000, mean=8, sd=3))



severe <- function(x){100/(0.3*(x-0.7)^2+1)}

moderate <- function(x){100/(0.05*(x-1.4)^2+1)}
mild <- function(x){100/(0.02*(x-3)^2+1)}


langhorne_plot <- ggplot()+
    stat_function(fun=mild, geom="area", aes(fill="asymptomatic"))+
  stat_function(fun=moderate, geom="area", aes(fill="moderate"))+
  stat_function(fun=severe, geom="area", aes(fill="severe"), n = 1000)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0.5, 20), expand = expansion(0,0))+
  ylab("percent of maximum")+
  xlab("age")+
  scale_fill_manual(breaks=c("severe", "moderate", "asymptomatic"),
                    values=c("#800000", "#E69F00", "#6488EA"))+
  ggtitle("malaria disease phenotypes")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank())


ggsave("~/postdoc/presentations/ashworth_symposium/immunity_sequence_schematic.png", height=4, width=6, dpi=644, bg="white")

# goncalves kaplan meyer plot ####

goncalves_data <- data.frame("Infection"=seq(0,14),
                             "Severe"=c(0, 32, 18, 20, 9, 9, 2, 2, 8, NA, 1, NA, NA, NA, 1),
                             "Infections"=c(0, 715, 540, 424, 330, 267, 220, 168, 131, NA, 90, NA, NA, NA, 25),
                             "Susceptible"=c(102, 70, 52, 32, 23, 14, 12, 10, 2, 2, 1, 1, 1, 1, 0),
                             "Severe_Moderate"=c(0, 76, 41, 31, 15, 12, 5, 8, 5, 1, 1, 3, NA, NA, 1))

goncalves_data$SevMod_Susceptible <- sum(goncalves_data$Severe_Moderate, na.rm = TRUE)-cumsum(ifelse(is.na(goncalves_data$Severe_Moderate),0, goncalves_data$Severe_Moderate))


decay <- function(x){0.69^x}
decay3a <- function(x){0.65^(x)}
# decay4a <- function(x){0.6^(x)}
# decay5a <- function(x){0.5^(x)}

kaplan_meyer_esque <- ggplot(goncalves_data, aes(x=Infection, y=Susceptible/102))+
    ylab("Percentage Susceptible to Severe Disease")+
    xlab("Order of Infection")+
    ggtitle("Severe Disease")+
    geom_function(fun=decay, inherit.aes = TRUE, colour="red")+
    #geom_function(fun=decay3a, inherit.aes = TRUE, colour="blue")+
    scale_y_continuous(label=scales::percent, expand = expansion(add = c(0.03,0.03)))+
    scale_x_continuous(breaks = seq(0,14, by=1),
                     #labels = c(seq(0,10,by=1), 12, 14),
                     expand = expansion(add = c(0.2,0)))+
    #geom_smooth(color="black", method = "loess", fullrange=TRUE)+
    geom_point(color="red")+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))

ggsave("~/postdoc/goncalves_modelling/png/kaplan_meyer_severe.png", kaplan_meyer_esque, height=4, width=3.5, bg="white", dpi=444)



kaplan_meyer_esque2 <- ggplot(goncalves_data, aes(x=Infection, y=SevMod_Susceptible/199))+
  ylab("Percentage Susceptible to Moderate / Severe Disease")+
  xlab("Order of Infection")+
  ggtitle("Moderate / Severe Disease")+
  scale_y_continuous(label=scales::percent, expand = expansion(add = c(0.03,0.03)))+
  scale_x_continuous(breaks = seq(0,14, by=1),
                     #labels = c(seq(0,10,by=1), 12, 14),
                     expand = expansion(add = c(0.2,0)))+
  #geom_function(fun=decay, inherit.aes = TRUE, colour="red")+
  geom_function(fun=decay3a, inherit.aes = TRUE, colour="blue")+
  geom_text(aes(x=10, y=0.6), label=expression(paste("f(x)"~"="~0.65^x)), colour="blue", size=6)+
  geom_text(aes(x=10, y=0.5), label=expression(paste("f(x)"~"="~0.69^x)), colour="red", size=6)+
  #geom_smooth(color="black", method = "loess", fullrange=TRUE)+
  geom_point(color="blue")+
  theme_minimal()+
  theme()+
  theme(plot.title = element_text(hjust=0.5))

ggsave("~/postdoc/goncalves_modelling/png/kaplan_meyer_mod_severe.png", kaplan_meyer_esque, height=4, width=3.5, bg="white", dpi=444)


combo_plot <- cowplot::plot_grid(kaplan_meyer_esque, kaplan_meyer_esque2, nrow=1)


ggsave("~/postdoc/goncalves_modelling/png/combo_kaplan_meyer.png", combo_plot, height=5, width=8, bg="white", dpi=444)






library(dplyr)


data <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/mic_drop_enrolment_march23.csv", header = TRUE, stringsAsFactors = FALSE)
data <- data %>%
  mutate(cumulative=cumsum(enrolled))

ggplot(data, aes(x=Week, y=cumulative))+
  geom_point(shape=21, fill="darkred")+
  geom_line()+
  ylab("Cumulative Enrollment")+
  theme_minimal()

data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDropSpecimensJan2023_withClinical.dta")

summary <- data %>%
  filter(immunol==1)%>%
  group_by(id)%>%
  summarise("n_immuno"=n())