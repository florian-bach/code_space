malaria <- read.csv("~/Downloads/malaria_report_2022_cases_deaths.csv")

# long_malaria <- malaria %>%
  
malaria%>%
  filter(Year>2012)%>%
ggplot(aes(x=Year))+
  geom_point(aes(y=Deaths), color="darkred")+
  geom_point(aes(y=Cases), color="darkblue")+
  geom_line(aes(y=Deaths), color="darkred")+
  geom_line(aes(y=Cases), color="darkblue")+
  scale_y_continuous(
    # Features of the first axis
    name = "Deaths (Thousands)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1, name="Cases (Millions)")
  ) + 
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "darkblue", size=13),
    axis.title.y.right = element_text(color = "darkred", size=13)
  )


mal_mort_plot <- malaria%>%
  filter(Year>2012)%>%
  ggplot(aes(x=Year))+
  geom_line(aes(y=Deaths*1000), color="darkred", linewidth=1.5, alpha=0.7)+
  geom_point(aes(y=Deaths*1000), fill="darkred", color="white", shape=21, size=4)+
  theme_minimal()+
  ggtitle("malaria deaths / year\n")+
  scale_y_continuous(labels = scales::label_comma())+
  scale_x_continuous(breaks = seq(2013, 2021))+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=15),
        #axis.title = element_text(size=14),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/funding_applications/stanford_embl/mal_mort_plot.png", mal_mort_plot, width = 7, height=4, dpi=444, bg="white")
  