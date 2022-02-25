dummy_data <- data.frame("x"=1:100,
                         "y"=dbinom(x=1:100, size=100, prob = 0.5)
                         )

ggplot(dummy_data, aes(x=x, y=y))+
  geom_point()+
  theme_minimal()



n_successes <- goncalves_data$Severe_Moderate[2:9]
m_trials <- goncalves_data$Supp_Infections[2:9]


flat_likelihood <- function(p){
  dbinom(n_successes, m_trials, p)
}

nlm(flat_likelihood, p=0.06, stepmax = 0.5)





constant_p <- dbinom(n_successes, m_trials, seq(0.001, 0.1, by=0.001))
x_axis <-  seq(0.001, 0.1, by=0.001)

plot_data <- data.frame("y"=constant_p,
                        "x"=x_axis)

ggplot(plot_data, aes(x=x_axis, y=constant_p))+
  geom_point()+
  theme_minimal()+
  geom_smooth()



