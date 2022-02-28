
goncalves_data <- data.frame("N_Infection"=seq(0,14),
                             "Severe"=c(0, 32, 18, 20, 9, 9, 2, 2, 8, NA, 1, NA, NA, NA, 1),
                             "Infections"=c(0, 715, 540, 424, 330, 267, 220, 168, 131, NA, 90, NA, NA, NA, 25),
                             "Susceptible"=c(102, 70, 52, 32, 23, 14, 12, 10, 2, 2, 1, 1, 1, 1, 0),
                             "Severe_Moderate"=c(0, 76, 41, 31, 15, 12, 5, 8, 5, 1, 1, 3, NA, NA, 1),
                             "Supp_Infections"=c(0, 715, 501, 369, 276, 213, 169, 120, 86, 62, 55, 47, NA, NA, 25))

goncalves_data <- goncalves_data%>%
  replace(is.na(.), 0)

# <- data.frame("Model", "P1", "SE(P1)", "alpha1", "SE(alpha1)",  "alpha2", "SE(alpha2)","logL", "AIC" )



n_successes <- goncalves_data$Severe_Moderate[2:12]
m_trials <- goncalves_data$Supp_Infections[2:12]


brute_force <- list()
for(i in 1:10){
  
  print(paste("model with kink after", i))
  
  nick_decay_likelihood <- function(p_i, alpha_1, alpha_2){
   #p_i=0.1; alpha_1=0.005; alpha_2=0.001
   
    #define stepwise function
    p_vec_1 <- p_i - alpha_1 * 0:10
    p_vec_2 <- p_i - alpha_1 * 0:10 - alpha_2 * (1:11-i)
    
    #combine based on nick position
    p_vec <- c(p_vec_1[1:i], p_vec_2[(i+1):length(n_successes)])
    #calc distribution
    res <- sum(
      dbinom(x = n_successes, size = m_trials, prob = p_vec, log = TRUE)
      )*-1
    return(res)}
  

  #print(nick_decay_likelihooderator(p_i=0.1, alpha_1=0.005, alpha_2=0.001))
  nick_decay_estimate <- stats4::mle(minuslog = nick_decay_likelihood,
                                     start = list(p_i=0.2,
                                                  alpha_1=0.001,
                                                  alpha_2=0.001)
                                     )
  
  brute_force <- c(brute_force, nick_decay_estimate)
  }

lapply(brute_force, function(x)AIC(x))

