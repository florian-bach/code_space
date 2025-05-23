head(cluster_freqs)


dummy_data <- subset(cluster_freqs, cluster_freqs$cluster_id=="resting  Vd2+ ")


#dummy_data <- subset(cluster_freqs, cluster_freqs$cluster_id==unique(cluster_freqs$cluster_id)[31])

simple_lm <- (ggplot(dummy_data,aes(x=as.numeric(factor(timepoint)), y=count))+
  geom_point(aes(color=volunteer), size = 3)+
  theme_minimal()+
  xlab("timepoint")+
  geom_smooth(formula = y~x, data=dummy_data, method = lm, se=FALSE))

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/simple_lm.png", simple_lm)       


p0 <- ggplot(dummy_data,aes(x=as.numeric(factor(timepoint)), y=count, group=volunteer))+
                geom_point(aes(color=volunteer), size = 3)+
                theme_minimal()+
                xlab("timepoint")

plm <- glm(formula = count ~ timepoint + volunteer, data=dummy_data)

simple_glm_vol <- p0+geom_smooth(method = "glm", mapping=aes(y=predict(plm,dummy_data), color=volunteer), se=FALSE)


mixed_glm_vol <- p0+geom_smooth(method = "glm", se=0, aes(color=volunteer))


ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/simple_glm_vol.png", simple_glm_vol)       

ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/mixed_glm_vol.png", mixed_glm_vol)






p1 <- ggplot(cluster_freqs,aes(x=as.numeric(factor(timepoint)), y=frequency, group=volunteer))+
  geom_point(size = 3, aes(shape=volunteer, color=timepoint))+
  theme_minimal()+
  xlab("timepoint")+
  geom_smooth(aes(x=as.numeric(factor(timepoint)), y=frequency, group=volunteer), data=cluster_freqs, method="lm", se = FALSE, formula= y~x+volunteer)



p1 <- ggplot(cluster_freqs, aes(x=factor(timepoint), y=count, group=volunteer))+
  geom_point(size = 3, aes(shape=volunteer, color=timepoint))+
  theme_minimal()+
  xlab("timepoint")+
  facet_wrap(~cluster_id, labeller=label_wrap_gen(), ncol=8, scales="free")

#p1+geom_smooth(aes(group=volunteer), method = "nls", se=FALSE, formula = y~a*x+b,  method.args = list(start = list(a = 0.1, b = 0.1)))

  p1+geom_smooth(aes(group=volunteer), method = "auto", se=FALSE, formula = y~x+b)


  
#gamma delta clustering busines stuff et cetera
  
  
  
  subset(merging_table1$new_cluster, grepl("^45$", merging_table1$old_cluster, fixed=F))
  
  vd2 <- as.character(subset(merging_table1$new_cluster, grepl("^35$", merging_table1$old_cluster, fixed=F)))
  
  
  vd2_old <- subset(merging_table1$old_cluster, grepl(vd2, merging_table1$new_cluster, fixed=T))
  
  
  play <- da_t6_count_box$data
  
  play <- subset(play, play$cluster_id %in% c(vd2_old))
  #play <- subset(play, play$timepoint %in% c("Baseline", "T6"))
  
  gamma_delta_count_plot <- ggplot(play, aes(x=timepoint, y=count))+
    geom_boxplot(aes(fill=timepoint))+
    geom_point(aes(shape=volunteer))+
    theme_minimal()+
    facet_wrap(~cluster_id, scales="free")+
    scale_y_log10()
  
  
  ggsave("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/figures_for_phil/gamma_delta_count_plot.png", gamma_delta_count_plot)
  
  
