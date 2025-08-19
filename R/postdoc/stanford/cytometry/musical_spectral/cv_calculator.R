meta30_freq_files <- list.files("~/Downloads", pattern="_freqs.csv", full.names = T)

list_of_dfs <- data.frame()

for(freq_file in meta30_freq_files){
  
  tmp = read.csv(freq_file)
  colnames(tmp)= c("File", "experiment", paste("meta_thirty_norm", seq(1,30), sep=""))
  tmp$workflow = basename(freq_file)
  
  list_of_dfs = bind_rows(list_of_dfs, tmp)
  
}

# old
# cv_summary <- list_of_dfs%>%
#   pivot_longer(cols=starts_with("meta_thirty"), names_to = "cluster", values_to = "freq")%>%
#   group_by(workflow, cluster, experiment)%>%
#   summarise("mean"=mean(freq, na.rm=T),
#             "sd"=sd(freq, na.rm = T),
#             "cv"=mean/sd)%>%
#   group_by(workflow, experiment)%>%
#   summarise("mean_cv"=mean(cv))

mean_summary <- list_of_dfs%>%
  pivot_longer(cols=starts_with("meta_thirty"), names_to = "cluster", values_to = "freq")%>%
  group_by(workflow, cluster, experiment)%>%
  summarise("mean"=mean(freq+0.0001, na.rm=T))


mean_summary%>%
  filter(workflow=="raw_freqs.csv")%>%
  ggplot(., aes(x=experiment, y=mean, fill = experiment))+
  geom_boxplot()+
  facet_wrap(~cluster)+
  theme_minimal()

cv_summary_experiment <- mean_summary%>%
  group_by(workflow, experiment)%>%
  summarise("mean_across_experiment"=mean(mean, na.rm=T),
            "sd_across_experiment"=sd(mean, na.rm = T),
            "cv_across_experiment"=sd_across_experiment/mean_across_experiment)

cv_summary <- mean_summary%>%
  group_by(workflow, cluster)%>%
  summarise("mean_across_workflow"=mean(mean, na.rm=T),
            "sd_across_workflow"=sd(mean, na.rm = T),
            "cv_across_workflow"=sd_across_workflow/mean_across_workflow)

cv_summary%>%
  ggplot(., aes(x=mean_across_workflow, y=cv))+
  geom_point()+
  facet_wrap(~workflow)+
  theme_minimal()


cv_summary%>%
  ggplot(., aes(x=workflow, y=cv_across_workflow, fill=workflow))+
  geom_violin(draw_quantiles = seq(0,1,0.25), color="black", size=0.21)+
  viridis::scale_fill_viridis(discrete = T, option="D")+
  theme_minimal()+
  theme(legend.position = "none")

# workflow       `sum(cv_across_workflow > 0.3)`
# <chr>                                    <int>
# 1 raw_freqs.csv                               23
# 2 w4m4_freqs.csv                              16
# 3 w4m6_freqs.csv                              15
# 4 w5m4_freqs.csv                              13
# 5 w5m6_freqs.csv                              12

cv_summary2 <- cv_summary%>%
  group_by(workflow)%>%
  summarise("mean_cv"=mean(cv_across_workflow))


  

  
