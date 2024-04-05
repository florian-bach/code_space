library(dplyr)
library(ggplot2)
library(patchwork)

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

comp_pal <- c("asymptomatic"="lightgrey",
              "uncomplicated"="black",
              "complicated"="orange",
              "severe"="darkred")


# mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDROP expanded database through December 31st 2022.dta")
mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2024_01/MICDROP expanded database through January 31st 2024.dta")


smear_positive <- mic_drop %>%
  filter(BSdich==1)%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_positive_smears") %>%
  arrange(compage) %>%
  mutate(n_smear = seq(1, max(total_positive_smears)))
  

n_infection_cols <- c("white", colorspace::sequential_hcl(n=max(smear_positive$n_smear), palette = "Lajolla"))

# overview plots ####

micro_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=pardens, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])


pcr_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=qPCRparsdens+0.0001, fill=factor(n_smear)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(n_smear)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Order of Positive Smear")+
  scale_y_log10()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values=n_infection_cols[-1])

dens_n_smear_combo <- cowplot::plot_grid(micro_dens, pcr_dens,  nrow = 1)






micro_age <- ggplot(smear_positive, aes(x=AGE, y=pardens, fill=factor(AGE)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("Microscopic Parasite Density")+
  xlab("Age (Weeks)")+
  geom_smooth()+
  # geom_vline(xintercept = 8)+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE, direction = -1)


pcr_age <- ggplot(smear_positive, aes(x=AGE, y=qPCRparsdens, fill=factor(AGE)))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3, aes(group=id))+
  geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ylab("qPCR Parasite Density")+
  xlab("Age (Weeks)")+
  scale_y_log10(limits=c(0.01, 10^6))+
  theme_minimal()+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(discrete = TRUE)

dens_age_combo <- cowplot::plot_grid(micro_age, pcr_age,  nrow = 1)


# n_infection ####
malaria <- mic_drop %>%
  filter(!is.na(mstatus) & mstatus !=0)%>%
  mutate(id=factor(id))%>%
  group_by(id) %>%
  add_count(name="total_malaria_episodes") %>%
  arrange(date) %>%
  mutate(n_infection = seq(1, max(total_malaria_episodes)))%>%
  mutate(mstatus = case_match(mstatus,
                              0~"asymptomatic",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"uncomplicated",
                              4~"uncomplicated"),
         )%>%
  mutate(comp_dich=if_else(mstatus=="complicated", 1, 0))

comp_cases_n <- malaria %>%
  arrange(desc(factor(mstatus)))%>%
  ggplot(., aes(x=n_infection, y=pardens))+
  geom_point(aes(color=factor(mstatus)))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  # geom_line(alpha=0.3, aes(group=id))+
  # geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ggtitle("")+
  ylab("qPCR Parasite Density")+
  xlab("Order of Infection")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  scale_x_continuous(limits = c(0,max(malaria$n_infection)), breaks = seq(1,max(malaria$n_infection)))+
  annotation_logticks(sides = "l", )+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.position="none")+
  scale_color_manual(values = c(comp_pal, "quinine"="purple"))

bins <- seq(0, 104, by=4)
nice_labels <- list(vector(length = length(bins)))

for(i in 1:(length(bins)-1)){
  nice_labels[[i]]=paste(bins[i], "to", bins[i+1])
}

# age ####

(comp_cases_age <- malaria %>%
  arrange(desc(factor(mstatus)))%>%
  mutate(agebins=cut(AGE, breaks = seq(0, 104, by=4), labels = nice_labels))%>%
  ggplot(., aes(x=factor(agebins), y=pardens+0.001))+
  geom_vline(xintercept = c("24 to 28",
                            "48 to 52",
                            "76 to 80"), linetype="dashed", color="grey")+
  geom_point(aes(color=factor(mstatus)))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  ggtitle("")+
  ylab("qPCR Parasite Density")+
  xlab("Age in weeks")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  # scale_x_continuous(limits = c(0,max(malaria$n_infection)), breaks = seq(1,max(malaria$n_infection)))+
  annotation_logticks(sides = "l", )+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_manual(values = c(comp_pal, "quinine"="purple"))
)




# linear  regression ####
comp_df <- malaria %>%
  group_by(n_infection)%>%
  summarise(comp_cases=sum(mstatus=="complicated"),
            comp_rate=comp_cases/n(),
            "total_cases"=sum(mstatus!="complicated"))


comp_model <- glm(comp_rate~n_infection+I(n_infection^2), weights = total_cases, data=comp_df, family = "binomial")
comp_model_integer <- glm(comp_dich~n_infection+I(n_infection^2), data=malaria, family = "binomial")


comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}


comp_model_plot <- ggplot(comp_df, aes(x=n_infection, y=comp_rate))+
  geom_point(color="darkred")+
  scale_y_continuous(limits = c(0,0.15), labels = scales::label_percent())+
  scale_x_continuous(limits = c(1,9), breaks = seq(1,9))+
  ylab("risk of complicated malaria")+
  xlab("n_infection")+
  geom_ribbon(data=model_visualiser(comp_model, "n_infection"),
              aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  theme_minimal()
  



# combo plot ####
mic_drop_comp <- comp_cases_age + comp_cases_n + comp_model_plot + plot_annotation(
    title = 'Complicated Disease in MIC-DRoP ',
    subtitle = 'Data inclusive of 01/2024')+
  theme(plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/complicated/mic_drop_comp_20240131.png", mic_drop_comp, width=15, height=5, bg="white", dpi=444)




# sandbox ####
# 
# labels <- list(colnames(mic_drop))
# # labels <- list(vector(length = ncol(mic_drop)))
# 
# for(i in 1:ncol(mic_drop)){
#   entry <- attr(mic_drop[[colnames(mic_drop[i])]], "label")
#   labels[i] <- entry
# }
# 
# names(labels) <- colnames(mic_drop)
# 
# grep("*resp*", unlist(labels), value = TRUE)



latecomers <- malaria %>%
  filter(n_infection==1 & AGE>52)%>%
  summarise(id)

malaria %>%
  arrange(desc(factor(mstatus)))%>%
  mutate(agebins=cut(AGE, breaks = seq(0, 104, by=4), labels = nice_labels))%>%
  mutate(arm=if_else(id %in% latecomers$id, "IPT12", "unknown"))%>%
  ggplot(., aes(x=agebins, y=n_infection))+
  geom_vline(xintercept = c("24 to 28",
                            "48 to 52",
                            "76 to 80"), linetype="dashed", color="grey")+
  geom_point(aes(color=factor(arm)), position = position_jitter(width=0.1))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  # ggtitle("")+
  ylab("Order of Infection")+
  xlab("Age in weeks")+
  scale_y_continuous(limits = c(0,max(malaria$n_infection)), breaks = seq(1,max(malaria$n_infection)))+
  # scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  # scale_x_continuous(limits = c(0,max(malaria$n_infection)), breaks = seq(1,max(malaria$n_infection)))+
  # annotation_logticks(sides = "l", )+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_manual(values=c("darkgrey", "black"))




para_temp_plot <- malaria %>%
  filter(total_malaria_episodes>3)%>%
  mutate(arm=if_else(id %in% latecomers$id, "IPT12", "unknown"))%>%
  arrange(arm)%>%
  ggplot(., aes(x=n_infection, y=pardens, group=id))+
  geom_point(aes(shape=arm, color=factor(round(temp))))+
  geom_line(alpha=0.2)+
  facet_wrap(~id, scales = "free_x")+
  ggtitle("")+
  ylab("qPCR Parasite Density")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  scale_x_continuous(limits = c(1,15))+
  theme_minimal()+
  geom_vline(xintercept = 37.5)+
  # facet_wrap(~mstatus)+
  theme(legend.title = element_blank())+
  scale_color_manual(values=c("darkgreen", "green", "yellow", "orange", "red", "darkred"))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/complicated/para_temp_plot.png", para_temp_plot, width=11, height=16, bg="white", dpi=444)

# ipt12 comparison 
# malaria %>%
#   filter(total_malaria_episodes>3)%>%
#   mutate(arm=if_else(id %in% latecomers$id, "IPT12", "unknown"))%>%
#   arrange(desc(arm))%>%
#   ggplot(., aes(x=factor(n_infection), y=pardens+0.01))+
#   # geom_point(aes(shape=arm, color=factor(round(temp))))+
#   geom_boxplot(aes(fill=arm))+
#   # geom_line(alpha=0.2)+
#   ggtitle("")+
#   ylab("qPCR Parasite Density")+
#   scale_y_log10(limits=c(0.1, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
#   theme_minimal()+
#   # facet_wrap(~mstatus)+
#   theme(legend.title = element_blank())
# 
# malaria %>%
#   # filter(total_malaria_episodes>3)%>%
#   mutate(arm=if_else(id %in% latecomers$id, "IPT12", "unknown"))%>%
#   group_by(arm, n_infection)%>%
#   summarise("n"=n())%>%
#   print(n=25)
# 
# 
# malaria %>%
#   filter(n_infection <= 4)%>%
#   mutate(arm=if_else(id %in% latecomers$id, "IPT12", "unknown"))%>%
#   arrange(desc(arm))%>%
#   ggplot(., aes(x=factor(n_infection), y=temp))+
#   geom_point(aes(color=arm), position=position_dodge(width=0.75))+
#   geom_boxplot(aes(fill=arm))+
#   # geom_line(alpha=0.2)+
#   ggtitle("")+
#   ylab("qPCR Parasite Density")+
#   theme_minimal()+
#   # facet_wrap(~mstatus)+
#   theme(legend.title = element_blank())
