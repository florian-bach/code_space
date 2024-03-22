library(dplyr)
library(tidyr)
library(ggplot2)

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
mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_10/MICDROP expanded database through October 31st 2023.dta")


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


pcr_dens <- ggplot(smear_positive, aes(x=factor(n_smear), y=qPCRparsdens, fill=factor(n_smear)))+
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
         )

comp_cases_plot <- malaria %>%
  arrange(desc(factor(mstatus)))%>%
  ggplot(., aes(x=n_infection, y=pardens))+
  geom_point(aes(color=factor(mstatus)))+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color="darkred")+
  # geom_line(alpha=0.3, aes(group=id))+
  # geom_boxplot(aes(group=factor(AGE)))+
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.8, color="white")+
  ggtitle("")+
  ylab("qPCR Parasite Density")+
  xlab("n_infection")+
  scale_y_log10(limits=c(0.9, 10^6), breaks=c(10^0, 10^2, 10^4, 10^6), labels=scales::label_number())+
  scale_x_continuous(limits = c(0,10), breaks = seq(1,10))+
  annotation_logticks(sides = "l", )+
  theme_minimal()+
  # facet_wrap(~mstatus)+
  theme(legend.title = element_blank())+
  scale_color_manual(values = c(comp_pal, "quinine"="purple"))

comp_df <- malaria %>%
  group_by(n_infection)%>%
  summarise(comp_cases=sum(mstatus=="complicated"),
            comp_rate=comp_cases/n(),
            "total_cases"=sum(mstatus!="complicated"))


comp_model <- glm(comp_rate~n_infection+I(n_infection^2), weights = total_cases, data=comp_df, family = "binomial")

comp_model_fun <- function(x){
  exp(comp_model$coefficients[1])*
    exp(comp_model$coefficients[2])^x*
    exp(comp_model$coefficients[3])^x^2}


comp_model_plot <- ggplot(comp_df, aes(x=n_infection, y=comp_rate))+
  geom_point(color="darkred")+
  scale_y_continuous(limits = c(0,0.2), labels = scales::label_percent())+
  scale_x_continuous(limits = c(1,8), breaks = seq(1,8))+
  ylab("modeled risk of\ncomplicated malaria")+
  xlab("n_infection")+
  geom_ribbon(data=model_visualiser(comp_model, "n_infection"),
              aes(x=n_infection, ymin = exp(lci), ymax = exp(uci)),
              alpha = 0.2, inherit.aes = FALSE)+
  geom_function(fun = comp_model_fun, colour="black")+
  theme_minimal()
  

visreg(comp_model, xvar="n_infection", xlim=c(1,7), scale="response")


mic_drop_comp <- comp_cases_plot + comp_model_plot + plot_annotation(
    title = 'Complicated Disease in MIC-DRoP ',
    subtitle = 'Data inclusive of 10/2023')+
  theme(plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/complicated/mic_drop_comp.png", mic_drop_comp, width=8, height=4, bg="white", dpi=444)


para_temp_plot <- complicated_data %>%
  filter(total_n_infection>4, n_infection < 15)%>%
  ggplot(., aes(x=n_infection, y=parsdens, group=id))+
  geom_point(aes(color=factor(round(temp))))+
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

ggsave("~/postdoc/stanford/clinical_data/PROMOTE/figures/para_temp_plot.png", para_temp_plot, width=11, height=8, bg="white", dpi=444)



# sandbox ####

labels <- list(colnames(mic_drop))
# labels <- list(vector(length = ncol(mic_drop)))

for(i in 1:ncol(mic_drop)){
  entry <- attr(mic_drop[[colnames(mic_drop[i])]], "label")
  labels[i] <- entry
}

names(labels) <- colnames(mic_drop)

grep("*resp*", unlist(labels), value = TRUE)
