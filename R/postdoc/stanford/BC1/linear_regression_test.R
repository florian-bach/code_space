


lm_ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  mutate(id=as.character(id))%>%
  nest() %>%
  mutate(cell_incidence_model=map(data, ~lm(conc ~ incidence_value, data=.)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

lm_ab_incidence_sigs <- lm_ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)





poiss_ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  mutate(id=as.character(id))%>%
  nest() %>%
  # lme4 poisson model
  mutate(cell_incidence_model=map(data, ~glm(data=., incidence_value ~ conc, family="poisson")))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

poiss_ab_incidence_sigs <- poiss_ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)


 

nb_ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  mutate(id=as.character(id))%>%
  nest() %>%
  # lme4 nbon model
  mutate(cell_incidence_model=map(data, ~glm.nb(data=., incidence_value ~ conc)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

nb_ab_incidence_sigs <- nb_ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)





mass_nb_ab_incidence_purf <- ab_clin %>%
  group_by(antigen, incidence_type, timepoint)%>%
  filter(!duplicated(id), incidence_type!="symp_6_12")%>%
  mutate(id=as.character(id))%>%
  nest() %>%
  # lme4 nbon model
  mutate(cell_incidence_model=map(data, ~glm.nb(data=., incidence_value ~ conc)))%>%
  mutate(cell_incidence_model_summary=map(cell_incidence_model, ~summary(.)))%>%
  mutate(cell_incidence_model_summary_p=map_dbl(cell_incidence_model_summary, ~coef(.)[8]))%>%
  ungroup()%>%
  group_by(timepoint)%>%
  mutate(cell_incidence_model_summary_padj= p.adjust(cell_incidence_model_summary_p))

mass_nb_ab_incidence_sigs <- mass_nb_ab_incidence_purf %>%
  dplyr::select(antigen, incidence_type, timepoint, cell_incidence_model_summary_padj)%>%
  filter(cell_incidence_model_summary_padj<0.1)



ab_poisson_incidence_sigs=poiss_ab_incidence_sigs



list_of_ab_poisson_plots <- list()
for(i in 1:nrow(poiss_ab_incidence_sigs)){
  
  plot_data <- ab_clin %>%
    # filter(!is.na(antigen), timepoint==3) %>%
    filter(antigen==paste(ab_poisson_incidence_sigs[i,1]) & incidence_type==paste(ab_poisson_incidence_sigs[i,2]) & timepoint==paste(ab_poisson_incidence_sigs[i,3]))%>%
    filter(!duplicated(id))
  
  plot <- ggplot(plot_data, aes(x=factor(incidence_value), y=conc, fill=factor(incidence_value)))+
    geom_boxplot()+
    geom_point(shape=21)+
    ggtitle(paste(unique(ab_poisson_incidence_sigs[i,3])))+
    scale_fill_manual(values=incidence_cols)+
    # facet_wrap(~timepoint)+
    # geom_point(fill=pc1_cols[i], alpha=1, shape=21)+
    # geom_smooth(method="lm")+
    # scale_y_log10()+
    # ggpubr::stat_cor(method = "spearman", label.y = -1.5, na.rm = TRUE, size=2)+
    ylab(ab_poisson_incidence_sigs[i,1])+
    xlab(ab_poisson_incidence_sigs[i,2])+
    theme_minimal()+
    theme(legend.position = "none")
  list_of_ab_poisson_plots[[i]] <- plot
}

sig_ab_poisson_plot <- cowplot::plot_grid(plotlist = list_of_ab_poisson_plots)
