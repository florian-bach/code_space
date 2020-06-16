#plot shtuff####
smol_time12_level <- ggplot(big_table, aes(x=UMAP1, y=UMAP2))+
    
    #facet-specfic min/max in color scheme
    #stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=12)+
    
    #same min/max across all facets in color scheme  
    #stat_density_2d(aes(fill = after_stat(level)), geom="polygon")+
    geom_point(shape="o", color="black")+
    stat_density_2d(contour=T, bins=n_bins, aes(color=after_stat(level)))+
    scale_color_viridis(option = "A")+
    xlim(c(-13, 10))+
    ylim(c(-10, 10))+
    theme_minimal()+
    facet_grid(.~timepoint)+
    theme(strip.text = element_text(size=14))



smol_time12_nlevel <- ggplot(big_table, aes(x=UMAP1, y=UMAP2, color=after_stat(nlevel)))+
  
  #facet-specfic min/max in color scheme
  #stat_density_2d(aes(fill = after_stat(nlevel)), geom="polygon", bins=12)+
  
  #same min/max across all facets in color scheme  
  stat_density_2d(contour=T, bins=n_bins, aes(color=after_stat(level)))+
  geom_point(shape="o", color="black")+
  scale_color_viridis(option = "A")+
  xlim(c(-13, 10))+
  ylim(c(-10, 10))+
  theme_minimal()+
  facet_grid(.~timepoint)+
  theme(strip.text = element_text(size=14))

#combine shtuff ####

plot_grid(smol_time12_level, smol_time12_nlevel, ncol=1)

plot_nlevel <- ggplot_build(smol_time12_level)
plot_level <- ggplot_build(smol_time12_nlevel)

head(plot_nlevel$data[[1]])
head(plot_level$data[[1]])

