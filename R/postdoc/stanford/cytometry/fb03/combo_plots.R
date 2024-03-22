fb02_data <- read.csv("~/postdoc/stanford/cytometry/attune/FB02/fb02_long_data.csv")
fb02_data$experiment <- "fb02"

fb03_data <- read.csv("~/postdoc/stanford/cytometry/attune/FB03/fb03_long_data.csv")
fb03_data$experiment <- "fb03"

combo_data <- rbind(fb02_data, fb03_data[,match(colnames(fb02_data), colnames(fb03_data))])

combo_data$origin <- ifelse(combo_data$experiment=="fb02", "Uganda", combo_data$origin)
combo_data$serum <- factor(combo_data$serum, levels=c("autologous", "human AB", "FBS", "no serum"))

combo_data_medians <- combo_data %>%
  filter(!is.na(lineage),
         !is.na(origin),
         stim %in% c("PMA", "media"),
         marker_combo %in% c("CD69 CD40L", "OX40 CD137"),
         lineage %in% c("CD4", "CD4neg"),
         quadrant == "sum(Q)")%>%
  group_by(lineage, stim, serum, cell_pop, marker_combo, quadrant)%>%
  summarise("median"=round(median(freq), digits = 1))


cd4_cd4neg_summary_data <- combo_data %>%
  filter(!is.na(lineage),
         !is.na(origin),
         stim %in% c("PMA", "media"),
         marker_combo %in% c("CD69 CD40L", "OX40 CD137"),
         lineage %in% c("CD4", "CD4neg"),
         quadrant == "sum(Q)")



(combo_summary_plot <- ggplot(cd4_cd4neg_summary_data, aes(x=quadrant, y=(freq+0.0001)/100))+
    geom_point(aes(color=factor(id), group=serum, shape=origin), alpha=0.6, position = position_dodge(width=0.75))+
    geom_boxplot(aes(fill=serum), outlier.shape = NA)+
    # facet_grid(lineage+marker_combo~stim, scales = "free_y")+
    facet_wrap(lineage+marker_combo~stim, scales = "free_y", ncol = 4)+
    geom_text(aes(x=quadrant, y=1, label=median, group=serum), position = position_dodge(width=0.75), data = combo_data_medians, size=2.8)+
    scale_y_log10(labels = scales::label_percent())+
    theme_minimal()+
    guides(color=guide_legend(title="id"))+
    scale_fill_manual(values=serum_pal)+
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_blank(),
          legend.title = element_text(hjust=0.5)))

ggsave("~/postdoc/stanford/cytometry/attune/FB03/figures/fb02_fb03_combo_summary_plot.png", combo_summary_plot, width=10*1.5, height=8*1.5, bg="white")
