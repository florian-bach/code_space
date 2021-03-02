library(dplyr)
library(tidyr)
library(ggplot2)


stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")

activated_clusters <- grep("activated", unique(stacked_bar_data$cluster_id), value=TRUE)

activated_barchart_data <- subset(stacked_bar_data, stacked_bar_data$cluster_id %in% activated_clusters)

activated_barchart_data$lineage <- substr(activated_barchart_data$cluster_id, nchar(activated_barchart_data$cluster_id)-2, nchar(activated_barchart_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg"), unique(activated_barchart_data$lineage))
activated_barchart_data$lineage <- stringr::str_replace_all(activated_barchart_data$lineage, lin_replacement)


summary <- activated_barchart_data %>%
  group_by(timepoint, cluster_id) %>%
  mutate(mean=base::mean(frequency), sd=stats::sd(frequency))

#tail(summary)

summary$lineage <- factor(summary$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))


#taken from khroma bright 6

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id


# Stacked Barchart Figure 1E ####


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")

summary$timepointf <- factor(summary$timepoint, levels=c("Baseline", "DoD", "T6", "C45"))

cluster_palette <- colcsv$colour
names(cluster_palette) <- colcsv$cluster_id

summary$n_infection <- ifelse(summary$volunteer %in% c("v313", "v315", "v320"), "first", "third")

activation_stacked_barchart <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=cluster_id))+
  geom_bar(stat="identity", position="stack")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=lineage_palette, labels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "NKT")))+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  coord_polar()+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")

ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/activation_stacked_barchart.png", activation_stacked_barchart, height=4, width=8)






summary_t6 <- subset(summary, summary$timepoint=="T6")

summary_t6$scaled_freq <- ifelse(summary_t6$n_infection=="third",summary_t6$frequency/6, summary_t6$frequency/3)

sig_summary_t6 <- subset(summary_t6, summary_t6$cluster_id %in% sig_t6_clusters)


sig_radial_plot_label_df <- data.frame(xpos=13, ypos=seq(1,5)/100, label_text=paste(seq(1,5), "%", sep=''))
all_radial_plot_label_df <- data.frame(xpos=18, ypos=seq(1,5)/100, label_text=paste(seq(1,5), "%", sep=''))


sig_activation_polar <- ggplot()+
  geom_bar(data=summary_t6,
           aes(x=cluster_id, y=scaled_freq/100, fill=cluster_id), stat="identity", position="stack")+
  theme_minimal()+
  facet_wrap(~n_infection, strip.position = "bottom", ncol=4)+
  # geom_text(sig_radial_plot_label_df, mapping=aes(x=xpos,
  #           y=ypos,
  #           label=label_text),
  #           size=2.5,
  #           colour="black")+
  geom_text(all_radial_plot_label_df, mapping=aes(x=xpos,
                                                  y=ypos,
                                                  label=label_text),
            size=2.5,
            colour="black")+
  ggtitle("Significant Activated T cells at T6")+
  #geom_text(data=all_radial_plot_label_df, aes(x=xpos, y=ypos, label=label_text), size=3, colour="darkgrey")+
  scale_fill_manual(values=cluster_palette)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  coord_polar()+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())
 
  

ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_activation_stacked_polar.png", sig_activation_polar, height=4, width=8)

#ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/sig_activation_stacked_polar.png", sig_activation_polar, height=4, width=9)





# 
# 