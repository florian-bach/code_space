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


# Stacked Barchart ####

sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1


lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "pink")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")

summary$timepointf <- factor(summary$timepoint, levels=c("Baseline", "DoD", "T6", "C45"))


colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")

col_pal <- colcsv$colour
names(col_pal) <- colcsv$cluster_id

cluster_order <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_order.csv")$x


summary$n_infection <- ifelse(summary$volunteer %in% c("v313", "v315", "v320"), "first", "third")

summary <- summary %>%
  filter(cluster_id %in% sig_clusters) %>%
  filter(timepoint=="T6")
  
  activation_stacked_barchart <- ggplot(summary, aes(x=volunteer, y=frequency/100, fill=factor(cluster_id, levels=c(cluster_order))))+
    geom_bar(stat="identity", position="stack")+
    #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
    theme_minimal()+
    #facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
    ggtitle("Overall T cell activation")+
    scale_fill_manual(values=col_pal)+
    # scale_fill_manual(values=lineage_palette, labels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "NKT")))+
    scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
    theme(plot.title = element_text(hjust=0.5, size=11),
          strip.text = element_text(hjust=0.5, size=10, face = "bold"),
          axis.title.x = element_blank(),
          legend.position="none",
          axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid.minor.y = element_blank(),
          strip.placement = "outside")
  
  ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/t6_activation_cluster_id_stack.png", activation_stacked_barchart, height=4, width=4)






summary_t6 <- subset(summary, summary$timepoint=="T6")

summary_t6$scaled_freq <- ifelse(summary_t6$n_infection=="third",summary_t6$frequency/6, summary_t6$frequency/3)

sig_summary_t6 <- subset(summary_t6, summary_t6$cluster_id %in% sig_t6_clusters)




sig_ter_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)
sig_prim_t6_clusters <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/prim_t6_df_edger.csv", header = TRUE, stringsAsFactors = FALSE)

sig_ter_t6_clusters <- subset(sig_ter_t6_clusters, sig_ter_t6_clusters$p_adj<0.05 & sig_ter_t6_clusters$logFC>1)$cluster_id
sig_prim_t6_clusters <- subset(sig_prim_t6_clusters, sig_prim_t6_clusters$p_adj<0.05 & sig_prim_t6_clusters$logFC>1)$cluster_id

write.table(unique(c(sig_prim_t6_clusters,sig_ter_t6_clusters)), "prim_ter_sig_clusters.txt", sep="\t", col.names = FALSE, row.names = FALSE)


prim_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="first" & summary_t6$cluster_id %in% sig_prim_t6_clusters)


cluster_order <- prim_summary_t6 %>%
  group_by(cluster_id) %>%
  mutate("sum"=sum(frequency)) %>%
  arrange(desc(sum)) %>%
  select(cluster_id)

cluster_order <- rev(unique(cluster_order$cluster_id))

write.csv(cluster_order, "cluster_order.csv", row.names = FALSE)

prim_radial_plot_label_df <- data.frame(xpos=length(unique(prim_summary_t6$cluster_id))*0.1, ypos=seq(1,5)/100, label_text=paste(seq(1,5), "%", sep=''))

prim_sig_activation_polar <- ggplot()+
  geom_hline(prim_radial_plot_label_df, mapping=aes(yintercept = ypos), color="darkgrey")+
  geom_bar(data=prim_summary_t6,aes(x=factor(cluster_id, levels=cluster_order), y=scaled_freq/100, fill=factor(cluster_id, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  shadowtext::geom_shadowtext(data=prim_radial_plot_label_df, mapping=aes(x=xpos, y=ypos, label=label_text), size=2.5,  colour="black", bg.color="white")+
  ggtitle("differentially UP at T6 first")+
  scale_fill_manual(values=cluster_palette)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0,0.05))+
  coord_polar()+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11, vjust=-6),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- get_legend(prim_sig_activation_polar)

prim_sig_activation_polar <- prim_sig_activation_polar+theme(legend.position = "none")


ter_summary_t6 <- subset(summary_t6, summary_t6$n_infection=="third" & summary_t6$cluster_id %in% sig_prim_t6_clusters)

ter_summary_t6$fill <- ifelse(ter_summary_t6$cluster_id %in% sig_ter_t6_clusters, ter_summary_t6$cluster_id, "blank")

#cluster_palette <- c(cluster_palette, "blank"=NA)

ter_radial_plot_label_df <- data.frame(xpos=length(unique(ter_summary_t6$cluster_id))*0.1, ypos=seq(1,5)/100, label_text=paste(seq(1,5), "%", sep=''))

ter_sig_activation_polar <- ggplot()+
  geom_hline(ter_radial_plot_label_df, mapping=aes(yintercept = ypos), color="darkgrey")+
  geom_bar(data=ter_summary_t6,aes(x=factor(cluster_id, levels=cluster_order), y=scaled_freq/100, fill=fill), stat="identity", position="stack", width=1)+
  theme_minimal()+
  shadowtext::geom_shadowtext(data=prim_radial_plot_label_df, mapping=aes(x=xpos, y=ypos, label=label_text), size=2.5,  colour="black", bg.color="white")+
  ggtitle("differentially UP at T6 third")+
  scale_fill_manual(values=cluster_palette, na.translate=TRUE, na.value=NA)+
  scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0,0.05))+
  coord_polar()+
  theme(plot.title = element_text(hjust=0.5, size=11, vjust=-6),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,1,1), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())
  


combo_radial <- plot_grid(prim_sig_activation_polar, ter_sig_activation_polar)

combo_radial <- plot_grid(combo_radial, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")
ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_activation_stacked_polar.png", combo_radial, height=4, width=6.5)

#ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/sig_activation_stacked_polar.png", sig_activation_polar, height=4, width=9)

#non circular alternatives


dodged_lollipop_bar <- ggplot(summary, aes(y=volunteer, x=frequency/100, fill=factor(cluster_id, levels=c(cluster_order))))+
  geom_bar(stat="identity", position="dodge")+
  #geom_text(aes(y=(total_cd3/100)+0.01, label=paste0(round(total_cd3, digits = 1), "%", sep='')))+
  theme_minimal()+
  #facet_wrap(~timepointf, strip.position = "bottom", ncol=4)+
  ggtitle("Overall T cell activation")+
  scale_fill_manual(values=col_pal)+
  # scale_fill_manual(values=lineage_palette, labels=rev(c("CD4", "Treg", "CD8", "MAIT", expression(paste(gamma, delta)), "DN", "NKT")))+
  #scale_y_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.placement = "outside")



prim_sig_activation_lolli <- ggplot()+
  geom_bar(data=prim_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(cluster_id, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 first")+
  scale_fill_manual(values=col_pal)+
  scale_x_reverse(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0.05, 0), breaks=seq(1,5)/100, expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- get_legend(prim_sig_activation_lolli)

prim_sig_activation_lolli <- prim_sig_activation_lolli+theme(legend.position = "none")


ter_sig_activation_lolli <-   ggplot()+
  geom_bar(data=ter_summary_t6,aes(y=factor(cluster_id, levels=cluster_order), x=scaled_freq/100, fill=factor(fill, levels=cluster_order)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("differentially UP at T6 third")+
  scale_fill_manual(values=col_pal)+
  scale_x_continuous(name = "Percentage of CD3+ T cells\n", labels=scales::percent_format(accuracy = 1), limits = c(0,0.05), breaks=seq(1,5)/100 ,expand = c(0,0))+
  guides(fill=guide_legend(reverse = TRUE, nrow = 5,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,5,1,0), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())



combo_lolli <- plot_grid(prim_sig_activation_lolli, ter_sig_activation_lolli)

combo_lolli <- plot_grid(combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")
ggsave("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/all_activation_lollipop.png", combo_lolli, height=4, width=7)




# 
# 