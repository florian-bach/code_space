first_deg <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/first_DE.csv")
third_deg <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/third_DE.csv")

sig_first <- subset(first_deg, p_val_adj<0.05)
sig_third <- subset(third_deg, p_val_adj<0.05)

table(sig_first$cell_type)
table(sig_third$cell_type)

old_cluster_counts <- read.csv("~/Downloads/final_cluster_counts.csv")
new_cluster_counts <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/big_data_cluster_counts.csv")



keep <- "Third"
baseline <- subset(old_cluster_counts, Timepoint=="Baseline" & N_Infection == keep)
T6 <- subset(old_cluster_counts, Timepoint=="T6" & N_Infection == keep)


old_prim_lolli <- ggplot()+
  geom_bar(data=baseline, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage, fill=factor(Volunteer)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("Baseline")+
  scale_x_reverse(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0),limits=c(0.3, 0))+
  guides(fill=guide_legend(label.position = "top", reverse = FALSE, nrow = 1,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- cowplot::get_legend(old_prim_lolli)

old_prim_lolli <- old_prim_lolli+theme(legend.position = "none")


ter_lolli <-   ggplot()+
  geom_bar(data=T6, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage, fill=factor(Volunteer)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("T6")+
  scale_x_continuous(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(0,0.3))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())


old_combo_lolli <- cowplot::plot_grid(old_prim_lolli, ter_lolli)
old_combo_lolli <- cowplot::plot_grid(old_combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")

# title <- cowplot::ggdraw()+cowplot::draw_label("First Infection", hjust = 0.5) 
title <- cowplot::ggdraw()+cowplot::draw_label(paste(keep, "Infection"), hjust = 0.5) 
old_combo_lolli <- cowplot::plot_grid(title, old_combo_lolli, ncol=1, rel_heights = c(0.1,1))

#ggsave(filename = "~/postdoc/edinburgh/scRNAseq/revision_figures/old_first_vol_lollipop.png", old_combo_lolli, width = 8, height = 5, bg="white", dpi=444)
ggsave(filename = "~/postdoc/edinburgh/scRNAseq/revision_figures/old_third_vol_lollipop.png", old_combo_lolli, width = 8, height = 5, bg="white", dpi=444)













keep <- "First"
baseline <- subset(new_cluster_counts, Timepoint=="Baseline" & N_Infection == keep)
T6 <- subset(new_cluster_counts, Timepoint=="T6" & N_Infection == keep)


new_prim_lolli <- ggplot()+
  geom_bar(data=baseline, aes(y=factor(Cluster_ID, levels = unique(baseline$Cluster_ID)), x=Percentage, fill=factor(Volunteer)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("Baseline")+
  scale_x_reverse(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0),limits=c(0.3, 0))+
  guides(fill=guide_legend(label.position = "top", reverse = FALSE, nrow = 1,keyheight = unit(2, "mm"), keywidth = unit(4, "mm")))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_blank())


cluster_leg <- cowplot::get_legend(new_prim_lolli)

new_prim_lolli <- new_prim_lolli+theme(legend.position = "none")


ter_lolli <-   ggplot()+
  geom_bar(data=T6, aes(y=factor(Cluster_ID, levels = unique(T6$Cluster_ID)), x=Percentage, fill=factor(Volunteer)), stat="identity", position="stack", width=1)+
  theme_minimal()+
  ggtitle("T6")+
  scale_x_continuous(name = "Percentage of CD4+ T cells\n", labels=scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(0,0.3))+
  theme(plot.title = element_text(hjust=0.5, size=11),
        strip.text = element_text(hjust=0.5, size=10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        plot.margin = unit(c(1,0,1,5), "mm"),
        legend.position = "none",
        legend.text = element_text(size=6),
        legend.title = element_blank())


new_combo_lolli <- cowplot::plot_grid(new_prim_lolli, ter_lolli)
new_combo_lolli <- cowplot::plot_grid(new_combo_lolli, cluster_leg, ncol=1, rel_heights = c(5,1), rel_widths = c(1,1,2),align="v", axis = "tbrl")

# title <- cowplot::ggdraw()+cowplot::draw_label("First Infection", hjust = 0.5) 
title <- cowplot::ggdraw()+cowplot::draw_label(paste(keep, "Infection"), hjust = 0.5) 
new_combo_lolli <- cowplot::plot_grid(title, new_combo_lolli, ncol=1, rel_heights = c(0.1,1))

ggsave(filename = "~/postdoc/edinburgh/scRNAseq/revision_figures/new_first_vol_lollipop.png", new_combo_lolli, width = 8, height = 5, bg="white", dpi=444)
#ggsave(filename = "~/postdoc/edinburgh/scRNAseq/revision_figures/new_third_vol_lollipop.png", new_combo_lolli, width = 8, height = 5, bg="white", dpi=444)
