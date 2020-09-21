# Panel A clinical presentation of volunteers undergoing vivax CHMI: plateltets, fever and adverse events ####


data <- read.csv("~/PhD/clinical_data/vac69a/symptoms_vac69a.csv", header = T, stringsAsFactors = F)
data$flo_timepoint <- gsub("am ", "am", data$flo_timepoint)
data$flo_timepoint <- gsub("pm ", "pm", data$flo_timepoint)
data$flo_timepoint <- gsub("DoD ", "Diagnosis", data$flo_timepoint)

long_data <- tidyr::gather(data, Symptom, Severity, colnames(data)[c(12, 14:ncol(data))])

long_data <- mutate(long_data, Volunteer=gsub("69010", "V", long_data$trial_number))

timepoint_levels <- list("Timepoints"=unique(long_data$flo_timepoint[gtools::mixedorder(long_data$flo_timepoint)]))

timepoint_levels$Timepoints[[length(timepoint_levels$Timepoints)+1]] <- timepoint_levels$Timepoints[32]

timepoint_levels <- timepoint_levels$Timepoints[-32]
# 
# myLoc <- (which(levels(long_data$flo_timepoint) == "DoD ") +
#             which(levels(long_data$flo_timepoint) == "C15 pm")) / 
#   2
# 
# symptom_heatmap <- ggplot(long_data, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=Symptom))+
#   geom_tile(aes(fill=factor(Severity), width=0.93, height=0.93), color=ifelse(grepl("DoD", long_data$flo_timepoint), "black", "lightgrey"))+
#   scale_fill_manual(values =  list("lightgrey", "yellow", "orange", "red"))+
#   facet_wrap(~Volunteer, scales="free")+
#   theme_minimal()+
#   guides(fill=guide_legend(title="Severity",
#                            override.aes = list(size = 0.5)))+
#   theme(axis.text.x = element_text(hjust=1, angle=45),
#         axis.title = element_blank(),
#         legend.title = element_blank())
# 
# ggsave("/home/flobuntu/PhD/clinical_data/vac69a/figures/symptom_heatmap.png", symptom_heatmap, height=8, width=14)

# make a figure for number of AEs per timepoint

library(dplyr)
library(tidyr)



long_data$flo_timepoint <- factor(long_data$flo_timepoint)
long_data$Volunteer <- tolower(long_data$Volunteer)

adverse_events <- long_data %>%
  filter(Severity > 0) %>%
  group_by(Volunteer, flo_timepoint, Severity) %>%
  summarise(ae_count = n())

fig1_theme <- theme(axis.title.x = element_blank(),
                    legend.title = element_text(size = 9), 
                    legend.text = element_text(size = 9),
                    axis.title=element_text(size=10))


(colored_stack <- ggplot(adverse_events,  aes(x=factor(flo_timepoint, levels=timepoint_levels), y=ae_count/6, fill=factor(Severity, levels=paste(rev(1:3)))))+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values =  list("1"="yellow", "2"="chocolate1", "3"="red"))+
  #facet_wrap(~Volunteer)+
  ylab("Mean # of AEs\nper Volunteer")+
  xlab("Timepoint")+
  ggtitle("Adverse Events")+
  geom_vline(aes(xintercept = 21.5))+
  scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_minimal()+
  fig1_theme+
  theme(axis.text.x = element_text(hjust=1, angle=45, size=5),
        plot.title = element_text(hjust=0.5),
        panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Severity",
                            override.aes = list(size = 0.1),
                           keywidth = 0.5,
                           keyheight = 0.5)))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/adverse_events_stacked.png", colored_stack, height=4, width=6)



fever <- subset(long_data, long_data$pyrexia_temp>37)

volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

fever_curves <- ggplot(fever, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=pyrexia_temp, color=Volunteer, group=Volunteer))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  ggtitle("Fever")+
  xlab("Timepoint")+
  ylab(expression(paste("Temperature (",degree,"C)",sep="")))+
  scale_y_continuous(limits=c(37.5, 40), breaks = seq(37.5, 40, by=0.5) )+
  theme_minimal()+
  fig1_theme+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8),
        legend.position = "none", 
        axis.title.x = element_blank())

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_curves.png", fever_curves)


# PLatelets/Lymphocytes


haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$trial_number <- gsub("69010", "v", haem_data$trial_number)

long_haem_data <- haem_data %>%
  # select(trial_number, timepoint, platelets, lymphocytes) %>%
  select(trial_number, timepoint, platelets) %>%
  filter(timepoint %in% c("_C_1", "_C1C7", "_EP", "_T6")) %>%
  gather(Cell, Frequency, c(platelets))

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

long_haem_data$Cell <- .simpleCap(long_haem_data$Cell)

long_haem_data$timepoint <- gsub("_C_1", "Baseline", long_haem_data$timepoint)
long_haem_data$timepoint <- gsub("_C1C7", "C7", long_haem_data$timepoint)
long_haem_data$timepoint <- gsub("_EP", "T1", long_haem_data$timepoint)
long_haem_data$timepoint <- gsub("_T6", "T6", long_haem_data$timepoint)


thrombos <- ggplot(long_haem_data, aes(x=timepoint, y=Frequency*1000, color=trial_number, group=trial_number))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=trial_number), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(expression(Cells~"/"~mu*L~blood))+
  guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
  ggtitle("Platelets")+
  fig1_theme+
  scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank())

# lgd1 <- get_legend(thrombos)
# 
thrombos <- thrombos+theme(legend.position = "none")
# 
# 
# lgd2 <- get_legend(thrombos)
# 
# thrombos <- thrombos+theme(legend.position = "none")




# Panel B Significant Plasma Analytes ####

library(tidyr)
library(dplyr)
library(ggplot2)

data3 <- read.csv("~/PhD/plasma/vac69a/big_plasma_table.csv")


data3[,3:ncol(data3)] <- log2(data3[,3:ncol(data3)])


data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T)) %>%
  mutate(timepoint = gsub("DoD", "Diagnosis", timepoint))

data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))


list_of_dfs_for_glm <- lapply(colnames(data3)[3:ncol(data3)], function(x) data.frame(select(data3, Volunteer, timepoint, x)))


list_of_models <- lapply(list_of_dfs_for_glm, function(x) glm(x[,3]~x[,2]+x[,1], data=x))


list_of_summaries <- lapply(list_of_models, function(x) cbind(summary(x)$coefficients, names(x$data)[3]))


df_of_model_results <- data.frame(do.call(rbind, list_of_summaries))
colnames(df_of_model_results) <- c("Estimate", "SE", "t_value", "raw_p", "Analyte")
df_of_model_results$Coefficient <- rownames(df_of_model_results)


df_of_model_results <- df_of_model_results[!grepl("Intercept", df_of_model_results$Coefficient),]
df_of_model_results <- df_of_model_results[!grepl("v0", df_of_model_results$Coefficient),]

df_of_model_results$p_adj <- p.adjust(as.numeric(as.character(df_of_model_results$raw_p)), method = "fdr")

sig_hits <- subset(df_of_model_results, df_of_model_results$p_adj<0.1)


sig_levels<- as.character(sig_hits[order(sig_hits$p_adj),]$Analyte)



siggy_hits <- sig_hits %>%
  group_by(Analyte) %>%
  top_n(n = -1, wt = p_adj)

# sig_levels <- as.character(siggy_hits[order(siggy_hits$p_adj),]$Analyte)
# # this next step is necessary because TGFbeta is undetectable/unchanged and top_n returns it three times because reasons
# sig_levels <- unique(sig_levels) 
# write.table(sig_levels, "analytes_sorted_by_padj.txt", row.names = FALSE, col.names = "Analyte")

long_data <- gather(data3, Analyte, Concentration, colnames(data3)[3:ncol(data3)])

long_data$timepoint <- gsub("DoD", "Diagnosis", long_data$timepoint)

sig_glm_data <- subset(long_data, long_data$Analyte %in% siggy_hits$Analyte)
sig_glm_data$AnalyteF <- factor(sig_glm_data$Analyte, levels=sig_levels)


sig_glm_data$p_adj <- sig_hits$p_adj[match(sig_glm_data$Analyte, sig_hits$Analyte)]

sig_glm_plot <- ggplot(sig_glm_data, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45")), y=2^Concentration, color=Volunteer, group=Volunteer))+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  facet_wrap(~ AnalyteF, scales = "free", ncol=4)+
  scale_y_continuous(trans = "log2", labels=scales::comma)+
  ylab("Plasma Concentration (pg / mL)")+
  theme_minimal()+
  guides(color=guide_legend(override.aes = list("size"=0.1)))+
  scale_color_manual(values=volunteer_palette)+
  fig1_theme+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1, angle=45, size=6),
        strip.background = element_rect(fill = "white", color = "white"))


ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/glm_sig_analytes_fdr_10-e-1.png", sig_glm_plot, width = 8, height=5)

clinical_graphs <- plot_grid(fever_curves, thrombos, colored_stack, ncol = 3, hjust=-2.2, vjust=0.8, rel_widths = c(0.8,1,1.7))
#ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/aes_and_fever_curves.png", clinical_graphs, height=3, width=8)

fig1 <- plot_grid(clinical_graphs, sig_glm_plot, nrow=2, ncol=1, hjust=-2.2, rel_heights = c(1.2,3))

ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/fig1_chunk.png", fig1, width = 8, height=6)

# Panel C  plasma heatmap ####

library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)


data3 <- read.csv("~/PhD/plasma/vac69a/big_plasma_table.csv")
data3$timepoint <- gsub("DoD", "Diagnosis", data3$timepoint)

data3[,3:ncol(data3)] <- log2(data3[,3:ncol(data3)])



data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T))

data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))
data3 <- arrange(data3, timepoint, Volunteer)

data3$Sample_ID <- paste(data3$Volunteer, data3$timepoint)

#censor volutneer 03
#data3 <- subset(data3, data3$Volunteer != "v03")



data4 <- data3
data4[, 1:2] <- NULL


trans_data <- data.frame("Sample_ID"=data4$Sample_ID, apply(data4[,1:ncol(data4)-1], 2, function(x) scale(x, center = TRUE, scale = TRUE)))


plasma_levels <- read.delim("~/PhD/plasma/vac69a/analytes_sorted_by_padj.txt")


t_trans_data <- t(trans_data)
colnames(t_trans_data) <- t_trans_data[1,]

plasma_matrix <- t_trans_data[2:nrow(t_trans_data),]

class(plasma_matrix) <- "numeric"


plasma_matrix <- plasma_matrix[match(plasma_levels$Analyte, rownames(plasma_matrix)),]


plasma_matrix <- ifelse(plasma_matrix>2.5, 2.5, plasma_matrix)
plasma_matrix <- ifelse(plasma_matrix < -2.5, -2.5, plasma_matrix)
plasma_matrix <- plasma_matrix[1:22,]

# top anno #

#define colors for timepoints
time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")
time_col[5] <- "peachpuff1"
Volunteer <- c("v02" = "#FB9A99","v03" = "#E31A1C","v05" = "#A6CEE3", "v06" = "#1F78B4", "v07" = "#F0E442")
Timepoint <- c("Baseline"=time_col[4], "Diagnosis"=time_col[2], "T6"=time_col[1], "C45"=time_col[5])


combo_top_anno <- HeatmapAnnotation(gap = unit(2, "mm"), annotation_name_side = "left",
                                    Volunteer = data3$Volunteer,
                                    Timepoint = data3$timepoint,
                                    col=list(Timepoint = c("Baseline"=time_col[4], "Diagnosis"=time_col[2], "T6"=time_col[1], "C45"=time_col[5]),
                                             Volunteer = c("v02" = "#FB9A99","v03" = "#E31A1C","v05" = "#A6CEE3", "v06" = "#1F78B4", "v07" =  "#F0E442")),
                                    simple_anno_size = unit(2.5, "mm"),
                                    annotation_legend_param = list(
                                      Volunteer = list(title = "Volunteer", at = names(Volunteer), labels_gp = gpar(fontsize=8),title_gp = gpar(fontsize=8), legend_gp = gpar(fill = unname(Volunteer), fontsize=8), title_position = "topleft"),
                                      Timepoint = list(title ="Timepoint",at = names(Timepoint), labels_gp = gpar(fontsize=8), title_gp = gpar(fontsize=8), legend_gp = gpar(fill = unname(Timepoint), fontsize=8), title_position = "topleft")
                                    )
)


# left anno #
number_of_hits <- 9
number_of_maybes <- 3

plasma_matrix <- head(plasma_matrix, n=18)


significant <-  c(rep("<0.05",number_of_hits), rep("<0.1", number_of_maybes), rep(">0.1", nrow(plasma_matrix)-number_of_hits-number_of_maybes))
sig <- c("<0.05"="darkgreen",
         "<0.1" = "lightgreen",
         ">0.1"= "lightgrey")


left_anno <-  rowAnnotation(gap = unit(5, "mm"),
                            #annotation_name_gp = gpar(angle=45),
                            show_annotation_name = FALSE,
                            "significant"=significant,
                            simple_anno_size = unit(2.5, "mm"), # width of the significance bar
                            col=list("significant" = c("<0.05"="darkgreen", "<0.1"="lightgreen", ">0.1" ="lightgrey")),
                            annotation_legend_param = list(significant = list(title ="FDR",
                                                                              at = rev(names(sig)),
                                                                              #title_gp=gpar(angle=45),
                                                                              legend_gp = gpar(fill = unname(sig), fontsize=8),
                                                                              title_gp = gpar(fontsize=10, fontface="bold"),
                                                                              labels_gp = gpar(fontsize=8),
                                                                              title_position = "topleft")
                            )
                            
)




# right anno #

#heatmap construction #


col_fun4 <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#0859C6", "black", "#FFA500"))



combo_map <- Heatmap(matrix = plasma_matrix,
                     cluster_rows = FALSE,
                     name = "Normalised Plasma Concentration",
                     cluster_columns = FALSE,
                     row_names_side = "left",
                     col = col_fun4,
                     column_names_gp = gpar(fontsize=10),
                     row_names_gp = gpar(fontsize=10),
                     #row_split = factor(rep(c("up", "down"), each = 12), levels = c("up", "down")),
                     rect_gp = gpar(col = "white"),
                     #row_title = c("",""),
                     top_annotation = combo_top_anno,
                     #right_annotation = combo_right_anno,
                     left_annotation = left_anno,
                     show_heatmap_legend = TRUE,
                     column_names_rot = 45,
                     heatmap_legend_param = list(col = col_fun4, title = "Z-Score", title_position = "topleft"),
                     #width = unit(16, "cm"),
                     #height = unit(16, "cm")
)



png("~/PhD/cytof/vac69a/final_figures_for_paper/plasma_zscore_heatmap.png", width=7, height=4.5, units = "in", res=400)
draw(combo_map,
     merge_legends = TRUE,
     #padding = unit(c(2, 20, 2, 2), "mm")
)
dev.off()
