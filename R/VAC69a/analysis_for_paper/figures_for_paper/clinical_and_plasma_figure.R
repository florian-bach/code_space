# Panel A clinical presentation of volunteers undergoing vivax CHMI: plateltets, fever lymphocytes, adverse events, parasitaemia ####

# adverse events
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

adverse_pies <- long_data %>%
  mutate(flo_timepoint=gsub(" pm", "", flo_timepoint)) %>%
  mutate(flo_timepoint=gsub(" am", "", flo_timepoint)) %>%
  group_by(Symptom, flo_timepoint)



fig1_theme <- theme(axis.title.x = element_blank(),
                    legend.title = element_text(size = 9), 
                    legend.text = element_text(size = 9),
                    axis.title=element_text(size=10))


(colored_stack <- ggplot(adverse_events,  aes(x=factor(flo_timepoint, levels=timepoint_levels), y=ae_count/6, fill=as.character(Severity)))+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values =  list("1"="#FACA0F", "2"="chocolate1", "3"="red"))+
  #facet_wrap(~Volunteer)+
  ylab("Mean # of AEs\nper Volunteer")+
  xlab("Timepoint")+
  ggtitle("Adverse Events")+
  geom_vline(aes(xintercept = 21.5))+
  scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_minimal()+
  fig1_theme+
  theme(axis.text.x = element_text(hjust=1, angle=45, size=5),
        plot.title = element_text(hjust=0.5, vjust=0, size=10),
        panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Severity",
                            override.aes = list(size = 0.1),
                           keywidth = 0.5,
                           keyheight = 0.5)))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/adverse_events_stacked.png", colored_stack, height=4, width=6)

  # fever
temp_data <- read.csv("~/PhD/clinical_data/vac69a/vac69a_body_temp.csv")

fever <- temp_data %>%
  filter(study_event_oid != "SE_VA69_C28") %>%
  select(trial_number, timepoint, temp) %>%
  #filter(temp>37.5)
  filter(timepoint %in% c("03_chall___01_am", "17_chall___08_pm", "24_chall___12_am","28_chall___14_am",
                          "56_diagnosis_or_c_21", "58_postchall_ep1", "59_postchall_ep2", "62_postchall_t___6"))


lousy_timepoints <- c("03_chall___01_am", "17_chall___08_pm", "24_chall___12_am","28_chall___14_am",
                       "56_diagnosis_or_c_21", "58_postchall_ep1", "59_postchall_ep2", "62_postchall_t___6")

good_timepoints <- c("Baseline", "C8", "C12","C14", "Diagnosis", "T1", "T2", "T6")

timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)
fever$timepoint <- stringr::str_replace_all(fever$timepoint, timepoint_replacement)

fever$volunteer <- gsub("69010", "v", fever$trial_number)


volunteer_colours <- list("v02" = "#FB9A99",
                          "v03" = "#E31A1C",
                          "v05" = "#A6CEE3",
                          "v06" = "#1F78B4",
                          "v07" = "#F0E442",
                          "v09" = "#E69F00")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

fever_curves <- ggplot(fever, aes(x=factor(timepoint, levels=good_timepoints), y=temp, color=volunteer, group=volunteer))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=volunteer), size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  ggtitle("Fever Curves")+
  geom_hline(yintercept = 37.5, linetype="dashed", color="black")+
  xlab("Timepoint")+
  ylab(expression(paste("Temperature (",degree,"C)",sep="")))+
  scale_y_continuous(breaks = seq(36.5, 40, by=0.5) )+
  theme_minimal()+
  fig1_theme+
  theme(plot.title = element_text(hjust=0.5, vjust = 0, size=10),
        axis.text.x = element_text(hjust=1, angle=45, size=8),
        legend.position = "none", 
        axis.title.x = element_blank())

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_curves.png", fever_curves, width=2, height=2)


# PLatelets/Lymphocytes


haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$trial_number <- gsub("69010", "v", haem_data$trial_number)

long_haem_data <- haem_data %>%
  select(trial_number, timepoint, platelets, lymphocytes) %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_EP", "_T6")) %>%
  gather(Cell, Frequency, c(platelets,lymphocytes))




long_haem_data$Cell <- paste(
  toupper(substr(long_haem_data$Cell, 1,1)),
  substr(long_haem_data$Cell, 2,nchar(long_haem_data$Cell)),
  sep="")


bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_EP", "_T6")
great_timepoints <- c("Baseline", "C8", "C14", "T1", "T6")

time_dic <- setNames(great_timepoints, bad_timepoints)



long_haem_data$timepoint <- stringr::str_replace_all(long_haem_data$timepoint, time_dic)


thrombos_lymphs <- ggplot(long_haem_data, aes(x=factor(timepoint, levels=unique(gtools::mixedsort(long_haem_data$timepoint))), y=Frequency*1000, color=trial_number, group=trial_number))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=trial_number), size=0.9)+
  geom_point(fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(~Cell, scales="free")+
  xlab("Timepoint")+
  ylab(expression(Cells~"/"~mu*L~blood))+
  guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
  fig1_theme+
  scale_y_continuous(label=scales::comma)+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=10))


ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/thrombos_lymphs.png", thrombos_lymphs, width=6, height=2.2)



# lgd1 <- get_legend(thrombos)
# 
thrombos <- thrombos+theme(legend.position = "none")
# 
# 
# lgd2 <- get_legend(thrombos)
# 
# thrombos <- thrombos+theme(legend.position = "none")
# adverse events color pies

haem_pies <- haem_data %>%
  select(trial_number, timepoint, grep("_ae", colnames(haem_data), value = T)) %>%
  gather(Symptom, Severity, grep("_ae", colnames(haem_data), value = T)) %>%
  filter(Symptom=="lymphocytes_ae") %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_EP", "_T6"))


bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_EP", "_T6")
great_timepoints <- c("Baseline", "C8", "C14", "T1", "T6")

time_dic <- setNames(great_timepoints, bad_timepoints)



haem_pies$timepoint <- stringr::str_replace_all(haem_pies$timepoint, time_dic)

haem_pies$Severity <- as.character(haem_pies$Severity)
haem_pies$Symptom <- substr(haem_pies$Symptom, 1,nchar(haem_pies$Symptom)-3)

haem_pies$Symptom <- paste(
  toupper(substr(haem_pies$Symptom, 1,1)),
  substr(haem_pies$Symptom, 2,nchar(haem_pies$Symptom)),
  sep="")

haem_pies$Symptom <- gsub("Wbc", "WBC", haem_pies$Symptom)
  
  lymphocyte_ae_pie <- ggplot(haem_pies,  aes(x="", fill=factor(Severity)))+
    geom_bar(stat="count", position ="stack")+
    scale_fill_manual(values =  list("0"="lightgrey", "1"="#FACA0F", "2"="chocolate1", "3"="red"))+
    facet_wrap(~timepoint, ncol=2)+
    coord_polar(theta = "y")+
    ggtitle("Adverse Events Lymphopenia\n")+
    #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
    theme_void()+
    fig1_theme+
    theme(plot.title = element_text(hjust=0.5, vjust=0, size=10),
          panel.grid.minor = element_blank(),
          axis.title = element_blank())+
    guides(fill=guide_legend(title="Severity",
                             override.aes = list(size = 0.1),
                             keywidth = 0.5,
                             keyheight = 0.5))
  
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lymphocyte_ae_pie.png", lymphocyte_ae_pie, height=2, width=3)
  
  
  adverse_pies <- long_data %>%
    mutate(flo_timepoint=gsub(" pm", "", flo_timepoint)) %>%
    mutate(flo_timepoint=gsub(" am", "", flo_timepoint)) %>%
    filter(flo_timepoint %in% c("Diagnosis", "T1 ","T2 ", "T6 ")) %>%
    filter(Symptom %notin% c("vomiting", "nausea","pyrexia", "diarrhoea", "back_pain", "arthralgia")) %>%
    select(Symptom, flo_timepoint, Severity)
  
  
  
  
  
  ae_pies <- ggplot(adverse_pies,  aes(x="", fill=factor(Severity)))+
    geom_bar(stat="count", position ="stack")+
    scale_fill_manual(values =  list("0"="lightgrey", "1"="#FACA0F", "2"="chocolate1", "3"="red"))+
    facet_grid(flo_timepoint~Symptom, switch="y")+
    coord_polar(theta = "y")+
    ggtitle("Adverse Events Symptoms\n")+
    #scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
    theme_void()+
    fig1_theme+
    theme(plot.title = element_text(hjust=0.5, vjust=0, size=10),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          strip.text.y.left = element_text(hjust=0.5, size=7, angle = 0))+
    guides(fill=guide_legend(title="Severity",
                             override.aes = list(size = 0.1),
                             keywidth = 0.5,
                             keyheight = 0.5))
  
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/ae_pies.png", ae_pies, width=7, height=4)



    # parasitaemias ####

data <- read.csv("~/PhD/clinical_data/vac69a/parasitaemia/better_vac69a_parasitaemia.csv", header=T)
parasitaemias <- gather(data, Timepoint, Genomes, colnames(data)[2:ncol(data)])
parasitaemias$Genomes <- as.numeric(parasitaemias$Genomes)
parasitaemias$Volunteer <- gsub("MVT-069010", "v",parasitaemias$Volunteer )
# get rid of garbage timepoints that mess up graph

parasitaemias$Timepoint <- ifelse(grepl(".5", parasitaemias$Timepoint, fixed = T), paste(parasitaemias$Timepoint, "pm"),  paste(parasitaemias$Timepoint, "am"))
parasitaemias$Timepoint <- gsub(".5", "", parasitaemias$Timepoint, fixed=T)
parasitaemias$Timepoint <- gsub("D", "C", parasitaemias$Timepoint, fixed=T)

dods <- data.frame(vol=c('v02', 'v03', 'v05', 'v06', 'v07', 'v09'),
                   paras=c(4907, 8054, 16733, 7464, 21870, 15051),
                   dod=c('D15.5', 'D12.5', 'D15.5', 'D15.5', 'D16', 'D17'))
dods$inter <- match(dods$dod, unique(parasitaemias$Timepoint))





parasitaemia_curves <- ggplot(data=parasitaemias[!is.na(parasitaemias$Genomes),], aes(x=factor(Timepoint, levels=unique(gtools::mixedsort(parasitaemias$Timepoint))), y=Genomes+1, group=factor(Volunteer)))+
  geom_line(aes(color=factor(Volunteer)), size=0.9)+
  geom_point(fill="white", stroke=1, size=0.9, shape=21, aes(color=factor(Volunteer)))+
  scale_color_manual(values=volunteer_palette)+
  theme_minimal()+
  ggtitle("Parasitaemia")+
  xlab("Day of Infection")+
  ylab("Genome Copies / mL")+
  scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000))+
  guides(color=guide_legend(override.aes = list(size = 0.1),
                           keywidth = 0.5,
                           keyheight = 0.5))+

  fig1_theme+
  theme(legend.title = element_blank(),
        plot.title = element_text(size=10, hjust=0.5, vjust=0),
        axis.text.x = element_text(hjust=1, angle=45, size=5))

#ggsave ("/Users/s1249052/PhD/oxford/vac69/parasitaemias_vac69.png", height = 8, width=10)
ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vac69a_parasitaemia.png", parasitaemia_curves, width=6, height=5)

fig1 <- cowplot::plot_grid(fever_curves, parasitaemia_curves, 
                           thrombos_lymphs, colored_stack, align = "h", axis = "tblr", rel_widths=c(1.6, 2))

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fig1_no_heatmap.png", fig1, height = 4, width=8)

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

#supp Panel B Significant Plasma Analytes ####

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


ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/supp_glm_sig_analytes_fdr_10-e-1.png", sig_glm_plot, width = 8, height=5)

