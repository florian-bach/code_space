# Panel ABC clinical presentation of volunteers undergoing vivax CHMI: plateltets, fever and adverse events ####


data <- read.csv("~/PhD/clinical_data/vac69a/symptoms_vac69a.csv", header = T, stringsAsFactors = F)
data$flo_timepoint <- gsub("am ", "am", data$flo_timepoint)
data$flo_timepoint <- gsub("pm ", "pm", data$flo_timepoint)
data$flo_timepoint <- gsub("DoD ", "Diagnosis", data$flo_timepoint)

long_data <- tidyr::gather(data, Symptom, Severity, colnames(data)[c(12, 14:ncol(data))])

long_data <- mutate(long_data, Volunteer=gsub("69010", "V", long_data$trial_number))

timepoint_levels <- list("Timepoints"=unique(long_data$flo_timepoint[gtools::mixedorder(long_data$flo_timepoint)]))

timepoint_levels$Timepoints[[length(timepoint_levels$Timepoints)+1]] <- timepoint_levels$Timepoints[32]

timepoint_levels <- timepoint_levels$Timepoints[-32]

myLoc <- (which(levels(long_data$flo_timepoint) == "DoD ") +
            which(levels(long_data$flo_timepoint) == "C15 pm")) / 
  2

symptom_heatmap <- ggplot(long_data, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=Symptom))+
  geom_tile(aes(fill=factor(Severity), width=0.93, height=0.93), color=ifelse(grepl("DoD", long_data$flo_timepoint), "black", "lightgrey"))+
  scale_fill_manual(values =  list("lightgrey", "yellow", "orange", "red"))+
  facet_wrap(~Volunteer, scales="free")+
  theme_minimal()+
  guides(fill=guide_legend(title="Severity"))+
  theme(axis.text.x = element_text(hjust=1, angle=45),
        axis.title = element_blank(),
        legend.title = element_blank())

ggsave("/home/flobuntu/PhD/clinical_data/vac69a/figures/symptom_heatmap.png", symptom_heatmap, height=8, width=14)

# make a figure for number of AEs per timepoint

library(dplyr)
library(tidyr)

long_data$flo_timepoint <- factor(long_data$flo_timepoint)


adverse_events <- long_data %>%
  filter(Severity > 0) %>%
  group_by(Volunteer, flo_timepoint, Severity) %>%
  summarise(ae_count = n())


colored_stack <- ggplot(adverse_events,  aes(x=factor(flo_timepoint, levels=timepoint_levels), y=ae_count/6, fill=factor(Severity, levels=paste(rev(1:3)))))+
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values =  list("1"="yellow", "2"="orange", "3"="red"))+
  #facet_wrap(~Volunteer)+
  ylab("Average # of AEs\nper Volunteer")+
  xlab("Timepoint")+
  ggtitle("Adverse Events")+
  geom_vline(aes(xintercept = 21.5))+
  guides(fill=guide_legend(title="Severity"))+
  scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, angle=45, size=5),
        plot.title = element_text(hjust=0.5),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/adverse_events_stacked.png", colored_stack, height=4, width=6)



fever <- subset(long_data, long_data$pyrexia_temp>37)

volunteer_colours <- list("V02" = "#FB9A99",
                          "V03" = "#E31A1C",
                          "V05" = "#A6CEE3",
                          "V06" = "#1F78B4",
                          "V07" = "#B2DF8A",
                          "V09" = "#33A02C")


volunteer_palette <- unlist(unname(volunteer_colours))
names(volunteer_palette) <- names(volunteer_colours)

fever_curves <- ggplot(fever, aes(x=factor(flo_timepoint, levels=timepoint_levels), y=pyrexia_temp, color=Volunteer, group=Volunteer))+
  scale_fill_manual(values=volunteer_palette)+
  scale_color_manual(values=volunteer_palette)+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  ggtitle("Fever")+
  xlab("Timepoint")+
  ylab(expression(paste("Temperature ",degree,"C",sep="")))+
  scale_y_continuous(limits=c(37.5, 40), breaks = seq(37.5, 40, by=0.5) )+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8),
        legend.position = "none", 
        axis.title.x = element_blank())

ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_curves.png", fever_curves)


# PLatelets/Lymphocytes


haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$trial_number <- gsub("69010", "V", haem_data$trial_number)

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
  ylab(expression(Cells~per~mu*L))+
  guides(color=guide_legend(title="Volunteer"))+
  ggtitle("Platelets")+
  scale_y_continuous(label=scales::label_scientific())+
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




# Panel C Significant Plasma Analytes ####

library(tidyr)
library(dplyr)
library(ggplot2)

data3 <- read.csv("~/PhD/plasma/vac69a/big_plasma_table.csv")


data3[,3:ncol(data3)] <- log2(data3[,3:ncol(data3)])


data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T))

data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "DoD", "T6", "C45"))


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

sig_glm_plot <- ggplot(sig_glm_data, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45")), y=Concentration, color=Volunteer, group=Volunteer))+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  facet_wrap(~ AnalyteF, scales = "free", ncol=4)+
  scale_y_log10()+
  ylab("Concentation (log2 pg/ mL)")+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1, angle=45),
        strip.background = element_rect(fill = "white", color = "white"))


ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/glm_sig_analytes_fdr_10-e-1.png", sig_glm_plot, width = 8, height=5)

clinical_graphs <- cowplot::plot_grid(fever_curves, thrombos, colored_stack, ncol = 3, rel_widths = c(0.8,1,1.7))
#ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/aes_and_fever_curves.png", clinical_graphs, height=3, width=8)

fig1 <- plot_grid(clinical_graphs, sig_glm_plot, nrow=2, ncol=1, label="auto", rel_heights = c(1,3))

ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/fig1_chunk.png", fig1, width = 8, height=6)
