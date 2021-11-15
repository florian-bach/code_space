# longitudinal haem ####
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)


vac69_colours <-list( "v02" = "#FB9A99",
                      "v03" = "#E31A1C",
                      "v05" = "#A6CEE3",
                      "v06" = "#1F78B4",
                      "v07" = "#F0E442",
                      "v09" = "#E69F00",
                      "v11" = "#D33F6A",
                      "v21" = "#E99A2C",
                      "v34" = "red", 
                      "v35" = "#6C3461",
                      "v37" = "#5E9B8A",
                      "v39" = "#FF08E8",
                      "v40" = "#61E160",
                      "v42" = "#6488EA",
                      "v44" = "darkblue")

vac69_colours <- unlist(unname(vac69_colours))
names(vac69_colours) <- names(vac69_colours)

# vac69a ####

haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")

haem_data$volunteer <- gsub("69010", "v", haem_data$trial_number)

supp_haem_data <- haem_data %>%
  select(volunteer, timepoint, wbc, neutrophils,  monocytes, eosinophils, lymphocytes, platelets) %>%
  #select(trial_number, timepoint, lymphocytes) %>%
  filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")) %>%
  gather(Cell, Frequency, c(wbc, neutrophils,  monocytes, eosinophils, lymphocytes, platelets))


supp_haem_data$Cell <- paste(
  toupper(substr(supp_haem_data$Cell, 1,1)),
  substr(supp_haem_data$Cell, 2,nchar(supp_haem_data$Cell)),
  sep="")

supp_haem_data$Cell <- gsub("Wbc", "Total White Cells", supp_haem_data$Cell)

supp_haem_data$Cell <- factor(supp_haem_data$Cell, levels=c("Total White Cells", "Lymphocytes", "Monocytes", "Neutrophils", "Eosinophils", "Platelets"))


bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")
great_timepoints <- c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")

time_dic <- setNames(great_timepoints, bad_timepoints)

supp_haem_data$timepoint <- stringr::str_replace_all(supp_haem_data$timepoint, time_dic)

vac69a_haem <- supp_haem_data

vac69a_haem$n_infection <- "First"

#vac69b ####

data <- read.csv("~/PhD/clinical_data/vac69b/haem.csv", header=T)


data_no_ae <- select(data, -c(colnames(data)[grep("_ae",colnames(data) ,fixed=T)]))


data_no_ae$timepoint <- data_no_ae$timepoint %>%
  gsub("B_C1", "Baseline", .) %>%
  gsub("B_", "", .) %>%
  gsub("C814", "C14", .) %>%
  gsub("C28DOD", "Diagnosis", ., fixed = TRUE) %>%
  gsub("C67", "C7", .)

data_no_ae$timepoint <- gsub("Baseline527", "Baseline", data_no_ae$timepoint)

long_data <- gather(data_no_ae, Cell, Frequency, colnames(data_no_ae)[7:13])

vac69b_data <- filter(long_data, long_data$flo_timepoint %in% c("Baseline", "C7", "C14", "DoD", "T6", "C96"))

vac69b_data$volunteer <- paste("v", substr(vac69b_data$trial_number, 6, 8), sep="")

vac69b_data$n_infection <- ifelse(vac69b_data$volunteer %in% c("v11", "v21"), "First", "Second")


#vac69c ####
haem <- read.csv("~/PhD/clinical_data/vac69c/haem.csv", header=TRUE, stringsAsFactors = FALSE)

haem$volunteer <- paste("v", substr(haem$trial_number, 6, 8), sep="")

lousy_haem_timepoints <- unique(haem$timepoint)
# [1] "C_RESCREEN" "C_C2"       "C_C17"      "C_C814"     "C_C28DIAG"  "C_T1"       "C_T3"      
# [8] "C_T6"       C_C56"  "C_EV"    "SCR"        "C_C96"

good_haem_timepoints <- c("ReScreening", "Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "EV", "Screening", "C96")

haem_timepoint_replacement <- setNames(good_haem_timepoints, lousy_haem_timepoints)
haem$timepoint <- stringr::str_replace_all(haem$timepoint, haem_timepoint_replacement)

haem$timepoint <- gsub("Baseline8DIAG", "Diagnosis", haem$timepoint)

haem$n_infection <- ifelse(haem$volunteer == "v07", "Third",
                           ifelse(haem$volunteer %in% c("v11", "v21"), "Second", "First"))

vac69c_haem <- haem %>%
  mutate("haemoglobin"=haemaglobin)%>%
  select(volunteer, timepoint, n_infection, date_performed, haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils) %>%
  pivot_longer(cols = c(haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils), names_to = "Cell", values_to = "Frequency")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))

# putting it all together ####
vac69a_haem <- vac69a_haem %>%
  filter(timepoint %in% c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6"),
         Cell %in% c("Lymphocytes", "Platelets"))

vac69b_haem <- vac69b_data %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6"),
         Cell %in% c("lymphocytes", "platelets"))

vac69c_haem <- vac69c_haem %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6"),
         Cell %in% c("lymphocytes", "platelets"))


big_table <- data.frame("volunteer"=c(vac69a_haem$volunteer, vac69b_haem$volunteer, vac69c_haem$volunteer),
                        "timepoint"=c(vac69a_haem$timepoint, vac69b_haem$timepoint, vac69c_haem$timepoint),
                        "n_infection"=c(vac69a_haem$n_infection, vac69b_haem$n_infection, vac69c_haem$n_infection),
                        "Cell"=c(as.character(vac69a_haem$Cell), vac69b_haem$Cell, vac69c_haem$Cell),
                        "Frequency"=c(vac69a_haem$Frequency, vac69b_haem$Frequency, vac69c_haem$Frequency))

big_table$timepoint <- gsub(" am", "", big_table$timepoint)
big_table$timepoint <- gsub(" pm", "", big_table$timepoint)

big_table$Cell <- gsub("platelets", "Platelets", big_table$Cell)
big_table$Cell <- gsub("lymphocytes", "Lymphocytes", big_table$Cell)


big_lymph <- filter(big_table, Cell=="Lymphocytes")
big_thromb <- filter(big_table, Cell=="Platelets")

lymph_box_plot <- ggplot(big_lymph, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency))+
  scale_color_manual(values=vac69_colours)+
  geom_boxplot(aes(fill=n_infection))+
  scale_fill_manual(values=n_infection_palette)+
  #geom_line(aes(color=volunteer, group=volunteer), size=0.9)+
  geom_point(aes(color=volunteer, group=n_infection), fill="white",  position = position_dodge(width = 0.75), stroke=1, shape=21, size=0.9)+
  ggtitle("Lymphocytes VAC69abc")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Frequency")+
  guides(color=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abc_lymph_box.png", lymph_box_plot, height=4, width=6, bg="white")



thromb_box_plot <- ggplot(big_thromb, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency))+
  scale_color_manual(values=vac69_colours)+
  geom_boxplot(aes(fill=n_infection), colour="black")+
  scale_fill_manual(values=n_infection_palette)+
  #geom_line(aes(color=volunteer, group=volunteer), size=0.9)+
  geom_point(aes(color=volunteer, group=n_infection), fill="white",  position = position_dodge(width = 0.75), stroke=1, shape=21, size=0.9)+
  ggtitle("Platelets VAC69abc")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Frequency")+
  guides(color=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abc_thromb_box.png", thromb_box_plot, height=4, width=6, bg="white")
  

big_combo <- rbind(big_lymph, big_thromb)

big_combo <- filter(big_combo, volunteer %in% c("v02", "v05", "v07", "v11", "v21"))


abc_line_plot <- ggplot(big_combo, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(aes(color=n_infection, group=n_infection), size=0.9)+
  geom_point(aes(color=n_infection), fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(Cell~volunteer, scales="free", ncol=2)+
  xlab("Timepoint")+
  ylab("Frequency")+
  guides(color=guide_legend(ncol=1))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abc_line_plot.png", abc_line_plot, height=10, width=5, bg="white")

# longitudinal biochem ####

# vac69a biochem ####
vac69a_alt <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")


colnames(vac69a_alt)[1] <- "volunteer"

vac69a_alt$volunteer <- paste("v", substr(vac69a_alt$volunteer, 6, 7), sep='')



lousy_timepoints <- unique(as.character(vac69a_alt$timepoint))

good_timepoints <- c("C7", "C14", "Screening","Baseline", "Diagnosis", "C28", "T1", "T6", "C90", "EV")

biochem_timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)
vac69a_alt$timepoint <- stringr::str_replace_all(vac69a_alt$timepoint, biochem_timepoint_replacement)

vac69a_alt <- select(vac69a_alt,  volunteer, timepoint, alt)

vac69a_alt <- vac69a_alt %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6", "C90"), volunteer %in% c("v02", "v05", "v07"))

vac69a_alt$n_infection <- "First"


#vac69b biochem ####



vac69b_alt <- read.csv("~/PhD/clinical_data/vac69b/biochem.csv", header=T)

vac69b_alt <- select(vac69b_alt, -c(colnames(vac69b_alt)[grep("_ae",colnames(vac69b_alt) ,fixed=T)]))

vac69b_alt <- data.frame("volunteer"=paste("v", substr(vac69b_alt$trial_number, 6,8), sep=""),
                          "timepoint"=gsub("DoD", "Diagnosis", vac69b_alt$flo_timepoint),
                          "alt"=vac69b_alt$alt,
                          "n_infection"=ifelse(vac69b_alt$trial_number %in% c(6901021, 6901011), "First", "Second"))


vac69b_alt <- filter(vac69b_alt, vac69b_alt$timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6", "C90"))


# vac69c biochem ####



biochem <- read.csv("~/PhD/clinical_data/vac69c/biochem.csv", header=TRUE, stringsAsFactors = FALSE)

biochem$volunteer <- paste("v", substr(biochem$trial_number, 6, 8), sep="")

lousy_biochem_timepoints <- unique(biochem$timepoint)
# [1] "C_RESCREEN" "C_C2"       "C_C17"      "C_C814"     "C_C28DIAG"  "C_T1"       "C_T3"      
# [8] "C_T6"       "C_EV"       "C_C56"      "SCR"        "C_C96"

good_biochem_timepoints <- c("ReScreening", "Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "EV", "C56", "Screening", "C96")

biochem_timepoint_replacement <- setNames(good_biochem_timepoints, lousy_biochem_timepoints)
biochem$timepoint <- stringr::str_replace_all(biochem$timepoint, biochem_timepoint_replacement)

biochem$timepoint <- gsub("Baseline8DIAG", "Diagnosis", biochem$timepoint)

biochem$n_infection <- ifelse(biochem$volunteer == "v07", "Third",
                              ifelse(biochem$volunteer %in% c("v11", "v21"), "Second", "First"))

quant_biochem <- biochem %>%
  select(volunteer, timepoint, n_infection, date_performed, sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin) %>%
  pivot_longer(cols = c(sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin), names_to = "Analyte", values_to = "Concentration")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


vac69c_alt <- quant_biochem %>%
  filter(Analyte=="alt", n_infection%in%c("Second", "Third")) %>%
  mutate("alt"=Concentration)%>%
  dplyr::select(volunteer, timepoint, n_infection, alt)


combo_alt <- rbind(vac69a_alt, vac69b_alt, vac69c_alt)

combo_alt <- filter(combo_alt, timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))

abc_alt_line_plot <- ggplot(combo_alt, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=alt))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(aes(color=n_infection, group=n_infection), size=0.9)+
  geom_point(aes(color=n_infection), fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(~volunteer, ncol=5)+
  xlab("Timepoint")+
  ylab("ALT")+
  guides(color=guide_legend(ncol=1))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abc_alt_line_plot.png", abc_alt_line_plot, height=3, width=8, bg="white")
