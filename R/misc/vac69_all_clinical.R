# preamble ####
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
                      "v44" = "darkblue",
                      "v47" = "black",
                      "v49" = "black",
                      "v50" = "black",
                      "v51" = "black")

vac69_colours <- unlist(unname(vac69_colours))
names(vac69_colours) <- names(vac69_colours)



n_infection_palette <- c(rgb(5,50,80, maxColorValue = 255),
                         rgb(250, 100, 0, maxColorValue = 255),
                         rgb(40,210,250, maxColorValue = 255))



# longitudinal haem ####

# vac69a ##

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

#vac69b ##

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


#vac69c ##
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

#vac69c ##


vac69d_haem <- read.csv("~/PhD/clinical_data/vac69d/vac69d_haem.csv", header=TRUE, stringsAsFactors = FALSE)
vac69d_haem$volunteer <- paste("v", substr(vac69d_haem$trial_number, 6, 8), sep="")


lousy_haem_timepoints <- unique(vac69d_haem$timepoint)
# [1] "D_RESCREEN" "DC2"        "D_C17"      "D_C21DIAG"  "D_C814"     "D_T1"       "D_T3"       "D_T6"      
# [9] "D_C56"      "D_CHALLEV

good_haem_timepoints <- c("ReScreening", "Baseline", "C7", "Diagnosis", "C14", "T1", "T3", "T6", "C56", "EV")


haem_timepoint_replacement <- setNames(good_haem_timepoints, lousy_haem_timepoints)
vac69d_haem$timepoint <- stringr::str_replace_all(vac69d_haem$timepoint, haem_timepoint_replacement)

vac69d_haem$n_infection <- ifelse(vac69d_haem$volunteer %in% c("v11", "v21"), "Third",
                                  ifelse(vac69d_haem$volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First")
                                  )

slim_vac69d <- vac69d_haem %>%
  mutate("haemoglobin"=haemaglobin)%>%
  select(volunteer, timepoint, n_infection, date_performed, haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils) %>%
  pivot_longer(cols = c(haemoglobin, wbc, platelets, haematocrit, red_cell_count, neutrophils, lymphocytes, monocytes, eosinophils, basophils), names_to = "Cell", values_to = "Frequency")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


# putting it all together ##
vac69a_haem <- vac69a_haem %>%
  filter(timepoint %in% c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6"),
         Cell %in% c("Lymphocytes", "Platelets"))

vac69b_haem <- vac69b_data %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6"),
         Cell %in% c("lymphocytes", "platelets"))

vac69c_haem <- vac69c_haem %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6"),
         Cell %in% c("lymphocytes", "platelets"))

vac69d_haem <- slim_vac69d %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6"),
         Cell %in% c("lymphocytes", "platelets"))

big_table <- data.frame("volunteer"=c(vac69a_haem$volunteer, vac69b_haem$volunteer, vac69c_haem$volunteer, vac69d_haem$volunteer),
                        "timepoint"=c(vac69a_haem$timepoint, vac69b_haem$timepoint, vac69c_haem$timepoint, vac69d_haem$timepoint),
                        "n_infection"=c(vac69a_haem$n_infection, vac69b_haem$n_infection, vac69c_haem$n_infection, vac69d_haem$n_infection),
                        "Cell"=c(as.character(vac69a_haem$Cell), vac69b_haem$Cell, vac69c_haem$Cell, vac69d_haem$Cell),
                        "Frequency"=c(vac69a_haem$Frequency, vac69b_haem$Frequency, vac69c_haem$Frequency, vac69d_haem$Frequency))

big_table$timepoint <- gsub(" am", "", big_table$timepoint)
big_table$timepoint <- gsub(" pm", "", big_table$timepoint)

big_table$Cell <- gsub("platelets", "Platelets", big_table$Cell)
big_table$Cell <- gsub("lymphocytes", "Lymphocytes", big_table$Cell)


big_lymph <- filter(big_table, Cell=="Lymphocytes")
big_thromb <- filter(big_table, Cell=="Platelets")

big_lymph <- filter(big_lymph, n_infection!="Third")

lymph_box_plot <- ggplot(big_lymph, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency*1e6))+
  #scale_color_manual(values=vac69_colours)+
  geom_boxplot(aes(fill=n_infection))+
  scale_fill_manual(values=n_infection_palette)+
  scale_y_continuous(labels=scales::scientific_format())+
  #geom_line(aes(color=volunteer, group=volunteer), size=0.9)+
  geom_point(aes(color=volunteer, group=n_infection), fill="white",  position = position_dodge(width = 0.75), stroke=1, shape=21, size=0.9)+
  ggtitle("Lymphocytes VAC69abcd")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab("Cells / mL")+
  guides(color=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.title =element_blank(),
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_lymph_box.png", lymph_box_plot, height=5, width=6, bg="white")


big_thromb <- filter(big_thromb, n_infection!="Third")

thromb_box_plot <- ggplot(big_thromb, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency*1000))+
  #scale_y_log10()+#scale_color_manual(values=vac69_colours)+
  geom_boxplot(aes(fill=n_infection), colour="black")+
  scale_fill_manual(values=n_infection_palette)+
  #geom_line(aes(color=volunteer, group=volunteer), size=0.9)+
  geom_point(aes(color=volunteer, group=n_infection), fill="white",  position = position_dodge(width = 0.75), stroke=1, shape=21, size=0.9)+
  ggtitle("Thrombocytes VAC69abcd")+
  theme_minimal()+
  xlab("Timepoint")+
  ylab(expression(Cells~"/"~mu*L))+
  guides(color=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        legend.title =element_blank(),
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_thromb_box.png", thromb_box_plot, height=5, width=6, bg="white")
  

big_combo <- rbind(big_lymph, big_thromb)

rechallenged_volunteers <- unique(subset(big_combo, n_infection=="Second")$volunteer)

big_combo <- filter(big_combo, n_infection != "Third")


abc_line_plot <- ggplot(big_combo, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=Frequency))+
  scale_color_manual(values=n_infection_palette)+
  geom_line(aes(color=n_infection, group=n_infection), size=0.9)+
  geom_point(aes(color=n_infection), fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  facet_wrap(volunteer~Cell, scales="free", ncol=4)+
  xlab("Timepoint")+
  ylab("Frequency")+
  guides(color=guide_legend(ncol=1))+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(hjust=1, angle=45, size=8), 
        axis.title.x = element_blank(),
        #legend.position = "none",
        strip.text = element_text(size=10))

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_line_plot.png", abc_line_plot, height=13, width=10, bg="white")


# longitudinal biochem ####

# vac69a biochem ##
vac69a_alt <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")


colnames(vac69a_alt)[1] <- "volunteer"

vac69a_alt$volunteer <- paste("v", substr(vac69a_alt$volunteer, 6, 7), sep='')



lousy_timepoints <- unique(as.character(vac69a_alt$timepoint))

good_timepoints <- c("C7", "C14", "Screening","Baseline", "Diagnosis", "C28", "T1", "T6", "C90", "EV")

biochem_timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)
vac69a_alt$timepoint <- stringr::str_replace_all(vac69a_alt$timepoint, biochem_timepoint_replacement)

vac69a_alt <- select(vac69a_alt,  volunteer, timepoint, alt)

vac69a_alt <- vac69a_alt %>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6", "C90"))
# 
# 
# vac69a_alt <- vac69a_alt %>%
#   filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6", "C90"), volunteer %in% c("v02", "v05", "v07"))

vac69a_alt$n_infection <- "First"


#vac69b biochem ##



vac69b_alt <- read.csv("~/PhD/clinical_data/vac69b/biochem.csv", header=T)

vac69b_alt <- select(vac69b_alt, -c(colnames(vac69b_alt)[grep("_ae",colnames(vac69b_alt) ,fixed=T)]))

vac69b_alt <- data.frame("volunteer"=paste("v", substr(vac69b_alt$trial_number, 6,8), sep=""),
                          "timepoint"=gsub("DoD", "Diagnosis", vac69b_alt$flo_timepoint),
                          "alt"=vac69b_alt$alt,
                          "n_infection"=ifelse(vac69b_alt$trial_number %in% c(6901021, 6901011), "First", "Second"))


vac69b_alt <- filter(vac69b_alt, vac69b_alt$timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T6", "C90"))


# vac69c biochem ##



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
  filter(Analyte=="alt") %>%
  mutate("alt"=Concentration)%>%
  dplyr::select(volunteer, timepoint, n_infection, alt)

# vac69d biochem ##

vac69d_biochem <- read.csv("~/PhD/clinical_data/vac69d/vac69d_biochem.csv", header=TRUE, stringsAsFactors = FALSE)

#get rid of non vac69d data, extra volunteers
vac69d_biochem <- subset(vac69d_biochem, substr(vac69d_biochem$timepoint, 1, 1)=="D")
vac69d_biochem <- subset(vac69d_biochem, trial_number!=6901054)

vac69d_biochem$volunteer <- paste("v", substr(vac69d_biochem$trial_number, 6, 8), sep="")


lousy_haem_timepoints <- unique(vac69d_biochem$timepoint)
# [1] "D_RESCREEN" "DC2"        "D_C17"      "D_C814"     "D_C21DIAG"  "D_T1"      
# [7] "D_T3"       "D_T6"       "D_C56"      "D_CHALLEV" 

good_haem_timepoints <- c("ReScreening", "Baseline", "C7", "C14",  "Diagnosis", "T1", "T3", "T6", "C56", "EV")

biochem_timepoint_replacement <- setNames(good_haem_timepoints, lousy_haem_timepoints)


vac69d_biochem$timepoint <- stringr::str_replace_all(vac69d_biochem$timepoint, biochem_timepoint_replacement)

vac69d_biochem$n_infection <- ifelse(vac69d_biochem$volunteer %in% c("v11", "v21"), "Third",
                                  ifelse(vac69d_biochem$volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First")
)

vac69d_quant_biochem <- vac69d_biochem %>%
  select(volunteer, timepoint, n_infection, date_performed, sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin) %>%
  pivot_longer(cols = c(sodium, potassium, urea, creatinine, bilirubin, alt, alkphos, albumin), names_to = "Analyte", values_to = "Concentration")%>%
  filter(timepoint %in% c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96"))


vac69d_alt <- vac69d_quant_biochem %>%
  filter(Analyte=="alt") %>%
  mutate("alt"=Concentration)%>%
  dplyr::select(volunteer, timepoint, n_infection, alt)


# putting it all together

combo_alt <- rbind(vac69a_alt, vac69b_alt, vac69c_alt, vac69d_alt)

combo_alt <- filter(combo_alt, timepoint %in% c("Baseline", "C7", "Diagnosis", "T1", "T3", "T6"))

combo_alt$sample_id <- paste(combo_alt$volunteer, combo_alt$n_infection, sep="_")


model <- lm(alt~timepoint+n_infection+volunteer, data=combo_alt)

mixed_model <- lme4::lmer(alt~timepoint+n_infection+(1|volunteer), data=combo_alt)

drop1(mixed_model,test="Chisq")

t6_first_vs_second_contrast <- t(matrix(c(0,0,0,1,1,rep(0,12))))
t6__contrast <- t(matrix(c(0,0,0,1,0, rep(0,12))))

t6_first_vs_second_contrast_mixed <- t(matrix(c(0,0,0,1,1,0)))


model_tests <- multcomp::glht(model, t6_first_vs_second_contrast)

mixed_model_tests <- multcomp::glht(mixed_model, t6_first_vs_second_contrast_mixed)

summary(model_tests) # first vs. second at T6; p=0.116
summary(mixed_model_tests) # first vs. second at T6; p=0.166
#



summary(model)

combo_alt <- filter(combo_alt, n_infection != "Third")

(abc_alt_line_plot <- ggplot(combo_alt, aes(x=factor(timepoint, levels=c("Baseline", "C7", "C14", "Diagnosis", "T1", "T3", "T6", "C56", "C96")), y=alt))+
  #geom_hline(yintercept = 35, color="orange", linetype="dotted")+
  geom_hline(yintercept = 45, color="orange", linetype="dashed")+
  #geom_hline(yintercept = 175, color="red", linetype="dotted")+
  geom_hline(yintercept = 225, color="red", linetype="dashed")+
  scale_fill_manual(values=n_infection_palette)+
  #geom_line(aes(color=volunteer, group=sample_id), size=0.9)+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(color=volunteer, group=n_infection),  position = position_dodge(width = 0.75), fill="white", stroke=1, shape=21, size=0.9)+
  theme_minimal()+
  scale_y_log10()+
  #facet_wrap(~n_infection, ncol=5)+
  xlab("Timepoint")+
  ylab("log(ALT)")+
  ggtitle(expression(paste('abnormal liver enzymes in first ',italic("P. vivax"),' infections', sep='')))+
  guides(color=guide_legend(ncol=2))+
  theme(plot.title = element_text(hjust=0.5, size=18),
        axis.text.x = element_text(hjust=1, angle=45, size=13), 
        axis.text.y = element_text(size=16), 
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(size=16),
        strip.text = element_text(size=10))
)

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_alt_line_plot.png", abc_alt_line_plot, height=5.5, width=6.5, bg="white")


# longitudinal fever / symptoms ####

#collate all temperatures recorded as pyrexia/fever from the symptoms/AE tab

# vac69a
temp_data <- read.csv("~/PhD/clinical_data/vac69a/vac69a_body_temp.csv")

fever <- temp_data %>%
  mutate("volunteer" = gsub("69010", "v", trial_number)) %>%
  filter(study_event_oid != "SE_VA69_C28") %>%
  select(volunteer, timepoint, temp) %>%
  #filter(temp>37.5)
  filter(timepoint %in% c("03_chall___01_am", "14_chall___07_am", "24_chall___12_am","28_chall___14_am",
                          "56_diagnosis_or_c_21", "58_postchall_ep1", "59_postchall_ep2", "62_postchall_t___6"))


lousy_fever_timepoints <- c("03_chall___01_am", "14_chall___07_am", "24_chall___12_am","28_chall___14_am",
                            "56_diagnosis_or_c_21", "58_postchall_ep1", "59_postchall_ep2", "62_postchall_t___6")

good_fever_timepoints <- c("Baseline", "C7", "C12","C14", "Diagnosis", "T1", "T2", "T6")

fever_timepoint_replacement <- setNames(good_fever_timepoints, lousy_fever_timepoints)
fever$timepoint <- stringr::str_replace_all(fever$timepoint, fever_timepoint_replacement)


symptom_data <- read.csv("~/PhD/clinical_data/vac69a/symptoms_vac69a.csv", header = T, stringsAsFactors = F)
symptom_data$flo_timepoint <- gsub("am ", "am", symptom_data$flo_timepoint)
symptom_data$flo_timepoint <- gsub("pm ", "pm", symptom_data$flo_timepoint)
symptom_data$flo_timepoint <- gsub("DoD ", "Diagnosis", symptom_data$flo_timepoint)
symptom_data$timepoint <- symptom_data$flo_timepoint

symptom_data$volunteer <- gsub("69010", "v", symptom_data$trial_number)


pyrexia <- symptom_data %>%
  mutate("temp"=pyrexia_temp) %>%
  select(volunteer, timepoint, temp) 

pyrexia <- na.omit(pyrexia)

pyrexia$timepoint <- c("T1", "C12", "T1", "C12", "Diagnosis", "T2", "C11", "Diagnosis", "T1", "Diagnosis")

fever <- rbind(fever, pyrexia)

fever$n_infection <- "First"

real_fever <- fever %>%
  group_by(volunteer, timepoint) %>%
  summarise("max_temp"=max(temp)) %>%
  ungroup() %>%
  group_by(volunteer) %>%
  summarise("max_temp"=max(max_temp))



fever <- fever %>%
  group_by(volunteer, timepoint) %>%
  summarise("max_temp"=max(temp)) %>%
  filter(timepoint %in% good_fever_timepoints)

# 


# 1 v02           38.5
# 2 v03           39.7
# 3 v05           37.1
# 4 v06           37.9
# 5 v07           38.7
# 6 v09           36.9

vac69a_fever <- subset(fever, fever$max_temp>37.5)
vac69a_fever$n_infection <- "First"

vac69b_data <- read.csv("~/PhD/clinical_data/vac69b/symptoms.csv")

vac69b_data <- subset(vac69b_data, !is.na(vac69b_data$pyrexia_temp))

vac69b_fever <- data.frame("volunteer" = gsub("69010", "v", vac69b_data$trial_number),
                           "timepoint" = c("T1", "T2", "Diagnosis", "T1", "T2"),
                           "max_temp"= vac69b_data$pyrexia_temp,
                           "n_infection"=ifelse(vac69b_data$trial_number=="6901007", "Second", "First"))


vac69c_data <- read.csv("~/PhD/clinical_data/vac69c/symptoms.csv")

vac69c_data <- subset(vac69c_data, !is.na(vac69c_data$pyrexia_temp))

vac69c_fever <- data.frame("volunteer" = gsub("69010", "v", vac69c_data$trial_number),
                           "timepoint" = c("Diagnosis", "T1", "C15", "T1", "T1", "Diagnosis", "C19", "T1", "T1", "T1"),
                           "max_temp"= vac69c_data$pyrexia_temp,
                           "n_infection"=ifelse(vac69c_data$trial_number=="6901007", "Third",
                                                ifelse(vac69c_data$trial_number %in% c("6901011", "6901021"), "Second", "First"))
)




vac69d_data <- read.csv("~/PhD/clinical_data/vac69d/vac69d_symptoms.csv")

vac69d_data <- subset(vac69d_data, !is.na(vac69d_data$pyrexia_temp))

vac69d_fever <- data.frame("volunteer" = gsub("69010", "v", vac69d_data$trial_number),
                           "timepoint" = c("C15", "T1", "C15", "T1", "C17", "T1", "C15", "Diagnosis"),
                           "max_temp"= vac69d_data$pyrexia_temp)

vac69d_fever$n_infection <-  ifelse(vac69d_fever$volunteer %in% c("v11", "v21"), "Third",
                      ifelse(vac69d_fever$volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First")
)




abc_fever <- rbind(vac69a_fever, vac69b_fever, vac69c_fever, vac69d_fever)

abc_max_fever <- abc_fever %>%
  group_by(volunteer, n_infection) %>%
  summarise("max_temperature"=max(max_temp))



#collate all information from the observations tab, which includes all non-fever temperatures taken
# at clinic. we then combine this with the fever temperatures and pick the highest value for each volunteer

vac69a_temp <- read.csv("~/PhD/clinical_data/vac69a/vac69a_observations.csv")
vac69a_temp <- vac69a_temp %>%
  group_by(trial_number) %>%
  summarise("max_temperature"=max(temp))
vac69a_temp$volunteer <-gsub("69010", "v", vac69a_temp$trial_number)
vac69a_temp$n_infection <- "First"


vac69b_temp <- read.csv("~/PhD/clinical_data/vac69b/vac69b_observations.csv")
vac69b_temp <- vac69b_temp %>%
  group_by(trial_number) %>%
  summarise("max_temperature"=max(temp))
vac69b_temp$volunteer <-gsub("69010", "v", vac69b_temp$trial_number)
vac69b_temp$n_infection <- ifelse(vac69b_temp$volunteer %in% c("v11", "v21"), "First", "Second")

vac69c_temp <- read.csv("~/PhD/clinical_data/vac69c/vac69c_observations.csv")
vac69c_temp <- vac69c_temp %>%
  group_by(trial_number) %>%
  summarise("max_temperature"=max(sc_temp, na.rm = TRUE))
vac69c_temp$volunteer <- gsub("69010", "v", vac69c_temp$trial_number)
vac69c_temp$n_infection <- ifelse(vac69c_temp$volunteer=="v07", "Third", ifelse(vac69c_temp$volunteer %in% c("v11", "v21"), "Second", "First")
)

vac69d_temp <- read.csv("~/PhD/clinical_data/vac69d/vac69d_observations.csv")
vac69d_temp <- vac69d_temp %>%
  group_by(trial_number) %>%
  summarise("max_temperature"=max(sc_temp, na.rm = TRUE))
vac69d_temp$volunteer <- gsub("69010", "v", vac69d_temp$trial_number)
vac69d_temp$n_infection <-  ifelse(vac69d_temp$volunteer %in% c("v11", "v21"), "Third",
                                   ifelse(vac69d_temp$volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First")
                                   )               


abc_temp <- rbind(vac69a_temp, vac69b_temp, vac69c_temp, vac69d_temp)

abc_max_temp <- abc_temp %>%
  group_by(volunteer, n_infection)%>%
  top_n(1, max_temperature)

fever_temp_combo <- rbind(abc_max_temp, abc_max_fever)

fever_temp_combo <- fever_temp_combo %>%
  group_by(volunteer, n_infection) %>%
  slice_max(order_by = max_temperature, with_ties = FALSE)


degrees_celsius <- scales::label_number(
  decimal.mark = ".",
  accuracy = 0.1,
  suffix = "°C"
)

fever_temp_combo <- filter(fever_temp_combo, n_infection != "Third")

abc_fever_plot <- ggplot(fever_temp_combo, aes(x=n_infection,
                             y=max_temperature))+
  geom_boxplot(aes(fill=n_infection), outlier.alpha = 0)+
  geom_dotplot(aes(fill=volunteer),
               binaxis = "y",
               colour="black",
               stackdir = "center",
               dotsize = 2,
               binpositions = "all",
               stackgroups = TRUE,
               binwidth = 0.03)+
  #geom_point(aes(color=volunteer),fill="white", stroke=1, shape=21, size=3, position = position_dodge2(width=0.2))+
  #geom_line(aes(color=volunteer, group=volunteer), size=0.9)+
  geom_hline(yintercept = 37.5, linetype="dashed", color="black")+
  #ylab(expression(paste("Temperature (",degree,"C)",sep="")))+
  ylab("Temperature")+
  # ggtitle(expression(paste("reduced fevers in \n",
  #                          italic("P. vivax"), " reinfection", sep='')))+
  scale_y_continuous(labels = degrees_celsius, breaks = seq(36, 40, by=0.5))+
  scale_fill_manual(values=c(n_infection_palette, scales::hue_pal()(19)))+
  theme_minimal()+
  theme(legend.position = "none", 
    #plot.title = element_text(hjust=0.5, vjust = 0, size=13),
    #plot.caption = element_custom(),
  
        axis.text.x = element_text(hjust=1, angle=45, size=11),
        axis.title.x = element_blank())

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_fever.png", abc_fever_plot, height=5, width=4, bg="white")

alt_theme <- theme(
      axis.text.x = element_text(hjust=1, angle=45, size=13), 
      axis.text.y = element_text(size=16), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=16),
      strip.text = element_text(size=10))


#
ggsave("~/PhD/clinical_data/vac69c/figures/abcd_combo_fever_alt.png", abc_fever_plot+alt_theme, height=5, width=4, bg="white")

# parasitaemias ####

# vac69a

vac69a_parasitaemia <- read.csv("~/PhD/clinical_data/vac69a/parasitaemia/better_vac69a_parasitaemia_ultimate_qc.csv", header=T)

vac69a_parasitaemia <- vac69a_parasitaemia %>%
  pivot_longer(colnames(vac69a_parasitaemia)[5:ncol(vac69a_parasitaemia)], names_to = "Timepoint", values_to = "Genomes") %>%
  select(Volunteer, Timepoint, Genomes) %>%
  mutate("n_infection"="First")

vac69a_parasitaemia$Timepoint <- gsub("D0", "Baseline", vac69a_parasitaemia$Timepoint)
vac69a_parasitaemia$Timepoint <- gsub("D", "C", vac69a_parasitaemia$Timepoint, fixed=T)

vac69a_parasitaemia$Timepoint <- ifelse(!grepl(".5", vac69a_parasitaemia$Timepoint, fixed = T),
                                        paste(vac69a_parasitaemia$Timepoint, "am"), 
                                        paste(vac69a_parasitaemia$Timepoint, "pm"))

vac69a_parasitaemia$Timepoint <- gsub(".5", "", vac69a_parasitaemia$Timepoint, fixed=T)

# vac69b

vac69b_parasitaemia <- read.csv("~/PhD/clinical_data/vac69b/parasitaemia/vac69b_parasitaemia.csv", header=T)
vac69b_parasitaemia <- subset(vac69b_parasitaemia, Vaccination == "Control") 

vac69b_parasitaemia <- vac69b_parasitaemia %>%
  pivot_longer(colnames(vac69b_parasitaemia)[5:ncol(vac69b_parasitaemia)], names_to = "Timepoint", values_to = "Genomes") %>%
  select(Volunteer, Timepoint, Genomes) %>%
  mutate("n_infection"=ifelse(.$Volunteer %in% c("v11", "v21"), "First", "Second"))

bad_names <- c("v7", "v2$", "v5", "v11", "v21")
good_names <- c("v07", "v02", "v05", "v11", "v21")

name_replacement <- setNames(good_names, bad_names)

vac69b_parasitaemia$Volunteer <- stringr::str_replace_all(vac69b_parasitaemia$Volunteer, name_replacement)

vac69b_parasitaemia$Timepoint <- ifelse(grepl(".5", vac69b_parasitaemia$Timepoint, fixed = T), paste(vac69b_parasitaemia$Timepoint, "pm"),  paste(vac69b_parasitaemia$Timepoint, "am"))
vac69b_parasitaemia$Timepoint <- gsub(".5", "", vac69b_parasitaemia$Timepoint, fixed=T)

# vac69c

vac69c_parasitaemia <- read.csv("~/PhD/clinical_data/vac69c/vac69c_unvaccinated_qpcr.csv", header=T)
vac69c_parasitaemia$X <- NULL
vac69c_parasitaemia$Volunteer <-  gsub("69010", "v",vac69c_parasitaemia$trial_code)

vac69c_parasitaemia <- vac69c_parasitaemia %>%
  pivot_longer(colnames(vac69c_parasitaemia)[3:(ncol(vac69c_parasitaemia)-1)], names_to = "Timepoint", values_to = "Genomes")%>%
  select(Volunteer, Timepoint, Genomes) %>%
  mutate("n_infection"= ifelse(.$Volunteer %in% c("v11", "v21"), "Third",
                               ifelse(.$Volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First"))
  )



vac69c_parasitaemia$Timepoint <- ifelse(grepl(".5", vac69c_parasitaemia$Timepoint, fixed = T), paste(vac69c_parasitaemia$Timepoint, "pm"),  paste(vac69c_parasitaemia$Timepoint, "am"))
vac69c_parasitaemia$Timepoint <- gsub(".5", "", vac69c_parasitaemia$Timepoint, fixed=T)
vac69c_parasitaemia$Timepoint <- gsub(".", "", vac69c_parasitaemia$Timepoint, fixed=T)

vac69c_parasitaemia$n_infection <- ifelse(vac69c_parasitaemia$Volunteer == "v07", "Third",
       ifelse(vac69c_parasitaemia$Volunteer %in% c("v11", "v21"), "Second", "First"))

# vac69d

vac69d_parasitaemia <- read.csv("~/PhD/clinical_data/vac69d/vac69d_parasitaemias.csv", header=T)
vac69d_parasitaemia$X <- NULL
vac69d_parasitaemia$Volunteer <-  gsub("69010", "v", vac69d_parasitaemia$trial_number)

vac69d_parasitaemia <- vac69d_parasitaemia %>%
  pivot_longer(colnames(vac69d_parasitaemia)[3:(ncol(vac69d_parasitaemia)-1)], names_to = "Timepoint", values_to = "Genomes")%>%
  select(Volunteer, Timepoint, Genomes) %>%
  mutate("n_infection"= ifelse(.$Volunteer %in% c("v11", "v21"), "Third",
                               ifelse(.$Volunteer %in% c("v34", "v35", "v37", "v39", "v40", "v42", "v44"), "Second", "First"))
  )
         


vac69d_parasitaemia$Timepoint <- ifelse(grepl(".5", vac69d_parasitaemia$Timepoint, fixed = T), paste(vac69d_parasitaemia$Timepoint, "pm"),  paste(vac69d_parasitaemia$Timepoint, "am"))
vac69d_parasitaemia$Timepoint <- gsub(".5", "", vac69d_parasitaemia$Timepoint, fixed=T)
vac69d_parasitaemia$Timepoint <- gsub(".", "", vac69d_parasitaemia$Timepoint, fixed=T)

combo_parasitaemia <- rbind(vac69a_parasitaemia, vac69b_parasitaemia, vac69c_parasitaemia, vac69d_parasitaemia)

combo_parasitaemia$sample_id <- paste(combo_parasitaemia$Volunteer, combo_parasitaemia$n_infection)
combo_parasitaemia$Timepoint <- gsub("D", "C", combo_parasitaemia$Timepoint)



# get rid of ugly timepoints that nobody wants
parasitaemia_levels <- unique(gtools::mixedsort(combo_parasitaemia$Timepoint))
combo_parasitaemia <- subset(combo_parasitaemia, Timepoint %in% parasitaemia_levels[1:56])

# get rid of post treatment timepoints: what this code does is assign an integer value to each timepoint like so:
combo_parasitaemia$FTimepoint <- as.numeric(factor(combo_parasitaemia$Timepoint, levels=c(parasitaemia_levels)))

# turn our tibble into a dataframe, not sure why this is important, but it is
combo_parasitaemia <- data.frame(combo_parasitaemia)

# split the dataframe according to each individual infection
list_of_curves <- split(combo_parasitaemia, combo_parasitaemia$sample_id)
# we make a dataframe that contains the timpoint of maximum parasitaemia and all associated data
list_of_max <- data.frame(t(sapply(list_of_curves, function(x) subset(x, Genomes==max(x$Genomes, na.rm = TRUE))[1,])))

# this monster cycles through the list of df's for each infection and subsets so at or previous to the timepoint of maximum parasitaemia 
# are retained
shorter_list_of_curves <- lapply(list_of_curves, function(x){
  subset(x, x$FTimepoint <= unique(
                              unname(
                                unlist(list_of_max$FTimepoint[match(x$sample_id, names(list_of_curves))]
                                       )
                                )
  )
  )
})

#put the split data frames back together as 1, like Jah would have wanted
pre_drug_parasitaemia <- do.call(rbind, shorter_list_of_curves)
pre_drug_parasitaemia$Timepoint <- gsub("am", "", pre_drug_parasitaemia$Timepoint)

parasitaemia_levels <- unique(gtools::mixedsort(pre_drug_parasitaemia$Timepoint))

# bam
vac69abcd_parasitaemia_plot <- ggplot(pre_drug_parasitaemia[!is.na(pre_drug_parasitaemia$Genomes),], aes(x=factor(Timepoint, levels=parasitaemia_levels), y=Genomes+1, group=sample_id))+
  geom_point(aes(color=n_infection))+
  geom_line(aes(color=n_infection))+
  scale_y_log10()+
  scale_color_manual(values=n_infection_palette)+
  theme_minimal()+
  scale_x_discrete(breaks=parasitaemia_levels[!grepl("*pm*", parasitaemia_levels)][-2])+
  ylab("Genomes / mL")+
  xlab("Timepoint")+
  #guides(color=guide_legend(ncol=2))+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        legend.text = element_text(size=17),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18))

ggsave("~/PhD/clinical_data/vac69c/figures/abcd_parasitaemia.png", vac69abcd_parasitaemia_plot, height = 6, width = 8, bg="white", dpi=444)
  