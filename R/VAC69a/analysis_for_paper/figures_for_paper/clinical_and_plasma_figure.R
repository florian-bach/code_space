  # Panel A clinical presentation of volunteers undergoing vivax CHMI: plateltets, fever lymphocytes, adverse events, parasitaemia ####
    
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    `%notin%` <- Negate(`%in%`)
    
    fig1_theme <- theme(axis.title.x = element_blank(),
                        legend.title = element_text(size = 9), 
                        legend.text = element_text(size = 9),
                        axis.title=element_text(size=10))
    
    
    volunteer_colours <- list("v02" = "#FB9A99",
                              "v03" = "#E31A1C",
                              "v05" = "#A6CEE3",
                              "v06" = "#1F78B4",
                              "v07" = "#F0E442",
                              "v09" = "#E69F00")
    
    
    volunteer_palette <- unlist(unname(volunteer_colours))
    names(volunteer_palette) <- names(volunteer_colours)
    
    
    
    # biochem ####
    
    
    
    biochem_data <- read.csv("~/PhD/clinical_data/vac69a/biochem.csv")
    
    colnames(biochem_data)[1] <- "volunteer"
    
    biochem_data$volunteer <- paste("v", substr(biochem_data$volunteer, 6, 7), sep='')
    
    
    
    lousy_timepoints <- unique(as.character(biochem_data$timepoint))
    
    good_timepoints <- c("C7 am", "C14 am", "Screening","Baseline", "Diagnosis", "C28", "T1", "T6", "C90", "EV")
    
    biochem_timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)
    biochem_data$timepoint <- stringr::str_replace_all(biochem_data$timepoint, biochem_timepoint_replacement)
    
    
    
    #biochem_data <- filter(biochem_data, timepoint %in% c("_C_1", "_T6", "_EP"))
  
    line_plot_data <- select(biochem_data,  volunteer, timepoint, bilirubin, alt, alkphos, albumin)
    
    long_line_plot_data <- line_plot_data %>%
      filter(timepoint %in% c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")) %>%
      gather(Analyte, Concentration, c(bilirubin, alt, alkphos, albumin))
    
    
    long_line_plot_data$Analyte <- paste(
      toupper(substr(long_line_plot_data$Analyte, 1,1)),
      substr(long_line_plot_data$Analyte, 2,nchar(long_line_plot_data$Analyte)),
      sep="")
    
    
    
      biochem_line_plot <- ggplot(long_line_plot_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Concentration, color=volunteer, group=volunteer))+
      scale_fill_manual(values=volunteer_palette)+
      scale_color_manual(values=volunteer_palette)+
      geom_line(aes(color=volunteer), size=0.9)+
      geom_point(fill="white", stroke=1, shape=21, size=0.9)+
      theme_minimal()+
      facet_wrap(~Analyte, scales="free", ncol = 4)+
      xlab("Timepoint")+
      #ylab(expression(Cells~"/"~mu*L~blood))+
      guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
      fig1_theme+
      #scale_y_continuous(label=scales::comma)+
      theme(plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(hjust=1, angle=45, size=8), 
            axis.title = element_blank(),
            legend.position = "none",
            strip.text = element_text(size=10))
    
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/biochem_line_plot.pdf", biochem_line_plot, width=6.4, height=2)
    
    
    
    biochem_data <- select(biochem_data,  volunteer, timepoint, grep("_ae", colnames(biochem_data), fixed = TRUE, value = TRUE))
    
    
    lousy_timepoints <- unique(as.character(biochem_data$timepoint))
    
    good_timepoints <- c("C7 am", "C14 am", "Screening","Baseline", "Diagnosis", "C28", "T1", "T6", "C90", "EV")
    
    biochem_timepoint_replacement <- setNames(good_timepoints, lousy_timepoints)
    biochem_data$timepoint <- stringr::str_replace_all(biochem_data$timepoint, biochem_timepoint_replacement)
    
    long_biochem_data <- gather(biochem_data, symptom, severity, grep("_ae", colnames(biochem_data), fixed = TRUE, value = TRUE))
    long_biochem_data$symptom <- gsub("_ae", "", long_biochem_data$symptom)
    
    
    
    
    
    long_biochem_data$symptom <- paste(
      toupper(substr(long_biochem_data$symptom, 1,1)),
      substr(long_biochem_data$symptom, 2,nchar(long_biochem_data$symptom)),
      sep="")
    
    
    # symptoms ####
    symptom_data <- read.csv("~/PhD/clinical_data/vac69a/symptoms_vac69a.csv", header = T, stringsAsFactors = F)
    symptom_data$flo_timepoint <- gsub("am ", "am", symptom_data$flo_timepoint)
    symptom_data$flo_timepoint <- gsub("pm ", "pm", symptom_data$flo_timepoint)
    symptom_data$flo_timepoint <- gsub("DoD ", "Diagnosis", symptom_data$flo_timepoint)
    symptom_data$timepoint <- symptom_data$flo_timepoint
    
    symptom_data$volunteer <- gsub("69010", "v", symptom_data$trial_number)
      
    
    
    long_symptom_data <- tidyr::gather(symptom_data, symptom, severity, colnames(symptom_data)[c(12, 15:ncol(symptom_data)-1)])
    
    long_symptom_data$timepoint <- as.character(symptom_data$flo_timepoint)
    
    
    long_symptom_data <- symptom_data %>%
      gather(symptom, severity, colnames(symptom_data)[c(12, 14:ncol(symptom_data)-1)]) %>%
      select(volunteer, timepoint, symptom, severity)
    
    long_symptom_data$symptom <- paste(
      toupper(substr(long_symptom_data$symptom, 1,1)),
      substr(long_symptom_data$symptom, 2,nchar(long_symptom_data$symptom)),
      sep="")
    
    
    
    # haematology ####
    
    haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")
    
    haem_data$volunteer <- gsub("69010", "v", haem_data$trial_number)
    
    long_haem_data <- haem_data %>%
      select(volunteer, timepoint, grep("_ae", colnames(haem_data), fixed = TRUE, value = TRUE))%>%
      gather(symptom, severity, grep("_ae", colnames(haem_data), fixed = TRUE, value = TRUE))
    
    long_haem_data <- na.omit(long_haem_data)
    
    lousy_haem_timepoints <- unique(long_haem_data$timepoint)
    good_haem_timepoints <- c("C90", "T6", "C7 am", "Baseline", "C28", "EV", "C14 am", "Diagnosis", "Screening", "T1")
    
    haem_time_dic <- setNames(good_haem_timepoints, lousy_haem_timepoints)
    
    long_haem_data$timepoint <- stringr::str_replace_all(long_haem_data$timepoint, haem_time_dic)
    
    long_haem_data$severity <- as.character(long_haem_data$severity)
    long_haem_data$symptom <- substr(long_haem_data$symptom, 1,nchar(long_haem_data$symptom)-3)
    
    long_haem_data$symptom <- paste(
      toupper(substr(long_haem_data$symptom, 1,1)),
      substr(long_haem_data$symptom, 2,nchar(long_haem_data$symptom)),
      sep="")
    
    long_haem_data$symptom <- gsub("Wbc", "WBC", long_haem_data$symptom)
    
    
    # putting it all together ####
    
    long_haem_data$source <- "Haematology"
    long_biochem_data$source <- "Biochemistry"
    long_symptom_data$source <- "Symptom"
    
    long_ae_data <- rbind(long_haem_data, long_biochem_data, long_symptom_data)
    
    long_ae_data$timepoint <- gsub("Baseline ", "Baseline", long_ae_data$timepoint)
    long_ae_data$timepoint <- gsub("DoD", "Diagnosis", long_ae_data$timepoint)
    
    long_ae_data$timepoint <- gsub("T6 ", "T6", long_ae_data$timepoint)
    long_ae_data$timepoint <- gsub("T1 ", "T1", long_ae_data$timepoint)
    long_ae_data$timepoint <- gsub("T2 ", "T2", long_ae_data$timepoint)
    long_ae_data$timepoint <- gsub("C28 ", "C28", long_ae_data$timepoint)
    
    long_ae_data$severity <- as.numeric(as.character(long_ae_data$severity))
    
    adverse_events <- long_ae_data %>%
      filter(severity >= 1) %>%
      filter(symptom != "Pyrexia_temp") %>%
      group_by(volunteer, timepoint, severity) %>%
      summarise(ae_count = n())
    
    
    timepoints <- unique(long_ae_data$timepoint[gtools::mixedorder(long_ae_data$timepoint)])
    
    timepoint_levels <- timepoints[c(1:31, 34, 37:length(timepoints))]
    
    
    adverse_events <- filter(adverse_events, timepoint %in% timepoint_levels)
    
    (all_ae_stack <- ggplot(adverse_events,  aes(x=factor(timepoint, levels=timepoint_levels), y=ae_count/6, fill=factor(severity, levels=rev(1:3))))+
        geom_bar(stat="identity", position = "stack")+
        scale_fill_manual(values =  list("1"="#FACA0F", "2"="chocolate1", "3"="red"))+
        ylab("Mean # of AEs per Volunteer")+
        xlab("Timepoint")+
        ggtitle("All Adverse Events")+
        geom_vline(xintercept = 23.5)+
        scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by=2))+
        theme_minimal()+
        fig1_theme+
        theme(axis.text.x = element_text(hjust=1, angle=45, size=5),
              plot.title = element_text(hjust=0.5, vjust=0, size=10),
              panel.grid.minor = element_blank())+
        guides(fill=guide_legend(title="Severity",
                                 override.aes = list(size = 0.1),
                                 keywidth = 0.5,reverse = T,
                                 keyheight = 0.5)))
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/all_adverse_events_stacked.pdf", all_ae_stack, height=4, width=6)
    
    
    symptom_plot_data <- long_symptom_data %>%
      filter(severity > 0) %>%
      group_by(volunteer, timepoint, severity) %>%
      summarise(ae_count = n())
    
    
    symptom_timepoints <- unique(symptom_plot_data$timepoint[gtools::mixedorder(symptom_plot_data$timepoint)])
    
    
    (symptom_stack <- ggplot(symptom_plot_data,  aes(x=factor(timepoint, levels=symptom_timepoints), y=ae_count/6, fill=as.character(severity)))+
        geom_bar(stat="identity", position = "stack")+
        scale_fill_manual(values =  list("1"="#FACA0F", "2"="chocolate1", "3"="red"))+
        ylab("Mean # of AEs\nper Volunteer")+
        xlab("Timepoint")+
        ggtitle("All Adverse Events")+
        geom_vline(aes(xintercept = 22.5))+
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
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/symptom_adverse_events_stacked.pdf", symptom_stack, height=4, width=6)
    
    
    
    # max ae score ####
    
    dods <- data.frame(vol=c('v02', 'v03', 'v05', 'v06', 'v07', 'v09'),
                       paras=c(4907, 8054, 16733, 7464, 21870, 15051),
                       dod=c('D15.5', 'D12.5', 'D15.5', 'D15.5', 'D16', 'D17'))
    
    
    v02_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v02", timepoint %in% timepoint_levels[c(26:29, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v02 <- nrow(v02_adverse_events)
    
    
    v03_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v03", timepoint %in% timepoint_levels[c(20:23, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v03 <- nrow(v03_adverse_events)
    
    
    v05_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v05", timepoint %in% timepoint_levels[c(26:29, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v05 <- nrow(v05_adverse_events)
    
    
    v06_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v06", timepoint %in% timepoint_levels[c(26:29, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v06 <- nrow(v06_adverse_events)
    
    
    v07_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v07", timepoint %in% timepoint_levels[c(27:30, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v07 <- nrow(v07_adverse_events)
    
    
    v09_adverse_events <- long_ae_data %>%
      filter(severity>0, volunteer=="v09", timepoint %in% timepoint_levels[c(28:31, 32:34)]) %>%
      group_by(timepoint) %>%
      summarise(ae_count = n())
    
    max_ae_v09 <- nrow(v09_adverse_events)
    
    
    
    
    # ae plots ####
    
    
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
    
    
    # fever
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
    
     
    pyrexia <- symptom_data %>%
      mutate("temp"=pyrexia_temp) %>%
      select(volunteer, timepoint, temp) 
      
    pyrexia <- na.omit(pyrexia)
    
    pyrexia$timepoint <- c("T1", "C12", "T1", "C12", "Diagnosis", "T2", "C11", "Diagnosis", "T1", "Diagnosis")
    
    fever <- rbind(fever, pyrexia)
    
    
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
    
    
    
    
    
    fever_curves <- ggplot(fever, aes(x=factor(timepoint, levels=good_fever_timepoints), y=max_temp, color=volunteer, group=volunteer))+
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
      guides(color=guide_legend(title = "Volunteer"))+
      theme(plot.title = element_text(hjust=0.5, vjust = 0, size=10),
            axis.text.x = element_text(hjust=1, angle=45, size=8),
            legend.position = "right", 
            axis.title.x = element_blank())
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fever_curves.pdf", fever_curves, width=2, height=2)
    
    
    # PLatelets/Lymphocytes
    
    
    haem_data <- data.table::fread("~/PhD/clinical_data/vac69a/haem.csv")
    
    haem_data$volunteer <- gsub("69010", "v", haem_data$trial_number)
    
    long_haem_data <- haem_data %>%
      select(volunteer, timepoint, platelets, lymphocytes) %>%
      #select(trial_number, timepoint, lymphocytes) %>%
      filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")) %>%
      gather(Cell, Frequency, c(platelets,lymphocytes))
    
    
    
    
    
    supp_haem_data <- haem_data %>%
      select(volunteer, timepoint, wbc, neutrophils,  monocytes, eosinophils, haemaglobin) %>%
      #select(trial_number, timepoint, lymphocytes) %>%
      filter(timepoint %in% c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")) %>%
      gather(Cell, Frequency, c(wbc, neutrophils,  monocytes, eosinophils, haemaglobin))
    
  
    supp_haem_data$Cell <- paste(
      toupper(substr(supp_haem_data$Cell, 1,1)),
      substr(supp_haem_data$Cell, 2,nchar(supp_haem_data$Cell)),
      sep="")
    
    supp_haem_data$Cell <- factor(supp_haem_data$Cell, levels=c("Wbc", "Neutrophils",  "Monocytes", "Eosinophils", "Haemaglobin"))
    
  
    bad_timepoints <- c("_C_1", "_C1_7", "_C8_14", "_C21", "_EP", "_T6", "_C90")
    great_timepoints <- c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")
    
    time_dic <- setNames(great_timepoints, bad_timepoints)
    
    long_haem_data$timepoint <- stringr::str_replace_all(long_haem_data$timepoint, time_dic)
    supp_haem_data$timepoint <- stringr::str_replace_all(supp_haem_data$timepoint, time_dic)
    
    
    supp_haem_data$Cell <- gsub("Haemaglobin", "Haemoglobin", supp_haem_data$Cell)
    
    haemoglobin_data <-  filter(supp_haem_data, Cell == "Haemoglobin")
    supp_haem_data <- filter(supp_haem_data, Cell != "Haemoglobin")
    
    
    
    lymphsss <- filter(long_haem_data, Cell=="lymphocytes")
    
    thrombos_lymphs <- ggplot(long_haem_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
      scale_fill_manual(values=volunteer_palette)+
      scale_color_manual(values=volunteer_palette)+
      geom_line(aes(color=volunteer), size=0.9)+
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
    
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/thrombos_lymphs.pdf", thrombos_lymphs, width=6, height=2.2)
    
    
    
    lymph_plot <- ggplot(lymphsss, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
      scale_fill_manual(values=volunteer_palette)+
      scale_color_manual(values=volunteer_palette)+
      geom_line(aes(color=volunteer), size=0.9)+
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
    
    
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lymphs.pdf", lymph_plot, width=3, height=2.2)
    
    
    supp_haem_data_plots <- ggplot(supp_haem_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency*1000, color=volunteer, group=volunteer))+
      scale_fill_manual(values=volunteer_palette)+
      scale_color_manual(values=volunteer_palette)+
      geom_line(aes(color=volunteer), size=0.9)+
      geom_point(fill="white", stroke=1, shape=21, size=0.9)+
      theme_minimal()+
      facet_wrap(~Cell, scales="free", nrow = 1)+
      xlab("Timepoint")+
      ylab(expression(Cells~"/"~mu*L~blood))+
      guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
      fig1_theme+
      scale_y_continuous(label=scales::comma)+
      theme(plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(hjust=1, angle=45, size=8), 
            axis.title.x = element_blank(),
            legend.position = "right",
            strip.text = element_text(size=10))
    
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/supp_haem_counts.pdf", supp_haem_data_plots, height=2, width=8)
    
    
    haemoglobin_data_plots <- ggplot(haemoglobin_data, aes(x=factor(timepoint, levels=c("Baseline", "C7 am", "C14 am", "Diagnosis", "T1", "T6", "C90")), y=Frequency, color=volunteer, group=volunteer))+
      scale_fill_manual(values=volunteer_palette)+
      scale_color_manual(values=volunteer_palette)+
      geom_line(aes(color=volunteer), size=0.9)+
      geom_point(fill="white", stroke=1, shape=21, size=0.9)+
      facet_grid(~Cell)+
      theme_minimal()+
      xlab("Timepoint")+
      ylab("mg / mL")+
      guides(color=guide_legend(title="Volunteer", override.aes = list(size=1)))+
      fig1_theme+
      scale_y_continuous(label=scales::comma)+
      theme(plot.title = element_text(hjust=0.5),
            axis.text.x = element_text(hjust=1, angle=45, size=8), 
            axis.title.x = element_blank(),
            legend.position = "none",
            strip.text = element_text(size=10))
    
    
    ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/haemoglobin_data_plots.pdf", haemoglobin_data_plots, height = 2, width=)
    
  #create max lymphopenia metric
  
  # wide_haem_data <- spread(long_haem_data, timepoint, Frequency)
  # wide_haem_data[,3:8] <- wide_haem_data[,3:8]/wide_haem_data[,3]
  # 
  # haem_fc <- wide_haem_data %>%
  #   gather(Timepoint, FC, colnames(wide_haem_data)[3:8]) %>%
  #   filter(Timepoint != "T6") %>%
  #   group_by(volunteer) %>%
  #   drop_na()%>%
  #   summarise("max_lymphopenia"=min(FC))
  
  # volunteer max_lymphopenia
  # <chr>               <dbl>
  #   1 v02                 0.389
  # 2 v03                 0.399
  # 3 v05                 0.306
  # 4 v06                 0.169
  # 5 v07                 0.345
  # 6 v09                 0.307
  
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
  great_timepoints <- c("Baseline", "C7", "C14", "T1", "T6")
  
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
    
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/lymphocyte_ae_pie.pdf", lymphocyte_ae_pie, height=2, width=3)
    
    
    adverse_pies <- long_ae_data %>%
      mutate(timepoint = gsub(" pm", "", timepoint, fixed = T)) %>%
      mutate(timepoint = gsub(" am", "", timepoint, fixed = T)) %>%
      mutate(symptom = gsub("Haemaglobin", "Haemoglobin", symptom, fixed = T)) %>%
      filter(timepoint %in% c("Diagnosis", "T1","T2", "T6")) %>%
      filter(symptom %notin% c("vomiting", "nausea","pyrexia", "diarrhoea", "back_pain", "arthralgia",
                               "Ggt", "Ast", "Albumin", "Alkphos", "Pyrexia_temp", "Bilirubin", "Creatinin")) 
    
    
  for(i in unique(adverse_pies$source)){  
    
    
    plottable_data <- filter(adverse_pies, source==i)
    plottable_data$severity <- factor(plottable_data$severity, levels=paste(seq(0,3)))
    
    ae_pies <- ggplot(plottable_data,  aes(x="", fill=severity))+
      geom_bar(stat="count", position ="stack")+
      scale_fill_manual(values =  list("0"="lightgrey", "1"="#FACA0F", "2"="chocolate1", "3"="red"))+
      facet_grid(timepoint~symptom, switch="y")+
      coord_polar(theta = "y")+
      ggtitle(paste(i, "Adverse Events\n"))+
      theme_void()+
      fig1_theme+
      theme(plot.title = element_text(hjust=0.5, vjust=0, size=10),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.title = element_blank(),
            strip.text.x = element_text(hjust=0.5, size=6, angle = 0),
            strip.text.y.left = element_text(hjust=0.5, size=6, angle = 0))#+
      # guides(fill=guide_legend(title="Severity",
      #                          override.aes = list(size = 0.1),
      #                          keywidth = 0.5,
      #                          keyheight = 0.5))
    
    ggsave(paste("~/PhD/cytof/vac69a/final_figures_for_paper/", i, "_ae_pies.pdf", sep=''), ae_pies, height=3.5, width=0.75*length(unique(plottable_data$symptom)))
  }
  
  
      # parasitaemias ####
  
  data <- read.csv("~/PhD/clinical_data/vac69a/parasitaemia/better_vac69a_parasitaemia_ultimate_qc.csv", header=T)
  parasitaemias <- gather(data, Timepoint, Genomes, colnames(data)[5:ncol(data)])
  parasitaemias$Genomes <- as.numeric(parasitaemias$Genomes)
  # get rid of garbage timepoints that mess up graph
  
  parasitaemias$Timepoint <- ifelse(grepl(".5", parasitaemias$Timepoint, fixed = T), paste(parasitaemias$Timepoint, "pm"),  paste(parasitaemias$Timepoint, "am"))
  parasitaemias$Timepoint <- gsub(".5", "", parasitaemias$Timepoint, fixed=T)
  parasitaemias$Timepoint <- gsub("D0", "Baseline", parasitaemias$Timepoint, fixed=T)
  
  parasitaemias$Timepoint <- gsub("D", "C", parasitaemias$Timepoint, fixed=T)
  
  #parasitaemias$Treatment <- factor(parasitaemias$Treatment, levels=c("before", "after"))
  # dods <- data.frame(vol=c('v02', 'v03', 'v05', 'v06', 'v07', 'v09'),
  #                    paras=c(4907, 8054, 16733, 7464, 21870, 15051),
  #                    dod=c('D15.5', 'D12.5', 'D15.5', 'D15.5', 'D16', 'D17'))
  # dods$inter <- match(dods$dod, unique(parasitaemias$Timepoint))
  
  parasitaemia_levels <- unique(gtools::mixedsort(parasitaemias$Timepoint))
  parasitaemias$Timepoint <- factor(parasitaemias$Timepoint,  levels=parasitaemia_levels)
  parasitaemias$Treatment <- factor(parasitaemias$Treatment, levels = c("before Treatment", "after Treatment"))
  
  before_treatment <- filter(parasitaemias, Treatment=="before Treatment")
  after_treatment <- filter(parasitaemias, Treatment=="after Treatment")
  
  para_linetype <- c("before Treatment"="solid", "after Treatment"="dotted")
  
  parasitaemia_curves <- ggplot()+
    geom_point(data=parasitaemias[!is.na(parasitaemias$Genomes),],
               aes(color=factor(Volunteer),
                   x=Timepoint,
                   y=Genomes+1,
               ),
               fill="white", stroke=1, size=0.9, shape=21)+
    geom_line(data=before_treatment[!is.na(before_treatment$Genomes),],
              aes(color=factor(Volunteer),
                  x=Timepoint,
                  y=Genomes+1,
                  group=factor(Volunteer),
                  linetype=Treatment
              ),
              size=0.9)+
    geom_line(data=after_treatment[!is.na(after_treatment$Genomes),],
              aes(color=factor(Volunteer),
                  x=Timepoint,
                  y=Genomes+1,
                  group=factor(Volunteer),
                  linetype=Treatment
              ),
              size=0.5)+
    geom_point(data=parasitaemias[!is.na(parasitaemias$Genomes),],
               aes(color=factor(Volunteer),
                   x=Timepoint,
                   y=Genomes+1,
               ),
               fill="white", stroke=1, size=0.9, shape=21)+
    scale_color_manual(values=volunteer_palette)+
    scale_x_discrete(breaks=parasitaemia_levels)+
    theme_minimal()+
    scale_linetype_manual(values = para_linetype, guide = NULL)+
    ggtitle("Parasitaemia")+
    xlab("Day of Infection")+
    geom_hline(yintercept = 20, linetype="dashed")+
    ylab("Genome Copies / mL")+
    scale_y_continuous(trans="log10", limits=c(1, 27000), breaks=c(10, 100, 1000, 10000))+
    fig1_theme+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=9),
          legend.position = "none",
          plot.title = element_text(size=10, hjust=0.5, vjust=0),
          axis.text.x = element_text(hjust=1, angle=45, size=5))
  
  
  
  #ggsave ("/Users/s1249052/PhD/oxford/vac69/parasitaemias_vac69.pdf", height = 8, width=10)
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/vac69a_parasitaemia.pdf", parasitaemia_curves, width=6, height=5)
  
  fig1 <- cowplot::plot_grid(parasitaemia_curves, all_ae_stack,
                             thrombos_lymphs, fever_curves, align = "h", axis = "tblr", rel_widths=c(2, 1.66), rel_heights = c(2,1.6))
  
  ggsave("~/PhD/cytof/vac69a/final_figures_for_paper/fig1_no_heatmap.pdf", fig1, height = 5, width=8)
  




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


plasma_levels <- scan("~/PhD/plasma/vac69a/analytes_sorted_by_padj_better.txt", what = ",")
plasma_levels <- append(plasma_levels[-20], "IL18", after = 7)


t_trans_data <- t(trans_data)
colnames(t_trans_data) <- t_trans_data[1,]

plasma_matrix <- t_trans_data[2:nrow(t_trans_data),]

class(plasma_matrix) <- "numeric"


plasma_matrix <- plasma_matrix[match(plasma_levels, rownames(plasma_matrix)),]


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
                                      Volunteer = list(title = "Volunteer", at = names(Volunteer), labels_gp = gpar(fontsize=8),title_gp = gpar(fontsize=10, fontface="bold"), legend_gp = gpar(fill = unname(Volunteer), fontsize=8), title_position = "topleft"),
                                      Timepoint = list(title ="Timepoint",at = names(Timepoint), labels_gp = gpar(fontsize=8), title_gp = gpar(fontsize=10, fontface="bold"), legend_gp = gpar(fill = unname(Timepoint), fontsize=8), title_position = "topleft")
                                    )
)


# left anno #
number_of_hits <- 12
number_of_maybes <- 0

plasma_matrix <- head(plasma_matrix, n=17)


significant <-  c(rep("<0.05",number_of_hits), rep("<0.1", number_of_maybes), rep(">0.05", nrow(plasma_matrix)-number_of_hits-number_of_maybes))
sig <- c("<0.05"="darkgreen",
         #"<0.1" = "lightgreen",
         ">0.05"= "lightgrey")


left_anno <-  rowAnnotation(gap = unit(5, "mm"),
                            #annotation_name_gp = gpar(angle=45),
                            show_annotation_name = FALSE,
                            "significant"=significant,
                            simple_anno_size = unit(2.5, "mm"), # width of the significance bar
                            col=list("significant" = c("<0.05"="darkgreen", ">0.05" ="lightgrey")),
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

rownames(plasma_matrix) <- gsub("DDimer", "D-Dimer", rownames(plasma_matrix))
rownames(plasma_matrix) <- gsub("IFNy", "IFNγ", rownames(plasma_matrix))





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
                     heatmap_legend_param = list(col = col_fun4, title = "Z-Score", title_position = "topleft", title_gp = gpar(fontsize=10, fontface="bold")),
                     #width = unit(16, "cm"),
                     #height = unit(16, "cm")
)



pdf("~/PhD/cytof/vac69a/final_figures_for_paper/plasma_zscore_heatmap.pdf", width=7, height=4.5)
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


data3[,3:ncol(data3)] <- log10(data3[,3:ncol(data3)])


data3 <- data3 %>%
  mutate(Volunteer = gsub("00", "0", Volunteer)) %>%
  mutate(timepoint = gsub("C-1", "Baseline", timepoint)) %>%
  mutate(timepoint = gsub("+", "", timepoint, fixed = T)) %>%
  mutate(timepoint = gsub("DoD", "Diagnosis", timepoint))

data3$timepoint <- factor(data3$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))

long_data3 <- tidyr::gather(data3, analyte, concentration, colnames(data3)[3:ncol(data3)])
long_data3[,1:3] <- lapply(long_data3[,1:3], as.character)

list_of_dfs_for_glm <- split(long_data3, long_data3$analyte)

#list_of_models <- lapply(list_of_dfs_for_glm, function(x) glm(x[,3]~x[,2]+x[,1], data=x))
list_of_models <- lapply(list_of_dfs_for_glm, function(x) nlme::lme(concentration~timepoint, random=~1|Volunteer, data=x))


#list_of_summaries <- lapply(list_of_models, function(x) cbind(summary(x)$coefficients, names(x$data)[3]))
list_of_summaries <- lapply(list_of_models, function(x)cbind(summary(x)$tTable, "analyte"=unique(x$data$analyte)))


df_of_model_results <- data.frame(do.call(rbind, list_of_summaries))
#colnames(df_of_model_results) <- c("Estimate", "SE", "t_value", "raw_p", "Analyte")
df_of_model_results$Coefficient <- rownames(df_of_model_results)


df_of_model_results <- df_of_model_results[!grepl("Intercept", df_of_model_results$Coefficient, fixed=TRUE),]
#df_of_model_results <- df_of_model_results[!grepl("v0", df_of_model_results$Coefficient),]

df_of_model_results$p_adj <- p.adjust(as.numeric(as.character(df_of_model_results$p.value)), method = "fdr")

sig_hits <- subset(df_of_model_results, df_of_model_results$p_adj<0.1)


sig_levels<- as.character(sig_hits[order(sig_hits$p_adj),]$analyte)



siggy_hits <- sig_hits %>%
  group_by(analyte) %>%
  top_n(n = -1, wt = p_adj)

# sig_levels <- as.character(siggy_hits[order(siggy_hits$p_adj),]$Analyte)
# # this next step is necessary because TGFbeta is undetectable/unchanged and top_n returns it three times because reasons
# sig_levels <- unique(sig_levels)
# write.table(sig_levels, "analytes_sorted_by_padj.txt", row.names = FALSE, col.names = "Analyte")

long_data <- gather(data3, Analyte, Concentration, colnames(data3)[3:ncol(data3)])

long_data$timepoint <- gsub("DoD", "Diagnosis", long_data$timepoint)

sig_glm_data <- subset(long_data, long_data$Analyte %in% siggy_hits$analyte)
sig_glm_data$AnalyteF <- factor(sig_glm_data$Analyte, levels=sig_levels)


sig_glm_data$p_adj <- sig_hits$p_adj[match(sig_glm_data$Analyte, sig_hits$Analyte)]

sig_glm_plot <- ggplot(sig_glm_data, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45")), y=2^Concentration, color=Volunteer, group=Volunteer))+
  geom_line(aes(color=Volunteer), size=1.1)+
  geom_point(fill="white", stroke=1, shape=21)+
  scale_y_continuous(labels=scales::comma)+
  facet_wrap(~ AnalyteF, scales = "free", ncol=6)+
  ylab("Plasma Concentration (pg / mL)")+
  theme_minimal()+
  guides(color=guide_legend(override.aes = list("size"=0.1)))+
  scale_color_manual(values=volunteer_palette)+
  fig1_theme+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust=1, angle=45, size=6),
        strip.background = element_rect(fill = "white", color = "white"))


ggsave(filename = "~/PhD/cytof/vac69a/final_figures_for_paper/supp_glm_sig_analytes_fdr_10-e-1.pdf", sig_glm_plot, width = 12, height=3.5)

