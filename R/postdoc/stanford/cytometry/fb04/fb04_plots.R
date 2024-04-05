library(ggplot2); theme_set(theme_minimal())
library(tidyr)
library(dplyr)


add_summary_rows <- function(.data, ...) {
  group_modify(.data, function(x, y) bind_rows(x, summarise(x, ...)))
}

id_colors <- c("darkblue", "orange")
comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")

serum_pal <- c("no serum"="lightgrey",
               "human AB"="#fec200",
               "autologous"="#D67500",
               "UG"="#D67500",
               "UG-T"="#D64500",
               "FBS"="#5A270F"
)

# activation percentages ####
flo_plate <- read.csv("~/postdoc/stanford/cytometry/attune/FB04/pop_freqs.csv")

flo_long_data <- {flo_plate %>%
    pivot_longer(cols=colnames(flo_plate)[2:(ncol(flo_plate)-1)], names_to="cell_pop", values_to="freq")%>%
    mutate("sample"=X)%>%
    select(sample, cell_pop, freq, -X.1)%>%
    mutate(cell_pop=gsub("WBC.Single.Cells.T.Cells.", "", .$cell_pop))%>%
    mutate("lineage"=case_when(
      grepl("Vd2*", .$cell_pop) ~ "Vd2",
      grepl("CD4neg*", .$cell_pop) ~ "CD4neg",
      grepl("^CD4.*", .$cell_pop) ~ "CD4",
      grepl("Monocytes*", .$cell_pop) ~ "Monocyte",
    ))%>%
    mutate(lineage=if_else(is.na(lineage), "CD3", lineage))%>%
    mutate(cell_pop=gsub("..Freq..of.Parent....", "", .$cell_pop))%>%
    mutate("marker1"=case_when(
      grepl("CD69*", .$cell_pop) ~ "CD69",
      grepl("OX40*", .$cell_pop) ~ "OX40",
      grepl("CD137*", .$cell_pop) ~ "CD137",
      grepl("CD40L*", .$cell_pop) ~ "CD40L",
      # grepl("CD4-*", .$cell_pop) ~ "CD4",
    ))%>%
    mutate("marker2"=case_when(
      # grepl("CD4-*", .$cell_pop) ~ "CD4",
      grepl("CD40L*", .$cell_pop) ~ "CD40L",
      grepl("CD137*", .$cell_pop) ~ "CD137",
      grepl("OX40*", .$cell_pop) ~ "OX40",
      grepl("CD69*", .$cell_pop) ~ "CD69",
    ))%>%
    mutate("marker_combo"=paste(marker1, marker2))%>%
    mutate(marker_combo = factor(marker_combo, levels=c("OX40 CD137", "CD69 CD40L")))%>%
    mutate("quadrant"=case_when(
      grepl("Q1..", .$cell_pop, fixed = TRUE) ~ "Q1",
      grepl("Q2..", .$cell_pop, fixed = TRUE) ~ "Q2",
      grepl("Q3..", .$cell_pop, fixed = TRUE) ~ "Q3",
      grepl("Q4..", .$cell_pop, fixed = TRUE) ~ "Q4",
      grepl("Q5..", .$cell_pop, fixed = TRUE) ~ "Q1",
      grepl("Q6..", .$cell_pop, fixed = TRUE) ~ "Q2",
      grepl("Q7..", .$cell_pop, fixed = TRUE) ~ "Q3",
      grepl("Q8..", .$cell_pop, fixed = TRUE) ~ "Q4",
      grepl("Q9..", .$cell_pop, fixed = TRUE) ~ "Q1",
      grepl("Q10..", .$cell_pop, fixed = TRUE) ~ "Q2",
      grepl("Q11..", .$cell_pop, fixed = TRUE) ~ "Q3",
      grepl("Q12..", .$cell_pop, fixed = TRUE) ~ "Q4",
      grepl("Q13..", .$cell_pop, fixed = TRUE) ~ "Q1",
      grepl("Q14..", .$cell_pop, fixed = TRUE) ~ "Q2",
      grepl("Q15..", .$cell_pop, fixed = TRUE) ~ "Q3",
      grepl("Q16..", .$cell_pop, fixed = TRUE) ~ "Q4",
      grepl("Q17..", .$cell_pop, fixed = TRUE) ~ "Q1",
      grepl("Q18..", .$cell_pop, fixed = TRUE) ~ "Q2",
      grepl("Q19..", .$cell_pop, fixed = TRUE) ~ "Q3",
      grepl("Q20..", .$cell_pop, fixed = TRUE) ~ "Q4",
    ))%>%
    mutate("stim"=case_when(
      grepl("3.fcs$", .$sample) ~ "media",
      grepl("6.fcs$", .$sample) ~ "media",
      grepl("9.fcs$", .$sample) ~ "media",
      
      grepl("4.fcs$", .$sample) ~ "PMA",
      grepl("7.fcs$", .$sample) ~ "PMA",
      grepl("10.fcs$", .$sample) ~ "PMA",
      
      grepl("5.fcs$", .$sample) ~ "blank",
      grepl("8.fcs$", .$sample) ~ "blank",
      grepl("11.fcs$", .$sample) ~ "blank",
      grepl("12.fcs$", .$sample) ~ "blank",
    ))%>%
    mutate("serum"=case_when(
      grepl("3.fcs$", .$sample) ~ "human AB",
      grepl("4.fcs$", .$sample) ~ "human AB",
      grepl("6.fcs$", .$sample) ~ "UG",
      grepl("7.fcs$", .$sample) ~ "UG",
      grepl("9.fcs$", .$sample) ~ "UG-T",
      grepl("10.fcs$", .$sample) ~ "UG-T",
    ))%>%
    mutate("id"=case_when(
      grepl("^A", .$sample) ~ "400",
      grepl("^B", .$sample) ~ "410",
      grepl("^C", .$sample) ~ "415",
      grepl("^D", .$sample) ~ "420",
      grepl("^E", .$sample) ~ "425",
      grepl("^F", .$sample) ~ "430",
      grepl("^G", .$sample) ~ "440",
      grepl("^H", .$sample) ~ "445"
    ))%>%
    filter(quadrant!="Q4")%>%
    # mutate("person"="Flo")%>%
    group_by(sample, lineage, marker_combo, stim,  serum, id)%>%
    add_summary_rows(quadrant="sum(Q)",
                     freq=sum(freq))}



long_data <- flo_long_data
# write.csv(long_data, "~/postdoc/stanford/cytometry/attune/FB03/fb03_long_data.csv")


long_data_medians <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         stim!="blank",
         marker_combo != "NA NA",
         lineage %in% c("CD4", "CD4neg"),
         quadrant != "Q4")%>%
  group_by(lineage, stim, serum, cell_pop, marker_combo, quadrant)%>%
  summarise("median"=round(median(freq), digits = 2))


#Q1 is single-positive for second marker (Y axis)
#Q2 is double positive
#Q3 is single positive for first marker (X axis)
# this plot allows the comparison of plasma conditions within stim conditions

cd4_cd4neg_summary_data <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         stim!="blank",
         marker_combo != "NA NA",
         # lineage %in% c("CD4", "CD4neg"),
         quadrant != "Q4")

  ##cd4 cd4neg plot ####

(cd4_cd4neg_summary_plot <- ggplot(cd4_cd4neg_summary_data, aes(x=quadrant, y=(freq+0.0001)/100))+
    geom_point(aes(color=id, group=serum), alpha=0.6, position = position_dodge(width=0.75))+
    geom_boxplot(aes(fill=serum), outlier.shape = NA)+
    facet_grid(lineage+marker_combo~stim, scales = "free_y")+
    geom_text(aes(x=quadrant, y=1, label=median, group=serum), position = position_dodge(width=0.75), data = long_data_medians, size=2.8)+
    scale_y_log10(labels = scales::label_percent())+
    theme_minimal()+
    scale_fill_manual(values=serum_pal)+
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_blank()))

ggsave("~/postdoc/stanford/cytometry/attune/FB04/figures/cd4_cd4neg_summary_plot.png", cd4_cd4neg_summary_plot, width=8*1.5, height=10*1.5, bg="white")

##gamma delta cd4neg plot ####
gd_cd4neg_plot <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         stim!="blank",
         marker_combo == "OX40 CD137",
         lineage %in% c("Vd2", "CD4neg"),
         quadrant != "Q4")%>%
  ggplot(., aes(x=quadrant, y=(freq+0.0001)/100))+
  geom_point(aes(color=id, group=serum), alpha=0.6, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=serum), outlier.shape = NA)+
  facet_wrap(lineage+marker_combo~stim, scales = "free_y")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  scale_fill_manual(values=serum_pal)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB04/figures/gd_cd4neg_plot.png", gd_cd4neg_plot, width=8, height=4, bg="white")

  ## cd4 plot ####
cd4_plot <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         stim!="blank",
         marker_combo == "CD69 CD40L",
         lineage %in% c("CD4"),
         quadrant != "Q4")%>%
  ggplot(., aes(x=quadrant, y=(freq+0.0001)/100))+
  geom_point(aes(color=id, group=serum), alpha=0.6, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=serum), outlier.shape = NA)+
  facet_grid(lineage+marker_combo~stim, scales = "free_y")+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  scale_fill_manual(values=serum_pal)+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB04/figures/cd4_plot.png", cd4_plot, width=6, height=4, bg="white")


# cell counts ####
flo_counts <- read.csv("~/postdoc/stanford/cytometry/attune/FB03/2024_FB03_fbs_ab_no_serum 2/flo/flo_cell_counts.csv")
julio_counts <- read.csv("~/postdoc/stanford/cytometry/attune/FB03/2024_FB03_fbs_ab_no_serum 2/julio/julio_cell_counts.csv")

flo_count_data<- flo_counts %>%
  mutate("sample"=X)%>%
  pivot_longer(cols=colnames(.)[2:(ncol(.)-1)], names_to="cell_pop", values_to="freq")%>%
  select(sample, cell_pop, freq, -X)%>%
  mutate(cell_pop=gsub("WBC.Singlets.", "", .$cell_pop))%>%
  mutate("lineage"=case_when(
    grepl("CD3....Count*", .$cell_pop) ~ "T cells",
    grepl("CD3neg...Count*", .$cell_pop) ~ "non T cells",
  ))%>%
  # mutate("marker1"=case_when(
  mutate("stim"=case_when(
    grepl("2.fcs$", .$sample) ~ "media",
    grepl("5.fcs$", .$sample) ~ "media",
    grepl("8.fcs$", .$sample) ~ "media",
    grepl("3.fcs$", .$sample) ~ "PMA",
    grepl("6.fcs$", .$sample) ~ "PMA",
    grepl("9.fcs$", .$sample) ~ "PMA",
    grepl("4.fcs$", .$sample) ~ "blank",
    grepl("7.fcs$", .$sample) ~ "blank",
    grepl("10.fcs$", .$sample) ~ "blank",
  ))%>%
  mutate("serum"=case_when(
    grepl("2.fcs$", .$sample) ~ "FBS",
    grepl("3.fcs$", .$sample) ~ "FBS",
    grepl("5.fcs$", .$sample) ~ "human AB",
    grepl("6.fcs$", .$sample) ~ "human AB",
    grepl("8.fcs$", .$sample) ~ "no serum",
    grepl("9.fcs$", .$sample) ~ "no serum",
  ))%>%
  mutate("id"=case_when(
    grepl("^B", .$sample) ~ "435",
    grepl("^C", .$sample) ~ "440",
    grepl("^D", .$sample) ~ "445",
    grepl("^E", .$sample) ~ "455",
    grepl("^F", .$sample) ~ "460",
    grepl("^G", .$sample) ~ "470",
  ))%>%
  mutate("person"="Flo")

julio_count_data <- julio_counts %>%
  mutate("sample"=X)%>%
  pivot_longer(cols=colnames(.)[2:(ncol(.)-1)], names_to="cell_pop", values_to="freq")%>%
  select(sample, cell_pop, freq, -X)%>%
  mutate(cell_pop=gsub("WBC.Singlets.", "", .$cell_pop))%>%
  mutate("lineage"=case_when(
    grepl("CD3....Count*", .$cell_pop) ~ "T cells",
    grepl("CD3neg...Count*", .$cell_pop) ~ "non T cells",
  ))%>%
  mutate("stim"=case_when(
    grepl("2.fcs$", .$sample) ~ "media",
    grepl("5.fcs$", .$sample) ~ "media",
    grepl("8.fcs$", .$sample) ~ "media",
    grepl("3.fcs$", .$sample) ~ "PMA",
    grepl("6.fcs$", .$sample) ~ "PMA",
    grepl("9.fcs$", .$sample) ~ "PMA",
    grepl("4.fcs$", .$sample) ~ "blank",
    grepl("7.fcs$", .$sample) ~ "blank",
    grepl("10.fcs$", .$sample) ~ "blank",
  ))%>%
  mutate("serum"=case_when(
    grepl("2.fcs$", .$sample) ~ "FBS",
    grepl("3.fcs$", .$sample) ~ "FBS",
    grepl("5.fcs$", .$sample) ~ "human AB",
    grepl("6.fcs$", .$sample) ~ "human AB",
    grepl("8.fcs$", .$sample) ~ "no serum",
    grepl("9.fcs$", .$sample) ~ "no serum",
  ))%>%
  mutate("id"=case_when(
    grepl("^B", .$sample) ~ "400",
    grepl("^D", .$sample) ~ "405",
    grepl("^F", .$sample) ~ "430",
  ))%>%
  mutate("person"="Julio")


count_data <- rbind(julio_count_data, flo_count_data)

stim_cell_count_plot <- count_data %>%
  filter(lineage %in% c("non T cells", "T cells"),
         !is.na(serum),
  )%>%
  ggplot(., aes(x=lineage, y=freq))+
  geom_boxplot(aes(fill=serum), outlier.shape = NA)+
  geom_point(aes(color=id, group=serum, shape=person), position=position_dodge(width=0.75), alpha=0.6)+
  # geom_line(aes(color=id, group=cell_pop), position=position_dodge(width=0.75))+
  # scale_y_log10()+
  scale_fill_manual(values=serum_pal)+
  ylab("cell count")+
  facet_grid(~stim)+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB03/figures/stim_cell_counts.png", stim_cell_count_plot, width=7, height=5, bg="white")


