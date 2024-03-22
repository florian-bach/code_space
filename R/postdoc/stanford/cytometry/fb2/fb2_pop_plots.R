library(ggplot2); theme_set(theme_minimal())
library(tidyr)
library(dplyr)


add_summary_rows <- function(.data, ...) {
  group_modify(.data, function(x, y) bind_rows(x, summarise(x, ...)))
}

id_colors <- c("darkblue", "orange")
comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")

data <- read.csv("~/postdoc/stanford/cytometry/attune/FB02/big_pop_table.csv")

long_data <- data %>%
  pivot_longer(cols=colnames(data)[2:(ncol(data)-1)], names_to="cell_pop", values_to="freq")%>%
  mutate("sample"=X)%>%
  select(sample, cell_pop, freq, -X.1)%>%
  mutate(cell_pop=gsub("Lymphocytes.Single.Cells.T.Cells.", "", .$cell_pop))%>%
  mutate("lineage"=case_when(
    grepl("Vd2*", .$cell_pop) ~ "Vd2",
    grepl("CD4.neg*", .$cell_pop) ~ "CD4neg",
    grepl("^CD4..*", .$cell_pop) ~ "CD4",
    grepl("Monocytes*", .$cell_pop) ~ "Monocyte",
    ))%>%
  mutate(lineage=if_else(is.na(lineage), "CD3", lineage))%>%
  mutate(cell_pop=gsub("..Freq..of.Parent....", "", .$cell_pop))%>%
  mutate("marker1"=case_when(
    grepl("CD69*", .$cell_pop) ~ "CD69",
    grepl("OX40*", .$cell_pop) ~ "OX40",
    grepl("CD137*", .$cell_pop) ~ "CD137",
    grepl("CD25*", .$cell_pop) ~ "CD25",
    grepl("ICOS*", .$cell_pop) ~ "ICOS",
    grepl("CD40L*", .$cell_pop) ~ "CD40L",
    grepl("CD200*", .$cell_pop) ~ "CD200",
    grepl("PDL1*", .$cell_pop) ~ "PDL1",
    # grepl("CD4-*", .$cell_pop) ~ "CD4",
  ))%>%
  mutate("marker2"=case_when(
    # grepl("CD4-*", .$cell_pop) ~ "CD4",
    grepl("PDL1*", .$cell_pop) ~ "PDL1",
    grepl("CD200*", .$cell_pop) ~ "CD200",
    grepl("CD40L*", .$cell_pop) ~ "CD40L",
    grepl("ICOS*", .$cell_pop) ~ "ICOS",
    grepl("CD25*", .$cell_pop) ~ "CD25",
    grepl("CD137*", .$cell_pop) ~ "CD137",
    grepl("OX40*", .$cell_pop) ~ "OX40",
    grepl("CD69*", .$cell_pop) ~ "CD69",
  ))%>%
  mutate("marker_combo"=paste(marker1, marker2))%>%
  # mutate(marker_combo = factor(marker_combo, levels=c("OX40 CD137", "CD69 PDL1", "ICOS CD40L", "CD25 CD200")))%>%
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
    grepl("Q21..", .$cell_pop, fixed = TRUE) ~ "Q1",
    grepl("Q22..", .$cell_pop, fixed = TRUE) ~ "Q2",
    grepl("Q23..", .$cell_pop, fixed = TRUE) ~ "Q3",
    grepl("Q24..", .$cell_pop, fixed = TRUE) ~ "Q4",
  ))%>%
  mutate("stim"=case_when(
    grepl("3.fcs$", .$sample) ~ "iRBC",
    grepl("7.fcs$", .$sample) ~ "iRBC",
    grepl("4.fcs$", .$sample) ~ "media",
    grepl("8.fcs$", .$sample) ~ "media",
    grepl("5.fcs$", .$sample) ~ "PMA",
    grepl("9.fcs$", .$sample) ~ "PMA",
  ))%>%
  mutate("serum"=case_when(
    grepl("3.fcs$", .$sample) ~ "FBS",
    grepl("4.fcs$", .$sample) ~ "FBS",
    grepl("5.fcs$", .$sample) ~ "FBS",
    grepl("7.fcs$", .$sample) ~ "autologous",
    grepl("8.fcs$", .$sample) ~ "autologous",
    grepl("9.fcs$", .$sample) ~ "autologous",
  ))%>%
  mutate("id"=case_when(
    grepl("^B", .$sample) ~ "2003",
    grepl("^C", .$sample) ~ "2004",
    grepl("^D", .$sample) ~ "2005",
    grepl("^E", .$sample) ~ "2008",
    grepl("^F", .$sample) ~ "2009",
    grepl("^G", .$sample) ~ "2012",
  ))%>%
  filter(quadrant!="Q4")%>%
  mutate("origin"="Uganda")%>%
  group_by(sample, lineage, marker_combo, stim,  serum, id)%>%
  add_summary_rows(quadrant="sum(Q)",
                   freq=sum(freq))
  
# write fcs file metadata
# meta_data <- long_data %>%
#   ungroup()%>%
#   select(sample, id, serum, stim)%>%
#   filter(!duplicated(sample))
#   
# write.csv(meta_data, "~/postdoc/stanford/cytometry/attune/FB02/sample_metadata.csv")
write.csv(long_data, "~/postdoc/stanford/cytometry/attune/FB02/fb02_long_data.csv")


long_data_medians <- long_data %>%
  filter(!is.na(stim))%>%
  group_by(lineage, stim, serum, cell_pop, marker_combo, quadrant)%>%
  summarise("median"=round(median(freq), digits = 2))


#Q1 is single-positive for second marker (Y axis)
#Q2 is double positive
#Q3 is single positive for first marker (X axis)
# this plot allows the comparison of plasma conditions within stim conditions
(cd4_cd4neg_summary_plot <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         marker_combo != "NA NA",
         # lineage %in% c("CD4", "CD4neg"),
         quadrant != "Q4")%>%
  ggplot(aes(x=quadrant, y=(freq+0.0001)/100))+
  geom_point(aes(color=id, group=serum), alpha=0.6, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=serum), outlier.shape = NA)+
  facet_grid(lineage+marker_combo~stim)+
  geom_text(aes(x=quadrant, y=1, label=median, group=serum), position = position_dodge(width=0.75), data = long_data_medians, size=2.8)+
  scale_y_log10(labels = scales::label_percent())+
  theme_minimal()+
  scale_fill_manual(values=c("#fec200", "#0047ab"))+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_blank()))

ggsave("~/postdoc/stanford/cytometry/attune/FB2/figures/cd4_cd4neg_summary_plot.png", cd4_cd4neg_summary_plot, width=8*1.5, height=8*1.5, bg="white")





# this plot allows the comparison of stim conditions within plasma conditions
cd4_cd4neg_summary_plot2 <- long_data %>%
  filter(!is.na(marker_combo),
         !is.na(lineage),
         !is.na(stim),
         marker_combo != "NA NA",
         # lineage %in% c("CD4", "CD4neg"),
         quadrant != "Q4")%>%
  ggplot(aes(x=quadrant, y=(freq+0.0001)/100))+
  geom_point(aes(color=id, group=stim), alpha=0.6, position = position_dodge(width=0.75))+
  geom_boxplot(aes(fill=stim), outlier.shape = NA)+
  geom_text(aes(x=quadrant, y=1, label=median, group=stim), position = position_dodge(width=0.75), data = long_data_medians, size=2.8)+
  facet_grid(lineage+marker_combo~serum)+
  scale_y_log10(labels = scales::label_percent())+
  theme_minimal()+
  scale_fill_manual(values=c("darkred", "darkgrey", "purple"))+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB2/figures/cd4_cd4neg_summary_plot2.png", cd4_cd4neg_summary_plot2, width=8*1.5, height=8*1.5, bg="white")


(cd4_cd4neg_summary_plot3 <- long_data %>%
    filter(marker_combo %in% c("CD69 CD40L", "OX40 CD137"),
           !is.na(lineage),
           !is.na(stim),
           # stim %in% c("media", "PMA"),
           marker_combo != "NA NA",
           lineage %in% c("CD4", "CD4neg"),
           quadrant != "Q4")%>%
    ggplot(aes(x=quadrant, y=(freq+0.0001)/100))+
    geom_point(aes(color=id, group=serum), alpha=0.6, position = position_dodge(width=0.75))+
    geom_boxplot(aes(fill=serum), outlier.shape = NA)+
    facet_wrap(lineage+marker_combo~stim, scales="free", nrow=2)+
    scale_y_continuous(labels = scales::label_percent())+
    theme_minimal()+
    scale_fill_manual(values=c("#fec200", "#0047ab"))+
    theme(axis.text.x = element_text(angle=90, hjust=1),
          axis.title = element_blank()))

ggsave("~/postdoc/stanford/cytometry/attune/FB02/figures/cd4_cd4neg_summary_plot3.png", cd4_cd4neg_summary_plot3, width=8*1.5, height=6*1.5, bg="white")

# comparison of cell counts in iRBC stim condition ####

cell_counts <- readxl::read_xlsx("~/postdoc/stanford/cytometry/attune/FB2/iRBC_stim_cell_counts.xlsx")

iRBC_stim_cell_count_plot <- cell_counts %>%
  pivot_longer(cols=c(t_cells, non_t_lymphs), names_to = "cell_pop", values_to="count")%>%
  mutate("id"=case_when(
    grepl("^B", .$sample) ~ "2003",
    grepl("^C", .$sample) ~ "2004",
    grepl("^D", .$sample) ~ "2005",
    grepl("^E", .$sample) ~ "2008",
    grepl("^F", .$sample) ~ "2009",
    grepl("^G", .$sample) ~ "2012",
  ))%>%
  mutate(cell_pop=if_else(cell_pop=="t_cells", "T cells", "non T cells"))%>%
  ggplot(aes(x=cell_pop, y=count))+
  geom_boxplot(aes(fill=condition), outlier.shape = NA)+
  geom_point(aes(color=id, group=condition), position=position_dodge(width=0.75), alpha=0.6)+
  # geom_line(aes(color=id, group=cell_pop), position=position_dodge(width=0.75))+
  # scale_y_log10()+
  scale_fill_manual(values=c("#fec200", "#0047ab"))+
  ylab("cell count")+
  # facet_grid(lineage+marker_combo~serum)+
  theme_minimal()+
  theme(axis.title.x = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB02/figures/iRBC_stim_cell_counts.png", iRBC_stim_cell_count_plot, width=3, height=3, bg="white")


aim_cell_counts <- read.csv("~/postdoc/stanford/cytometry/attune/FB02/counts.csv")

iRBC_stim_aim_cell_counts <- aim_cell_counts %>%
  select(-T.Cells, -Non.T.cells)%>%
  pivot_longer(cols=c(CD4.neg.Q1., CD4.neg.Q2, CD4.neg.Q3, CD4.Q1, CD4.Q2, CD4.Q3), names_to = "cell_pop", values_to="count")%>%
  mutate("id"=case_when(
    grepl("^B", .$File) ~ "2003",
    grepl("^C", .$File) ~ "2004",
    grepl("^D", .$File) ~ "2005",
    grepl("^E", .$File) ~ "2008",
    grepl("^F", .$File) ~ "2009",
    grepl("^G", .$File) ~ "2012",
  ))%>%
  mutate("condition"=case_when(
    grepl("3.fcs$", .$File) ~ "FBS",
    grepl("4.fcs$", .$File) ~ "FBS",
    grepl("5.fcs$", .$File) ~ "FBS",
    grepl("7.fcs$", .$File) ~ "autologous",
    grepl("8.fcs$", .$File) ~ "autologous",
    grepl("9.fcs$", .$File) ~ "autologous",
  ))%>%
  mutate("stim"=case_when(
    grepl("3.fcs$", .$File) ~ "iRBC",
    grepl("7.fcs$", .$File) ~ "iRBC",
    grepl("4.fcs$", .$File) ~ "media",
    grepl("8.fcs$", .$File) ~ "media",
    grepl("5.fcs$", .$File) ~ "PMA",
    grepl("9.fcs$", .$File) ~ "PMA",
  ))%>%
  mutate("lineage"=case_when(
    grepl("^CD4.Q", .$cell_pop) ~ "CD4+",
    grepl("^CD4.n", .$cell_pop) ~ "CD4neg",
  ))%>%
  group_by(lineage, stim,  condition, id)%>%
  add_summary_rows(cell_pop="quadrant sum",
                   count=sum(count))


  # filter(stim=="iRBC")%>%
iRBC_stim_aim_cell_counts_plot <- ggplot(iRBC_stim_aim_cell_counts, aes(x=cell_pop, y=count+0.1))+
  geom_boxplot(aes(fill=condition), outlier.shape = NA)+
  geom_point(aes(color=id, group=condition), position=position_dodge(width=0.75), alpha=0.6)+
  # geom_line(aes(color=id, group=cell_pop), position=position_dodge(width=0.75))+
  scale_y_log10(labels=scales::label_log())+
  scale_fill_manual(values=c("#fec200", "#0047ab"))+
  ylab("cell count")+
  facet_wrap(lineage~stim, nrow = 2, scales="free_x")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1))

ggsave("~/postdoc/stanford/cytometry/attune/FB02/figures/iRBC_stim_aim_cell_count_plot.png", iRBC_stim_aim_cell_counts_plot, width=8, height=6, bg="white")


# # parasitaemia plots
# 
# mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/MICDROP all visit database through July 31st 2023.dta")
# 
# # merge parasitemia data so that qPCR takes precedent when both slide and qPCR are present
# mic_drop <- mic_drop %>%
#   mutate(mstatus = case_match(mstatus,
#                               0~"no malaria",
#                               1~"uncomplicated",
#                               2~"complicated",
#                               3~"quinine for AL failure",
#                               4~"Q/AS failure"))%>%
#   mutate(visit_id = paste(id, date, sep=""))%>%
#   mutate("parasitaemia_method" = if_else(qPCRdich==1, "qPCR", if_else(BSdich==1, "smear", "dunno")))%>%
#   mutate(any_parsdens = if_else(is.na(qPCRparsdens) & !is.na(pardens), pardens, qPCRparsdens))%>%
#   mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))
# 
# parasitaemia_plot <- mic_drop %>%
#   filter(id %in% c(10460, 10958),
#          !is.na(any_parsdens),
#          # ageinwks < 25
#   )%>%
#   ggplot(., aes(x=ageinwks, y=as.numeric(any_parsdens)+0.001))+
#   geom_point(aes(color=factor(mstatus, #levels=c("0",
#                               #       "1",
#                               #      "2",
#                               #     "3")
#   )))+
#   geom_line(alpha=0.3, aes(group=id))+
#   facet_wrap(~ id)+
#   ylab("qPCR parasites / μl\n")+
#   xlab("age in weeks")+
#   scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
#   theme_minimal()+
#   geom_vline(xintercept = 24, linetype="dashed")+
#   scale_color_manual(values=comp_pal)+
#   # scale_shape_manual(values=c(16,15))+
#   guides(color=guide_legend(title=""))+
#   theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))
# 
# ggsave("~/postdoc/stanford/cytometry/attune/FB1/figures/parasitaemia_plot.png", parasitaemia_plot, width=8, height=4, bg="white")

