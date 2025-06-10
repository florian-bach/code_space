library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(emmeans)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/tr1_paper/revised_baseline_clean_musical_combo_with_metadata.csv")%>%
  mutate(timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))

wide_data <- clean_data%>%
  filter(infectiontype %in% c("A", "A2", "S", "S2"), timepoint=="baseline")%>%
  
  pivot_wider(names_from = c(targetName), values_from = concentration, id_cols = c("id", "bino_infectiontype"))

write.csv(wide_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/wide_data.csv", row.names=FALSE)

wide_data%>%
  group_by(id)%>%
  filter(n()>1)%>%
  distinct(id)

clean_data%>%
  filter(infectiontype %in% c("A", "A2", "S", "S2"),
timepoint=="baseline",
#targetName %in% c("PTX3", "TNF", "FLT4", "MMP9")
targetName %in% c("KNG1", "TAFA5", "BMP7", "AGRP", "NGF", "CCL11"))%>%
  mutate("bino_infectiontype"=ifelse(grepl("^A", infectiontype), "A", "S"))%>%
  ggplot(., aes(x=bino_infectiontype, y=concentration))+
  geom_line(aes(group=id))+
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")
