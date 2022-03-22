# first import some stuff that so that we don't have to keep defining everything, even if we want to just run a chunk
stacked_bar_data <-  read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_all_cluster_freqs.csv")
stacked_bar_data$volunteer <- factor(stacked_bar_data$volunteer, levels=c("v313", "v315", "v320", "v306", "v301", "v308", "v305", "v304", "v310"))
stacked_bar_data$cluster_id <- gsub("activated CD8lo Effector", "activated DN", stacked_bar_data$cluster_id)
stacked_bar_data$timepoint <- gsub("DoD", "Diagnosis", stacked_bar_data$timepoint)
stacked_bar_data$timepoint <- factor(stacked_bar_data$timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45"))

# add column for lineage, relying on cluster_id names
stacked_bar_data$lineage <- substr(stacked_bar_data$cluster_id, nchar(stacked_bar_data$cluster_id)-2, nchar(stacked_bar_data$cluster_id))
#[1] "CD4" " gd" "AIT" "NKT" "CD8" " DN" "reg"

lin_replacement <- setNames(c("CD4", "gd", "MAIT", "NKT", "CD8", "DN", "Treg", "gd"), unique(stacked_bar_data$lineage))
stacked_bar_data$lineage <- stringr::str_replace_all(stacked_bar_data$lineage, lin_replacement)
stacked_bar_data$lineage <- factor(stacked_bar_data$lineage, levels=rev(c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "NKT")))





stacked_bar_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  mutate("scaled_freq"=frequency/sum)%>%
  ungroup()%>%
  filter(volunteer %in% c("v313", "v315", "v320"), cluster_id =="Naive CD4", timepoint %in% c("Baseline", "T6")) %>%
  summarise(mean(scaled_freq))


stacked_bar_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  mutate("scaled_freq"=frequency/sum)%>%
  ungroup()%>%
  filter(volunteer %in% c("v313", "v315", "v320"), lineage=="CD4", timepoint %in% c("Baseline", "T6")) %>%
  filter(grepl("*CM *", cluster_id, fixed=FALSE))%>%
  group_by(volunteer, timepoint)%>%
  summarise(sum(scaled_freq))



df <- stacked_bar_data %>%
  group_by(lineage, volunteer, timepoint) %>%
  mutate("sum"=sum(frequency)) %>%
  ungroup()%>%
  select(lineage, timepoint, volunteer, sum)%>%
  filter(lineage %in% c("CD4", "CD8", "gd" )) %>%
  summarise(lineage, timepoint, volunteer, sum)


df <- df[!duplicated(df$sum),]
