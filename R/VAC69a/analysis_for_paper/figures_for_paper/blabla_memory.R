
all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)
all_t6_data <- subset(all_t6_data, all_t6_data$timepoint == "Baseline")
cd4_t6_data <- subset(all_t6_data, grepl("CD4 ", all_t6_data$cluster_id))
cd4_memory_t6_data <- subset(cd4_t6_data, !grepl("Naïve", cd4_t6_data$cluster_id))


cd4_memory_t6_summary <- cd4_memory_t6_data %>%
  group_by(volunteer) %>%
  mutate("CD4_memory_percentage"=sum(frequency)) %>%
  select(volunteer, CD4_memory_percentage)

cd4_memory_t6_summary <- cd4_memory_t6_summary[!duplicated(cd4_memory_t6_summary), ]

cd4_bar_data <- cd4_t6_data %>%
  group_by(volunteer) %>%
  summarise("cd4_freq" = frequency/cd4_memory_t6_summary$CD4_memory_percentage[match(volunteer, cd4_memory_t6_summary$volunteer)]) %>%
  ungroup()


features <- c("CD4 Memory % Baseline", "activated CD27- cytotoxic CD4 EM")
try <- subset(big_fc_table, rownames(big_fc_table) %in% features)

ggplot()+
  geom_point(aes(x=unname(try[1,]), y=unname(try[2,]), color=names(try)))+
  xlab(features[2])+
  ylab(features[1])+
  theme_minimal()+
  scale_color_manual(values = volunteer_palette)
