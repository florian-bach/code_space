marker_info <- data.frame(
  channel_name = paste0("channel", sprintf("%03d", 1:20)), 
  marker_name = paste0("marker", sprintf("%02d", 1:20)), 
  marker_class = factor(c(rep("type", 10), rep("state", 10)), 
                        levels = c("type", "state", "none")), 
  stringsAsFactors = FALSE
)


d_se <- prepareData(vac69a, ei, panel)

# Generate clusters
d_se <- generateClusters(d_se)

#calculate medians
d_medians <- calcMedians(d_se)

#extract state names
state_names <- names(assays(d_medians))[markers_to_test]

metadata(d_medians)$id_state_markers
