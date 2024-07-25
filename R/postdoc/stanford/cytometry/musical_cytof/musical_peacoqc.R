library(flowCore)
library(FlowSOM)
library(ggplot2) 


filter_plot <- function(ff_pre, ff_post, title, channel_x, channel_y){
  df <- data.frame(x = exprs(ff_pre)[,channel_x],
                   y = exprs(ff_pre)[,channel_y])
  i <- sample(nrow(df), 10000)
  if (!"Original_ID" %in% colnames(exprs(ff_pre))) {
    ff_pre@exprs <- cbind(ff_pre@exprs,
                          Original_ID = seq_len(nrow(ff_pre@exprs)))
  }
  p <- ggplot(df[i,], aes(x = x, y = y)) +
    geom_point(size = 0.5,
               color = ifelse(exprs(ff_pre)[i,"Original_ID"] %in%
                                exprs(ff_post)[,"Original_ID"], 'blue', 'red')) +
    xlab(GetMarkers(ff_pre, channel_x)) +
    ylab(GetMarkers(ff_pre, channel_y)) +
    theme_minimal()+
    theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}

dir_raw <- "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/"
dir_prepr <- "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/"

files <- list.files(dir_raw, pattern="*.fcs", full.names = FALSE, recursive = FALSE)
dir_QC <- "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/peacoQC"

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")
channels_of_interest <- musical_panel$fcs_colname[-c(1,2,5,40:42, 44:48)]

#CD19
reference_marker <- "Nd142Di"

# run peacoQC
for (file in files){ 
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE) 
  # ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest) 
  # ff_c <- flowCore::compensate(ff_m, compensation_matrix) 
  # ff_t <- flowCore::transform(ff, translist)
  # ff_s <- PeacoQC::RemoveDoublets(ff_t, channel1 = "Residual", channel2 = "Event_length") 
  #selected_live <- filter(ff_s, live_gate) 
  #ff_l <- ff_s[selected_live@subSet,] 
  PQC <- PeacoQC::PeacoQC(ff = ff, 
                          channels = channels_of_interest, 
                          plot = TRUE,
                          IT_limit=0.7,
                          MAD=6,
                          remove_zeros=TRUE,
                          time_units=500000,
                          save_fcs = TRUE,
                          suffix_fcs = "peacoQC",
                          output_directory = dir_QC) 
  # write.FCS(PQC$FinalFF, 
  #           file = paste0(dir_prepr, file)) 
  to_plot <- list(
                  # list(ff_pre = ff, 
                  #      ff_post = ff_m, 
                  #      title = "Removed margin events", 
                  #      channel_x = "AmCyan-A", 
                  #      channel_y = "APC-Cy7-A"), 
                  # list(ff_pre = ff_t, 
                  #      ff_post = ff_l, 
                  #      title = "Removed doublets", 
                  #      channel_x = "FSC-A", 
                  #      channel_y = "FSC-H"), 
                  list(ff_pre = ff, 
                       ff_post = PQC$FinalFF, 
                       title = "Removed low quality events", 
                       channel_x = "Time",
                       # CD3
                       channel_y = "Nd150Di"),
                  # list(ff_pre = ff_s,
                  #      ff_post = ff_l, 
                  #      title = "Removed debris and dead cells", 
                  #      channel_x = "FSC-A", 
                  #      channel_y = "APC-Cy7-A"), 
                  list(ff_pre = ff, 
                       ff_post = PQC$FinalFF, 
                       title = "Removed low quality events", 
                       channel_x = "Time", 
                       channel_y = "Ir191Di"))
  
  plot_list <- list()
  for (plot in to_plot) {
    plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre, 
                                                      ff_post = plot$ff_post, 
                                                      title = plot$title, 
                                                      channel_x = plot$channel_x, 
                                                      channel_y = plot$channel_y) } 
  png(paste0(dir_QC, sub("fcs", "png", file)), width = 1920)
  print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
  dev.off()
}
