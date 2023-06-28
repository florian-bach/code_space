library(flowCore)
library(FlowSOM)
library(ggplot2) 

dir_prepr <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/Preprocessed/"
dir_raw <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/"
files <- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/", pattern="*.fcs", full.names = FALSE, recursive = FALSE)
dir_QC <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/peacoQC"

tfh_markers <- c("CD3", "CD4", "CD45RA", "CXCR5", "PD1", "CXCR3", "CCR4", "CCR6")
reference_marker <- "PE-A"


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


premessa_table <- premessa::read_parameters(list.files(dir_raw, pattern="*.fcs", full.names = TRUE, recursive = FALSE))
panel <- data.frame("fcs_colname"=rownames(premessa_table[1]),
                    "marker_name"=premessa_table[[1]],
                    "marker_class"="type")
panel$antigen <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)


channels_of_interest <-panel$fcs_colname[panel$antigen%in%tfh_markers]

#GetChannels(object = tfhs[[1]], markers = markers_of_interest, exact = FALSE)
compensation_matrix <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/comp_matrix.csv", check.names = FALSE, row.names = 1)


reference_file <- read.FCS(paste(dir_raw, files[1], sep=""))

ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix))
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
translist <- c(translist, 
               transformList("SSC-A", flowCore::linearTransform(a = SSCA_a, 
                                                                b = SSCA_b))) 






# run peacoQC
for (file in files){ 
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE, ) 
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest) 
  ff_c <- flowCore::compensate(ff_m, compensation_matrix) 
  ff_t <- flowCore::transform(ff_c, translist) 
  ff_l <- PeacoQC::RemoveDoublets(ff_t) 
  #selected_live <- filter(ff_s, live_gate) 
  #ff_l <- ff_s[selected_live@subSet,] 
  PQC <- PeacoQC::PeacoQC(ff = ff_l, 
                          channels = channels_of_interest, 
                          plot = TRUE, save_fcs = FALSE, 
                          output_directory = dir_QC) 
  write.FCS(PQC$FinalFF, 
            file = paste0(dir_prepr, file)) 
  to_plot <- list(list(ff_pre = ff, 
                       ff_post = ff_m, 
                       title = "Removed margin events", 
                       channel_x = "AmCyan-A", 
                       channel_y = "APC-Cy7-A"), 
                  list(ff_pre = ff_t, 
                       ff_post = ff_l, 
                       title = "Removed doublets", 
                       channel_x = "FSC-A", 
                       channel_y = "FSC-H"), 
                  list(ff_pre = ff_l, 
                       ff_post = PQC$FinalFF, 
                       title = "Removed low quality events", 
                       channel_x = "Time", 
                       channel_y = "Pacific Blue-A"),
                  # list(ff_pre = ff_s,
                  #      ff_post = ff_l, 
                  #      title = "Removed debris and dead cells", 
                  #      channel_x = "FSC-A", 
                  #      channel_y = "APC-Cy7-A"), 
                  list(ff_pre = ff_l, 
                       ff_post = PQC$FinalFF, 
                       title = "Removed low quality events", 
                       channel_x = "Time", 
                       channel_y = "APC-A"))
  
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