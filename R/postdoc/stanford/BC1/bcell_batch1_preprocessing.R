library(flowStats)
library(flowCore)
library(FlowSOM)

# something is wrong with the fcs headers, the channel descriptions are messed up, the original id
# one is causing problems for some files in the peacoqc function

dir_prepr <- "~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/preprocessed/"
dir_raw <- "~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/"
files <- list.files(dir_raw, pattern="*.fcs", full.names = FALSE, recursive = FALSE)
dir_QC <- "~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/peacoQC/"
comp_dir <- "~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/comp_matrix.csv"

bcell_markers <- c("CD10", "CD19", "CD20", "CD21", "CD27", "CD38")
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





#GetChannels(object = tfhs[[1]], markers = markers_of_interest, exact = FALSE)
compensation_matrix <- read.csv(comp_dir, check.names = FALSE, row.names = 1)



reference_file <- read.FCS(paste(dir_raw, files[1], sep=""))
good_marker_names <- markernames(reference_file)


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
for (file in files[-6]){ 
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE) 
  markernames(ff) <- good_marker_names
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
                       channel_y = "FITC-A"))
  
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




# compensation ####
# ss_file_paths <- list.files("~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/comp_controls", pattern="*.fcs", full.names = TRUE)
# ss_files <- list.files("~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/comp_controls/", pattern="*.fcs", full.names = FALSE)
# 
# 
# frames <- lapply(ss_file_paths, function(x)read.FCS(x))
# # extract stain names from file names
# names(frames) <- substr(ss_files, regexpr("Controls_", ss_files)+9, regexpr("* Control.fcs", ss_files)-9)
# names(frames) <- paste(names(frames), "-A", sep="")
# names(frames) <- gsub("U-A", "UNSTAINED", names(frames))
# channel_names <-names(frames)
# 
# frames <- as(frames, "flowSet")
# frames <- fsApply(frames, function(x)x[,c("FSC-H", "SSC-H", channel_names[-7])])
# 
# comp <- spillover(frames,
#                   unstained="UNSTAINED",
#                   #patt="fluo_",
#                   fsc = "FSC-H",
#                   ssc = "SSC-H",
#                   stain_match = "regex")
# 
# write.csv(comp, "~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/comp_matrix.csv", quote = FALSE, row.names = TRUE)
# 
# # make diagnostic plot 
# 
# premessa_table <- premessa::read_parameters(list.files("~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/", pattern="*.fcs", full.names = TRUE, recursive = FALSE)[1])
# panel <- data.frame("fcs_colname"=rownames(premessa_table[1]),
#                     "marker_name"=premessa_table[[1]])
# 
# bcell_markers <- panel$marker_name[c(7:10, 12:14)]
# 
# 
# 
# 
# fcs_file_paths <- list.files("~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/", pattern="*.fcs", full.names = TRUE)
# b_cells <- flowCore::read.flowSet(as.character(fcs_file_paths), truncate_max_range = FALSE)
# b_cells_comp <- compensate(b_cells_comp, comp)
# 
# 
# for(i in bcell_markers[-3]){
#   
#   marker_biaxial <- ggcyto(b_cells_comp, aes(x = `CD19`, y = {{i}}))+
#     scale_x_logicle()+
#     scale_y_logicle()+
#     # ggcyto_par_set(limits = list(x = c(1,10e5), y= c(1,10e5)))+
#     geom_hex(bins=128)+
#     theme_minimal()
#   
#   ggsave(paste("~/postdoc/stanford/clinical_data/BC1/b_cell_data/raw_batch1/figures/comp/compensated_", i, "_cd19_biaxial.png", sep = ""), marker_biaxial , height=20, width=40, bg="white")
# }

