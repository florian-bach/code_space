library(flowStats)
library(flowCore)

ss_file_paths <- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/comp_controls/", pattern="*.fcs", full.names = TRUE)
ss_files <- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/comp_controls/", pattern="*.fcs", full.names = FALSE)
frames <- lapply(ss_file_paths, function(x)read.FCS(x))

# extract stain names from file names
names(frames) <- substr(files, regexpr("Controls_", files)+9, regexpr("* Control.fcs", files)-9)
names(frames) <- paste(names(frames), "-A", sep="")
names(frames) <- gsub("U-A", "UNSTAINED", names(frames))
frames <- as(frames, "flowSet")
frames <- fsApply(frames, function(x)x[,c(2,5,7:14)])


comp <- spillover(frames,
                  unstained="UNSTAINED",
                  #patt="fluo_",
                  fsc = "FSC-H",
                  ssc = "SSC-H",
                  stain_match = "regex")

write.csv(comp, "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/comp_matrix.csv", quote = FALSE, row.names = TRUE)

# make diagnostic plot 

panel <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/tfh_panel", header = T, stringsAsFactors = F)


tfh_markers <- panel$antigen[1:8]




fcs_file_paths <- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/", pattern="*.fcs", full.names = TRUE)
tfhs <- flowCore::read.flowSet(as.character(fcs_file_paths), truncate_max_range = FALSE)
tfhs_comp <- compensate(tfhs, comp)


for(i in tfh_markers[-1]){
  
  marker_biaxial <- ggcyto(tfhs_comp, aes(x = `CD3`, y = {{i}}))+
    scale_x_flowjo_biexp(pos=4.5, maxValue = 262144)+
    scale_y_flowjo_biexp(pos=4.5, maxValue = 262144)+
    ggcyto_par_set(limits = list(x = c(1,10e5), y= c(1,10e5)))+
    geom_hex(bins=512)+
    theme_minimal()
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/comp/compensated_", i, "_cd3_biaxial.png", sep = ""), marker_biaxial , height=20, width=40, bg="white")
}

