library(flowCore)
library(ggcyto)

md <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/metadata.csv", header=T, stringsAsFactors = TRUE)

cleaned_tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
cleaned_tfh$idt <- paste("1", substr(cleaned_tfh$idt, 15, 18), sep='')
drop_removed_tfh <- subset(cleaned_tfh, tfh_drop=="No")

md <- subset(md, id %in% drop_removed_tfh$idt)

panel <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/tfh_panel", header = T, stringsAsFactors = F)


tfh_markers <- panel$antigen[1:8]

tfhs <- flowCore::read.flowSet(as.character(md$file_paths), truncate_max_range = FALSE)

channel_names <- markernames(tfhs[[1]])

markernames(tfhs) <- channel_names

for(i in tfh_markers){
  
  marker_histo <- ggcyto(tfhs, aes(x = {{i}}))+scale_x_flowCore_fasinh()+geom_histogram(bins=444)

  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", i, "_histo.png", sep = ""), marker_histo, height=20, width=40)
}

for(i in tfh_markers[-1]){
  
  marker_biaxial <- ggcyto(tfhs, aes(x = `CD3`, y = {{i}}))+
    scale_x_logicle()+
    scale_y_logicle()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", i, "_cd3_biaxial.png", sep = ""), marker_biaxial , height=20, width=40)
}

for(i in tfh_markers){
  
  marker_biaxial <- ggcyto(tfhs, aes(x = `Time`, y = {{i}}))+
    scale_x_logicle()+
    scale_y_flowCore_fasinh()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", i, "_time_biaxial.png", sep = ""), marker_biaxial , height=20, width=40)
}


  biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CD4"))+
    scale_x_flowCore_fasinh()+
    scale_y_flowCore_fasinh()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CD4", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
  

  
  biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CCR6"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CCR6", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
  
  biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CXCR5"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CXCR5", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
  
  biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "PD1"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "PD1", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
  
  biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CXCR3"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CXCR3", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
  
   biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CCR4"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CCR4", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
 
   biaxial <- ggcyto(tfhs, aes(x = "CD3", y = "CD45RA"))+
    scale_x_log10()+
    scale_y_log10()+
    geom_hex(bins=128)
  
  ggsave(paste("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/", "CD45RA", "_cd3_biaxial.png", sep = ""), biaxial, height=20, width=40)
 
# 
# channel_names <- fsApply(tfhs, function(x) getChannelMarker(frm = x, name = "FJComp-AmCyan-A"), simplify=TRUE)
# channel_df <- as.data.frame(unlist(channel_names))
# colnames(channel_df) <- "value"
# channel_df$class <- substr(rownames(channel_df), 44, 48)
# channel_df$file <- substr(rownames(channel_df), 1,42)
# channel_df <- tidyr::pivot_wider(channel_df, names_from = class, values_from = value, id_cols = file)
# 
# nas <- unlist(subset(channel_df, is.na(desc), select = file))
# 
# test <- tfhs[1]
# markernames(test)[1] <- "yo_mama"
# 
# thfs2 <- fsApply(tfhs, function(x) assign(markernames(x)["FJComp-AmCyan-A"], "CD3"))

channel_names <- markernames(tfhs[[1]])

markernames(tfhs) <- channel_names






markernames(tfhs[nas])

cxcr5_histo <- ggcyto(tfhs, aes(x = `CXCR5`))+scale_x_flowjo_biexp()+geom_histogram()
ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/cxcr5_histo.png", cxcr5_histo, height=20, width=40)



grep("export_T_FH_PANEL_B1_1165_Single_Cells.fcs", md$file_names)


# ??????
# for (file in md$file_paths){
#   ff <- read.FCS(file, truncate_max_range = FALSE)
#   translist <- estimateLogicle(ff, colnames(ff)[7:14])
#   ff_t <- flowCore::transform(ff, translist)
#   write.FCS(ff, filename = paste0(substr(file, 95, 98), "_transformed.fcs", sep=""))
# }


