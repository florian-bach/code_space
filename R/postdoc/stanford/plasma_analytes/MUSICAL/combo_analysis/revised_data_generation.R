library(tidyr)
library(dplyr)

`%notin%` <- Negate(`%in%`)

# read in pilot data ####

pilot_nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/nulisa_data.csv")

long_pilot_nulisa <- pilot_nulisa %>%
  pivot_longer(cols = colnames(nulisa)[2:ncol(nulisa)], names_to = "plasma1_barcode", values_to = "concentration")%>%
  mutate(plasma1_barcode = case_when(plasma1_barcode=="D1_PN8A"~"D1PN8A",
                                    plasma1_barcode=="DIJLGS"~"D1JLGS",
                                    plasma1_barcode=="D1_KWT2"~"D1KWT2",
                                    plasma1_barcode=="DEF4A"~"D1EF4A",
                                    .default = plasma1_barcode))%>%
  filter(grepl("^D1*", plasma1_barcode))%>%
  mutate(plate_number="pilot")

# read in big batch ####
list_of_files <- list.files(path = "~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/", pattern = "Report.csv", full.names = TRUE)

list_of_batches <- list(vector(length = length(list_of_files)))
# list_of_batches <- list()

system.time(for(i in 1:length(list_of_files)){
  
  tmp <- read.csv(list_of_files[i], header = TRUE)
  tmp$plate_number <- paste("plate", i)
  list_of_batches[[i]] <- tmp
  
})

list_of_long_batches <- lapply(list_of_batches, function(x) x %>%
                                 pivot_longer(cols= starts_with("X"),
                                              names_to = "sample_id",
                                              values_to = "concentration")%>%
                                 mutate(plasma1_barcode = substr(sample_id, nchar(sample_id)-5, nchar(sample_id)))%>%
                                 select(targetName, plasma1_barcode, concentration, plate_number))


# combine all batches and pilot ####
big_batch <- do.call(rbind, list_of_long_batches)

combo_batch <- rbind(big_batch, long_pilot_nulisa)

# add metadata ####
# one sample that we ran during NULISA pilot is missing from the big metadata file: 363 S day14 D1XM73; need to fix
metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

combo_with_meta <- inner_join(combo_batch, metadata, by = "plasma1_barcode")

# add z scores, tidy up
unclean_data <- combo_with_meta %>%
  mutate("timepoint"= case_when(timepoint_imm==-2~"bad_baseline",
                             timepoint_imm==-1~"baseline",
                             timepoint_imm==0~"day0",
                             timepoint_imm==7~"day7",
                             timepoint_imm==14~"day14",
                             timepoint_imm==28~"day28"))%>%
  mutate(sample_id=paste(id, timepoint, infectiontype))%>%
  group_by(targetName) %>%
  mutate(z_conc=scale(concentration))%>%
  ungroup()%>%
  group_by(plasma1_barcode)%>%
  mutate(mean_z_conc=mean(z_conc), median_z_conc=median(z_conc))%>%
  mutate(id=paste(plate_number, id, sep="_"),
         timepoint=factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")), 
         qpcr=as.numeric(qpcr),
         log_qpcr = log10(qpcr+0.01),
         qpcr_cat=case_when(qpcr==0~"0", 
                            qpcr>0 & qpcr<10~ ">1",
                            qpcr>10 & qpcr<100~ ">10",
                            qpcr>100 & qpcr<1000~ ">10e2",
                            qpcr>1000 & qpcr<10000~ ">10e3",
                            qpcr>10000 & qpcr<100000~ ">10e4",
                            qpcr>100000 & qpcr<1000000~ ">10e5",
                            qpcr>1000000 & qpcr<10000000~ ">10e6",
                            qpcr>10000000 ~ ">10e7"),
         qpcr_cat=factor(qpcr_cat, levels=c("0", ">1", ">10", ">10e2", ">10e3", ">10e4", ">10e5", ">10e6", ">10e7")),
         temperature = as.numeric(temperature),
         temperature_cat = case_when(temperature < 38 ~ "<38",
                                     temperature >=38 & temperature <39 ~ ">=38",
                                     temperature >39 & temperature <40 ~ ">39",
                                     temperature >40 & temperature <41 ~ ">40",
                                     temperature >41 & temperature <42 ~ ">41"))
           
clean_data <- unclean_data %>%
  # filter(mean_z_conc > (-1.3))%>%
  # filter(plasma1_barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  # filter(sample_id %notin% c("X317_S.1_D1RTFU", "X316_S.1_D12FNR"))%>%
  filter(sample_id %notin% c("164 day0 NM",
                             "316 baseline S",
                             "317 baseline S",
                             "410 basleine A",
                             "164 day7 NM",
                             "132 day0 S",
                             "495 day14 A",
                             "219 baseline S",
                             "176 day7 S",
                             "495 day0 S",
                             "176 day14 S",
                             "323 day14 A",
                             "324 day7 S",
                             "744 baseline A",
                             "384 baseline S")
                             )


# barcodes that are in the new combodata that werent in the previous: not sure why were missing
# unique(clean_combo_with_meta2$plasma1_barcode[clean_combo_with_meta2$plasma1_barcode %notin% clean_data$barcode])
# [1] "D1RTFU" "D12FNR" "D1KWT2" "D1PN8A" "D1JLGS"

# barcodes that were in clean data but aren't in new one: these make sense, need to fix D1XM73
# unique(clean_data$barcode[clean_data$barcode %notin% clean_combo_with_meta2$plasma1_barcode])
# [1] "D1_KWT2" "D1_PN8A" "D1XM73"  "DIJLGS" 


write.csv(clean_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv", row.names = FALSE)

write.csv(unclean_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/unclean_musical_combo_with_metadata.csv", row.names = FALSE)
