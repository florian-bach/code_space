library(tidyr)
library(dplyr)

`%notin%` <- Negate(`%in%`)

# data generation ####
## reading in pilot data ####
musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/MASTER_METADATA.csv")
random_codes <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/id_date_code.csv")
random_codes$plasma.barcode <- gsub("D1PN8A", "D1_PN8A", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1JLGS", "DIJLGS", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1KWT2", "D1_KWT2", random_codes$plasma.barcode)
random_codes$plasma.barcode <- gsub("D1EF4A", "DEF4A", random_codes$plasma.barcode)
random_codes$date <- as.Date(random_codes$date)

slim_musical_metadata <- musical_metadata %>%
  mutate(day_annotation=if_else(day_annotation==84, -1, day_annotation))%>%
  select(combined_id, combined_date, enrolltype, day_annotation, gender_categorical, ageyrs, qpcr, TEMP)%>%
  mutate(id=combined_id, date=combined_date, class=enrolltype, timepoint=paste("t", day_annotation, sep=""), temperature=TEMP)%>%
  select(-combined_id, -combined_date, -enrolltype, -day_annotation)%>%
  mutate("study"="MUSICAL", date=as.Date(date), plate_number="pilot")


combo_frame2 <- merge(slim_musical_metadata, random_codes, by=c("id", "date"))

nulisa <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/pilot/nulisa_data.csv")

wide_nulisa <- nulisa %>%
  pivot_longer(cols = colnames(nulisa)[2:ncol(nulisa)], names_to = "plasma.barcode", values_to = "concentration")

wide_nulisa <- inner_join(wide_nulisa, combo_frame2, by="plasma.barcode")

pilot_nulisa <- wide_nulisa %>%
  mutate("time_class"=paste(class, timepoint, sep='_'),
         "age_class"=if_else(.$id %in% c(268, 324, 137, 176, 353, 161, 363, 571, 10766, 10794, 10842), "child", "adult"),
         id=factor(id))%>%
  filter(age_class=="child")%>%
  mutate(timepoint=case_when(timepoint=="t-1"~"baseline",
                             timepoint=="t0"~"day0",
                             timepoint=="t14"~"day14",
                             timepoint=="t7"~"day7"),
         timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")),
         barcode=plasma.barcode,
         sample_id=paste("X", id, time_class, barcode))%>%
  group_by(targetName) %>%
  mutate("z_conc"=scale(concentration, center = TRUE, scale = TRUE))%>%
  ungroup()%>%
  group_by(sample_id)%>%
  mutate("mean_z_conc"=mean(z_conc), "median_z_conc"=median(z_conc), study="pilot")%>%
  select(targetName, sample_id, concentration, barcode, id, gender_categorical, class, timepoint, z_conc, mean_z_conc, ageyrs, study, qpcr, temperature, plate_number)


## reading in big batch data ####
metadata <- readxl::read_excel("~/postdoc/stanford/plasma_analytes/MUSICAL/big_data/Immunology List _MUSICAL.xlsx")

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
                                 select(targetName, sample_id, concentration, plate_number))

big_df <- do.call(rbind, list_of_long_batches)

big_df1 <- big_df %>%
  mutate(barcode= substr(sample_id, nchar(sample_id)-5, nchar(sample_id)),
         id=factor(metadata$id[match(barcode, metadata$plasma1_barcode)]),
         gender_categorical=factor(metadata$gender_categorical[match(barcode, metadata$plasma1_barcode)]),
         class=metadata$infectiontype[match(barcode, metadata$plasma1_barcode)],
         ageyrs=metadata$ageyrs[match(barcode, metadata$plasma1_barcode)],
         timepoint_imm=metadata$timepoint_imm[match(barcode, metadata$plasma1_barcode)],
         qpcr=as.numeric(metadata$qpcr[match(barcode, metadata$plasma1_barcode)]),
         temperature=as.numeric(metadata$temperature[match(barcode, metadata$plasma1_barcode)])
  )%>%
  filter(barcode != "Repeat")%>%
  mutate("timepoint"=case_when(timepoint_imm==-2~"bad_baseline",
                               timepoint_imm==-1~"baseline",
                               timepoint_imm==0~"day0",
                               timepoint_imm==7~"day7",
                               timepoint_imm==14~"day14",
                               timepoint_imm==28~"day28")
         
  )%>%
  mutate(timepoint=factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  ungroup()%>%
  group_by(targetName) %>%
  mutate(z_conc=scale(concentration))%>%
  ungroup()%>%
  group_by(sample_id)%>%
  mutate(mean_z_conc=mean(z_conc), median_z_conc=median(z_conc), study="big_batch")%>%
  ungroup()%>%
  select(-timepoint_imm)



## putting in parasitemia data####

combo_data <- bind_rows(big_df1, pilot_nulisa)%>%
  mutate(id=paste(study, id, sep="_"),
         # timepoint=factor(timepoint, levels=c("baseline", "day0", "day7", "day14")), 
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


qpcr_cat_at_day0 <- combo_data %>%
  filter(class %in% c("A", "S"), !is.na(qpcr_cat))%>%
  group_by(class, id)%>%
  reframe("day0_qpcr_cat"=qpcr_cat[timepoint=="day0"],
          "day0_qpcr"=log10(qpcr+0.1)[timepoint=="day0"])%>%
  distinct()

unclean_data <- combo_data%>%
  left_join(., qpcr_cat_at_day0, by=c("class", "id"))

## putting it all together ####

write.csv(unclean_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/unclean_musical_combo_with_metadata.csv", row.names = FALSE)

clean_data <- unclean_data %>%
  filter(mean_z_conc > (-1.3))%>%
  # filter(sample_id %notin% c("X384_S.1_D1V93U",
  #                            "X744_A.1_D1KT5Z",
  #                            "X323_A14_D1DZ2D",
  #                            "X496_NM7_D1CAYS",
  #                            "X316_S.1_D12FNR",
  #                            "X667_NM7_D1GBYA",
  #                            "X 176 S_t14 D1FXRJ",
  #                            "X164_NM0_D1WA6Q",
  #                            "X 176 S_t7 D1VFF7",
  #                            "X219_S.1_D1C4ZM",
  #                            "X132_S0_D1AUZD"))%>%
  filter(barcode %notin% c("D19E2G", "D1FK67", "D1SSPJ"))%>%
  filter(sample_id %notin% c("X317_S.1_D1RTFU", "X316_S.1_D12FNR"))


write.csv(clean_data, "~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv", row.names = FALSE)

