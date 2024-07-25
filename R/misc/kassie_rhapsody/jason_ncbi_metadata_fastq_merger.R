old <- read.csv("/Volumes/fbach/jason_ncbi/old_sra_metadata.csv")
old$library_ID <- gsub("_S", "_T", old$library_ID)

old$sample_name <- gsub("3528_T2", "3528_T1", old$sample_name)
old$sample_name <- gsub("3528_T3", "3528_T2", old$sample_name)

old$library_ID <- gsub("3528_T2", "3528_T1", old$library_ID)
old$library_ID <- gsub("3528_T3", "3528_T2", old$library_ID)


fastq_locations <- read.csv("/Volumes/fbach/jason_ncbi/fixed_fastq_locations.csv")
fastq_locations$library_ID <- paste(fastq_locations$subject_id, fastq_locations$timepoint, fastq_locations$library_type, sep="_")

#3528_T1 to be removed
#3481_T4_TCR to be removed
#3394 BAM only
#only include 3481_T4_TCR, remove 3481_T4_GEX, 3481_T4_HASH

# add these to metadata and we're gucci
# /labs/prasj/Jason/CD4_T_Cell_Sequencing/Raw_Sequencing/3410_T1-T4_RNA-TCR/GEX_HASH_Raw
# /labs/prasj/Jason/CD4_T_Cell_Sequencing/Raw_Sequencing/3410_T1-T4_RNA-TCR/TCR_Raw

#change metadata for 3528 so that T3 becomes T2; T2 is the real T1; original T1 was crap
new <- fastq_locations %>%
mutate("sample_name" = paste(subject_id, timepoint, sep="_"),
       "library_ID" = ifelse(timepoint!="", paste(subject_id, timepoint, library_type, sep="_"), paste(subject_id, library_type, sep="_")),
       "title"="RNA-TCR-seq of homo sapien: CD4 T cells",
       "library_strategy"="RNA-Seq",
       "library_source"="TRANSCRIPTOMIC SINGLE CELL",
       "library_selection"="polyA",
       "library_layout"="paired",
       "platform"="Illumina",
       "instrument_model"=old$instrument_model[match(fastq_locations$library_ID, old$library_ID)],
       "design_description"="Library preparation performed using reagents from 10X genomics according to their published protocol (CG000330 Rev D)",
       "filetype"="fastq.gz",
       "assembly"="")

#which libraries aren't in the old metadata
unique(new$library_ID)[unique(new$library_ID) %notin% unique(old$library_ID)]
# [1] "3354_S_GEX"  "3354_K_HASH" "3354_S_HASH" "3354_S_TCR"  "3481_T4_TCR" "_"

# remove
new_small <- subset(new, library_ID %notin% unique(new$library_ID)[unique(new$library_ID) %notin% unique(old$library_ID)])



#which old metadata files aren't in the new metadata
unique(old$library_ID)[unique(old$library_ID) %notin% unique(new_small$library_ID)]
# "3394_T1_GEX"  "3394_T2_GEX"  "3394_T3_GEX" <- BAM
# "3394_T1_TCR" # "3394_T2_TCR"  "3394_T3_TCR"  "3394_T1_HASH" "3394_T2_HASH" "3394_T3_HASH" <- BAM
# "3410_T1_GEX"  "3410_T2_GEX"  "3410_T3_GEX"  "3410_T4_GEX" <- add new folder
# "3410_T1_HASH" "3410_T2_HASH" "3410_T3_HASH" "3410_T4_HASH" <- add new folder


#### metadata file should have one sample per row, and then forward & reverse reads as additional columns
write.csv(new_small, "/Volumes/fbach/jason_ncbi/new_sra_metadata.csv", row.names = FALSE)


# read completed metadata

new_complete <- read.csv("/Volumes/fbach/jason_ncbi/complete_new_sra_metadata2.csv")

# new_complete_wide <- new_complete %>%
#   group_by(library_ID)%>%
#   mutate(ffile=paste("file_name", seq(1, n()), sep=""))%>%
#   pivot_wider(names_from = ffile, values_from = file_path)%>%
#   mutate("organism"="Homo sapiens")
# 
new_complete_paths_only <- new_complete$file_path
write.csv(new_complete_paths_only, "/Volumes/fbach/jason_ncbi/new_complete_fastq_paths_only.csv", row.names = FALSE, col.names = NULL, quote = FALSE)


new_complete_wide <- new_complete %>%
  # mutate(sample_name=library_ID)%>%
  # dplyr::select(-library_ID)%>%
  mutate("filename"=file_path)%>%
  dplyr::select(-file_path)%>%
  group_by(library_ID)%>%
  mutate(ffile=paste("filename", seq(1, n()), sep=""))%>%
  pivot_wider(names_from = ffile, values_from = filename)%>%
  mutate("filename"=filename1, .before = "filename2")%>%
  dplyr::select(-filename1)

#rename sorted shit
new_complete_wide$sample_name <- gsub("Sort_non_", "Sort_non", new_complete_wide$sample_name)
new_complete_wide$sample_name <-   gsub("Sort_Tr1_", "Sort_Tr1", new_complete_wide$sample_name)
  
  
  
write.table(data, "/Volumes/fbach/jason_ncbi/new_complete_wide_na_removed_fastq_instrument2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(new_complete_wide, "/Volumes/fbach/jason_ncbi/new_complete_wide.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



# 
# data <- read.table("/Volumes/fbach/jason_ncbi/new_complete_wide_na_removed_fastq_instrument2.tsv", header=TRUE, sep="\t")
# data$instrument_model <- gsub("GRCh38-v7.0.0", "Illumina NovaSeq 6000", data$instrument_model)
# 
# write.table(data, "/Volumes/fbach/jason_ncbi/new_complete_wide_na_removed_fastq_instrument2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

old_table <- read.table("/Volumes/fbach/jason_ncbi/new_complete_wide_na_removed_fastq_instrument.tsv", sep = "\t")


