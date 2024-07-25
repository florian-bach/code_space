missing_files <- read.table("~/Downloads/SUB14403133_missing_files_table.tsv", sep='\t')

simon_kattria <- grep("*Simon_Kattria*", missing_files$V1, value=TRUE)
t1_t3 <- grep("*3481_T1-T3_RNA-TCR*", missing_files$V1, value=TRUE)

fixed_simon_kattria <- gsub("T1", "T5", simon_kattria)
fixed_simon_kattria <- gsub("T2", "T6", fixed_simon_kattria)

fixed_t1_t3 <- gsub("T1-T3", "T1-T7", t1_t3)

missing_files$fixed <- c(fixed_simon_kattria, fixed_t1_t3)
colnames(missing_files)[1]="old"


big_table <- read.table("/Volumes/fbach/jason_ncbi/new_complete_wide_na_removed_fastq_instrument.tsv", sep = "\t", header=TRUE)

long_big_table <- big_table %>%
  pivot_longer(cols=starts_with("filename"), names_to="file_number", values_to="file_name")

long_big_table$file_name <- ifelse(long_big_table$file_name %in% missing_files$old, 
                                   missing_files$fixed[match(long_big_table$file_name, missing_files$old)],
                                   long_big_table$file_name)


old_old <- read.table("~/Downloads/new_complete_wide_na_removed_fastq_instrument2.tsv", sep='\t', header = TRUE)


wide_big_table <- long_big_table %>%
  pivot_wider(names_from = file_number, values_from = file_name)%>%
  mutate(sample_name=gsub(" ", "_", sample_name),
         sample_name=gsub("_$", "", sample_name),
         "instrument_model"=old_old$instrument_model[match(library_ID, old_old$library_ID)])





# write.csv(wide_big_table, "~/Downloads/wide_fixed_table.csv", quote = FALSE, row.names = FALSE)
# write.table(wide_big_table, "/Volumes/fbach/jason_ncbi/wide_sra_metadata_fixed.tsv", quote = FALSE, row.names = FALSE, sep="\t")

write.table(wide_big_table, "~/Downloads/wide_fixed_table.tsv", quote = FALSE, row.names = FALSE, sep="\t")






fixed_fastq_file_paths <- missing_files$fixed
write.table(fixed_fastq_file_paths, "/Volumes/fbach/jason_ncbi/fixed_missing_fastqs.tsv", quote = FALSE, row.names = FALSE, sep="\t")

fixed_fastq_file_paths <- read.table("/Volumes/fbach/jason_ncbi/fixed_missing_fastqs.tsv", sep="\t")
# grep("HT6KWDSX5", fixed_fastq_file_paths$V1)
more_fixed_fastq_file_paths <- gsub("HT6KWDSX5", "HT2KWDSX5", fixed_fastq_file_paths$V1)
write.table(more_fixed_fastq_file_paths, "/Volumes/fbach/jason_ncbi/more_fixed_missing_fastqs.tsv", quote = FALSE, row.names = FALSE, sep="\t")















big_table <- read.table("~/Downloads/wide_fixed_table.tsv", sep="\t", header = TRUE)

long_big_table <- big_table %>%
  pivot_longer(cols=starts_with("filename"), names_to="file_number", values_to="file_name")

long_big_table$file_name <- gsub("T6_HASH_CKDL230006713-1A_HT6KWDSX5", "T6_HASH_CKDL230006713-1A_HT2KWDSX5", long_big_table$file_name) 
long_big_table$file_name<- gsub("T5_HASH_CKDL230006712-1A_HT6KWDSX5", "T5_HASH_CKDL230006712-1A_HT2KWDSX5", long_big_table$file_name)


old_old <- read.table("~/Downloads/new_complete_wide_na_removed_fastq_instrument2.tsv", sep='\t', header = TRUE)


wide_big_table <- long_big_table %>%
  pivot_wider(names_from = file_number, values_from = file_name)%>%
  mutate(sample_name=gsub(" ", "_", sample_name),
         sample_name=gsub("_$", "", sample_name),
         "instrument_model"=old_old$instrument_model[match(library_ID, old_old$library_ID)])

write.table(wide_big_table, "/Volumes/fbach/jason_ncbi/wide_fixed_table2.tsv", quote = FALSE, row.names = FALSE, sep="\t")


new_fixed <- read.table("/Volumes/fbach/jason_ncbi/wide_fixed_table2.tsv", sep="\t", header = TRUE)

long_big_table2 <- new_fixed %>%
  pivot_longer(cols=starts_with("filename"), names_to="file_number", values_to="file_name")

final_fixed_paths <- long_big_table2$file_name[c(
  grep("T5_HASH_CKDL230006712-1A_HT6KWDSX5", long_big_table$file_name),
  grep("T6_HASH_CKDL230006713-1A_HT6KWDSX5", long_big_table$file_name))]

write.table(final_fixed_paths, "/Volumes/fbach/jason_ncbi/final_fixed_paths.tsv", quote = FALSE, row.names = FALSE, sep="\t")






# final fix?

all_files <- data.frame("full_paths"=unique(long_big_table2$file_name))
all_files$basenames <- noquote(basename(all_files$full_paths))
all_files$duplicated <- duplicated(all_files$basenames)

all_files_there <- as.data.frame(read.delim("/Volumes/fbach/jason_ncbi/all_files.txt", sep="\t", header = FALSE))
all_files_there$V1 <- gsub(" ", "", all_files_there)
all_files_there$V2 <- paste(all_files_there$V1)

missing_files <- all_basenames$basename[all_basenames$basename %notin% all_files_there$V2]

duplicates <- subset(all_files, all_files$duplicated)

all_files$basenames[match(duplicates$basenames, all_files$basenames)]

# these files are uploaded but aren't supposed to be 
# all_files_there$V2[!(all_files_there$V2 %in% all_basenames$basename)]
# [1] "H2_CKDL230013093-1A_H7377DSX7_S3_L003_I1_001.fastq.gz"
# [2] "H2_CKDL230013093-1A_H7377DSX7_S3_L003_I2_001.fastq.gz"
# [3] "H2_CKDL230013093-1A_H7377DSX7_S3_L003_R1_001.fastq.gz"
# [4] "H2_CKDL230013093-1A_H7377DSX7_S3_L003_R2_001.fastq.gz"

multiplication <- all_files[duplicated(all_files$basenames)|duplicated(all_files$basenames, fromLast = TRUE),]

multiplication$new_paths <- ifelse(grepl("Resequencing_Round1", multiplication$full_paths), paste(dirname(multiplication$full_paths), "/Resequence_1_", basename(multiplication$full_paths), sep=""),
                                   ifelse(grepl("Resequencing_Round2", multiplication$full_paths), paste(dirname(multiplication$full_paths), "/Resequence_2_", basename(multiplication$full_paths), sep=""),
                                          ifelse(grepl("gex_possorted_bam.bam", multiplication$full_paths), paste(dirname(multiplication$full_paths), "/A",1:3, "_", basename(multiplication$full_paths), sep=""),
                                                 ifelse(grepl("possorted_genome_bam", multiplication$full_paths), paste(dirname(multiplication$full_paths), "/T",1:3, "_", basename(multiplication$full_paths), sep=""), ":)")
                                          )))


write.csv(multiplication[, c(1,4)], "/Volumes/fbach/jason_ncbi/renamed_file_paths.csv", row.names = FALSE, quote = FALSE)


# rewriting SRA metadata to reflect the renamed files, & basenames ####

big_table <- read.table("/Volumes/fbach/jason_ncbi/wide_fixed_table2.tsv", header=TRUE, sep="\t")

new_big_table <- big_table %>%
  pivot_longer(cols=starts_with("filename"), names_to="file_number", values_to="file_name")%>%
  mutate(file_name = ifelse(file_name %in% multiplication$full_paths,
                          basename(multiplication$new_paths[match(file_name, multiplication$full_paths)]),
                           basename(file_name)))%>%
  pivot_wider(names_from = file_number, values_from = file_name)
  # 

write.table(new_big_table, "/Volumes/fbach/jason_ncbi/sra_metadata_basenames.tsv", quote = FALSE, row.names = FALSE, sep="\t")

write.table(whole_table$new_paths, "/Volumes/fbach/jason_ncbi/twenty_two_missing.tsv", quote = FALSE, row.names = FALSE, sep="\t")

# remove duplicate files

doubles <- read.csv("/Volumes/fbach/jason_ncbi/duplicate_submission_files.csv", header=TRUE)


long_doubles <- doubles %>%
  pivot_longer(cols=starts_with("Files"), names_to = "file_number", values_to = "file_names")

#dim(long_doubles) 32  4

shorter_doubles <- long_doubles[!grepl("Resequence_2", long_doubles$file_names),]
shorter_doubles <- shorter_doubles[!grepl("GEX_5_CKDL230013086-1A_H7377DSX7_S2", shorter_doubles$file_names),]
shortest_doubles <- shorter_doubles[!grepl("T2_CKDL230013088-1A_H7377DSX7_S2", shorter_doubles$file_names),]

#change sra metadata to remove duplicate files

big_table <- read.table("/Volumes/fbach/jason_ncbi/sra_metadata_basenames.tsv", header=TRUE, sep="\t")
long_big_table <- big_table %>%
  pivot_longer(cols=starts_with("filename"), names_to="file_number", values_to="file_name")
# > dim(long_big_table)
# [1] 1440   17
no_doubles_big_table <- long_big_table %>%
  filter(file_name %notin% shortest_doubles$file_names)%>%
  pivot_wider(names_from = file_number, values_from = file_name)
# 

write.table(no_doubles_big_table, "/Volumes/fbach/jason_ncbi/sra_metadata_no_doubles.tsv", quote = FALSE, row.names = FALSE, sep="\t")

write.csv(shortest_doubles$file_names, "/Volumes/fbach/jason_ncbi/duplicate_files_to_remove.csv")
