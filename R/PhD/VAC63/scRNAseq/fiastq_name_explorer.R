# slicing fastq file names 
files <- read.delim("~/postdoc/all_file_names.txt", sep="\t")

files <- subset(unlist(files, use.names = FALSE), grepl("fastq.gz", unlist(files, use.names = FALSE)))

name_df <- data.frame("lane"=substr(files, nchar(files)-29, nchar(files)-26), 
                      "s"=substr(files, nchar(files)-33, nchar(files)-31),
                      "AK"=substr(files, nchar(files)-40, nchar(files)-34),
                      "lib"=substr(files, 1, 4)
                      )

name_df$s <- gsub("_", "", name_df$s)
name_df$AK <- gsub("_", "", name_df$AK)
name_df$AK <- gsub("-", "", name_df$AK)
name_df$lib <- gsub("-", "", name_df$lib)

# making sample names for cellranger multi
txt_files <- list.files(path = "~/postdoc/scRNAseq/metadata/", pattern = ".txt", full.names = TRUE)

for(text_file in txt_files){

sample_names <- readLines(text_file, warn = FALSE)
sample_names <- scan(text=sample_names, what='', sep=' ')
sample_names <- gsub("[", "", sample_names, fixed = TRUE)
sample_names <- gsub("]", "", sample_names, fixed = TRUE)
sample_names <- gsub(",", "", sample_names)

sample_names <- sample_names[!duplicated(sample_names)]

new_name <- gsub(".txt", ".csv", text_file)


write.csv(sample_names, new_name, row.names = FALSE, col.names = NULL)}


# mkdir for each sample

L003 <- read.csv("/home/flobuntu/postdoc/scRNAseq/metadata/L003_gex_sample_names.csv")
print(L003$x, sep = " ", quote = FALSE)


# account for hashing in cellranger multi config csv

library(dplyr)

config <- read.csv("~/postdoc/scRNAseq/all_sample_names.csv")

config$sample <- ifelse(substr(config$fastq_id, 1,4) %in% c("GEX1", "GEX2", "VDJ1", "VDJ2"), "primary", "tertiary")

primary_config <- config %>%
  filter(sample=="primary") %>%
  select(-sample)

tertiary_config <- config %>%
  filter(sample=="tertiary") %>%
  select(-sample)


write.csv(primary_config, "~/postdoc/scRNAseq/primary_config_scratch.csv", row.names = FALSE)
write.csv(tertiary_config, "~/postdoc/scRNAseq/tertiary_config_scratch.csv", row.names=FALSE)


# cell surface names

# word_soup <- c("CS1-AK509_S43_L003_R1_001.fastq.gz,CS1-AK7617_S45_L003_R1_001.fastq.gz,CS1-X001_S47_L003_R1_001.fastq.gz,CS1-X042_S48_L003_R1_001.fastq.gz,CS2-AK14577_S25_L003_R1_001.fastq.gz,CS2-AK14578_S26_L003_R1_001.fastq.gz,CS2-AK14579_S27_L003_R1_001.fastq.gz,CS2-AK5796_S5_L003_R1_001.fastq.gz,CS3-AK13187_S2_L003_R1_001.fastq.gz,CS3-AK13188_S21_L003_R1_001.fastq.gz,CS3-AK7465_S44_L003_R1_001.fastq.gz,CS3-AK8767_S46_L003_R1_001.fastq.gz,CS4-AK12528_S6_L003_R1_001.fastq.gz,CS4-AK12529_S7_L003_R1_001.fastq.gz,CS4-AK12530_S8_L003_R1_001.fastq.gz,CS4-AK12531_S9_L003_R1_001.fastq.gz -R2 CS1-AK509_S43_L003_R2_001.fastq.gz,CS1-AK7617_S45_L003_R2_001.fastq.gz,CS1-X001_S47_L003_R2_001.fastq.gz,CS1-X042_S48_L003_R2_001.fastq.gz,CS2-AK14577_S25_L003_R2_001.fastq.gz,CS2-AK14578_S26_L003_R2_001.fastq.gz,CS2-AK14579_S27_L003_R2_001.fastq.gz,CS2-AK5796_S5_L003_R2_001.fastq.gz,CS3-AK13187_S2_L003_R2_001.fastq.gz,CS3-AK13188_S21_L003_R2_001.fastq.gz,CS3-AK7465_S44_L003_R2_001.fastq.gz,CS3-AK8767_S46_L003_R2_001.fastq.gz,CS4-AK12528_S6_L003_R2_001.fastq.gz,CS4-AK12529_S7_L003_R2_001.fastq.gz,CS4-AK12530_S8_L003_R2_001.fastq.gz,CS4-AK12531_S9_L003_R2_001.fastq.gz,CS1-AK509_S43_L004_R1_001.fastq.gz,CS1-AK7617_S45_L004_R1_001.fastq.gz,CS1-X001_S47_L004_R1_001.fastq.gz,CS1-X042_S48_L004_R1_001.fastq.gz,CS2-AK14577_S25_L004_R1_001.fastq.gz,CS2-AK14578_S26_L004_R1_001.fastq.gz,CS2-AK14579_S27_L004_R1_001.fastq.gz,CS2-AK5796_S5_L004_R1_001.fastq.gz,CS3-AK13187_S2_L004_R1_001.fastq.gz,CS3-AK13188_S21_L004_R1_001.fastq.gz,CS3-AK7465_S44_L004_R1_001.fastq.gz,CS3-AK8767_S46_L004_R1_001.fastq.gz,CS4-AK12528_S6_L004_R1_001.fastq.gz,CS4-AK12529_S7_L004_R1_001.fastq.gz,CS4-AK12530_S8_L004_R1_001.fastq.gz,CS4-AK12531_S9_L004_R1_001.fastq.gz -R2 CS1-AK509_S43_L004_R2_001.fastq.gz,CS1-AK7617_S45_L004_R2_001.fastq.gz,CS1-X001_S47_L004_R2_001.fastq.gz,CS1-X042_S48_L004_R2_001.fastq.gz,CS2-AK14577_S25_L004_R2_001.fastq.gz,CS2-AK14578_S26_L004_R2_001.fastq.gz,CS2-AK14579_S27_L004_R2_001.fastq.gz,CS2-AK5796_S5_L004_R2_001.fastq.gz,CS3-AK13187_S2_L004_R2_001.fastq.gz,CS3-AK13188_S21_L004_R2_001.fastq.gz,CS3-AK7465_S44_L004_R2_001.fastq.gz,CS3-AK8767_S46_L004_R2_001.fastq.gz,CS4-AK12528_S6_L004_R2_001.fastq.gz,CS4-AK12529_S7_L004_R2_001.fastq.gz,CS4-AK12530_S8_L004_R2_001.fastq.gz,CS4-AK12531_S9_L004_R2_001.fastq.gz")
# word_soup2 <- gsub(",", "", "", word_soup)

oligo_names <- t(read.csv("~/postdoc/scRNAseq/metadata/cs_names.csv", stringsAsFactors = FALSE, header = FALSE))
oligo_names <- data.frame(oligo_names)

oligo_names$sample <- ifelse(substr(oligo_names$oligo_names, 1,3) %in% c("CS1", "CS2"), "primary", "tertiary")
oligo_names$read <- ifelse(grepl("R1", oligo_names$oligo_names), "R1", "R2")


primary_oligo_R1 <- oligo_names %>%
  filter(sample=="primary", read=="R1") %>%
  select(oligo_names)

primary_oligo_R2 <- oligo_names %>%
  filter(sample=="primary", read=="R2") %>%
  select(oligo_names)

tertiary_oligo_R1 <- oligo_names %>%
  filter(sample=="tertiary", read=="R1") %>%
  select(oligo_names)

tertiary_oligo_R2 <- oligo_names %>%
  filter(sample=="tertiary", read=="R2") %>%
  select(oligo_names)


write.table(primary_oligo_R1, "~/postdoc/scRNAseq/primary_cs_R1_names.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(tertiary_oligo_R1, "~/postdoc/scRNAseq/tertiary_cs_R1_names.txt", row.names = FALSE, col.names = FALSE,  sep = "\t", quote = FALSE)

write.table(primary_oligo_R2, "~/postdoc/scRNAseq/primary_cs_R2_names.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(tertiary_oligo_R2, "~/postdoc/scRNAseq/tertiary_cs_R2_names.txt", row.names = FALSE, col.names = FALSE,  sep = "\t", quote = FALSE)
