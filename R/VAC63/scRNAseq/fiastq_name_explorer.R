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
# 