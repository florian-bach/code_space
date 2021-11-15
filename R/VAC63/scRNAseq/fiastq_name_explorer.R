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
