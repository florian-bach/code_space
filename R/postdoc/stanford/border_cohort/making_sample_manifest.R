files <- read.delim("/Volumes/fbach/all_fastq_files.txt", header = FALSE)
files <- gsub("./", "/labs/prasj/BIG_Flo/usftp21.novogene.com/01.RawData/", files$V1, fixed = TRUE)

df <- data.frame("read1-file-name"=files[grep("*1.fq.gz", files)],
                 "read2-file-name"=files[grep("*2.fq.gz", files)])

df$group <- substr(df$read1.file.name, 53, 61)
df$group <- gsub("/", "", df$group)

colnames(df)[3] <- "read.group.line"

colnames(df) <- gsub(".", "-", colnames(df), fixed = TRUE )

write.table(df, file = "/Volumes/fbach/sample_manifest.tsv", sep="\t", quote = FALSE, col.names = FALSE)
