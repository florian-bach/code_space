library(dplyr)

samer <- read.table("~/Downloads/test1.csv")

non_empty <- samer %>%
  select(Gene) %>%
  filter(Gene!="") 

grep("#VALUE!", non_empty, fixed=TRUE, value = TRUE)

write.table(non_empty, "samer_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
