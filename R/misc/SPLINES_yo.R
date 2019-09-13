if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("fission")


library(DESeq2)
library("fission")
data("fission")

se <- fission

ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("minute","strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))

yeast1 <- data.frame(spline(x=filter(fiss$strain=="wt")$minute, y=fiss$count, n=100))
yeast2 <- data.frame(spline(x=fiss$minute, y=fiss$count, n=100))

yeasts <- list(yeast1, yeast2)




yeet <- split(fiss, fiss$strain)
applied_yeet <- lapply(yeet, FUN = function(x){spline(x$minute, x$count, n=600)})  

ggplot()+
  geom_point(data=fiss, aes(x = minute, y = count, color = strain, group = strain))+
  #scale_y_log10()+
  #geom_line(data=hello, aes(x=x, y=y))
  #geom_line(data=dplyr::bind_rows(yeasts, .id="df"), aes(x=x, y=y, color=df))
  geom_line(data=dplyr::bind_rows(applied_yeet, .id="df"), aes(x=x, y=y, color=df))
  



hello1 <- data.frame(spline(x=fiss$minute, y=fiss$count))
