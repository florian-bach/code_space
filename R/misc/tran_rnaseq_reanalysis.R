library(DESeq2)
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE52166", "file=GSE52166_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
ref_tbl <- tbl 

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID


gsms <- paste0("XXXXXXXXXXXXXXXX00X2X2X0X0X2X2X01X2X1X0X2X2X1XXXX0",
               "XXX0X1X0X0X0X2X0XXX2X0X2X1X1X2X2X1X1X2X0X1X2XXX12X",
               "1X1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# filter lowly expressed genes 
# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
gs <- factor(sml)
groups <- make.names(c("immune","early fever","delayed fever"))
levels(gs) <- groups

keep <- rowSums(tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

# didn't see any DEG when variance selecting across whole data... let's redo on subset
tbl_df <- data.frame(tbl) 
normed_tbl_df <- as.data.frame(edgeR::cpm(tbl_df))
normed_tbl_df$variance <- apply(normed_tbl_df[,1:ncol(tbl)], 1, function(x)var(x))
normed_tbl_df$mean <- apply(normed_tbl_df[,1:ncol(tbl)-1], 1, function(x)mean(x))

normed_tbl_df$normed_variance <- normed_tbl_df$variance/normed_tbl_df$mean

# tbl_df$sd <- apply(tbl_df[,1:ncol(tbl)], 1, function(x)sd(x))
# tbl_df$min <- apply(tbl_df[,1:ncol(tbl)], 1, function(x)min(x))
# tbl_df$max <- apply(tbl_df[,1:ncol(tbl)], 1, function(x)max(x))
# make sure the root of the variance is the sd
# tbl_df$test <- sqrt(tbl_df$variance)

ordered_tbl_df <- normed_tbl_df[order(normed_tbl_df$variance, decreasing = TRUE),]

most_variable <- ordered_tbl_df[1:6834,]

subset_tbl <- subset(tbl, rownames(tbl)%in%rownames(most_variable))




# group membership for samples
gs <- factor(sml)
groups <- make.names(c("immune","early fever","delayed fever"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(subset_tbl))

ds <- DESeqDataSetFromMatrix(countData=subset_tbl, colData=sample_info, design= ~Group)

cts <- list(c("Group",groups[1],groups[2]),
            c("Group",groups[1],groups[3]),
            c("Group",groups[2],groups[3]))  # contrasts of interest

ds <- DESeq(ds, test="LRT", reduced = ~ 1)  # Use LRT for all-around gene ranking

# extract results for top genes table
r <- results (ds, alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:250],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","stat","baseMean","Symbol","Description"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

plotDispEsts(ds, main="GSE52166 Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "GSE52166 Frequencies of padj-values")

# Wald test to obtain contrast-specific results
ds <- DESeq(ds, test="Wald", sfType="poscount")
r <- results(ds, contrast=cts[[1]], alpha=0.05, pAdjustMethod = "fdr")

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(groups[1], "vs", groups[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
# with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0.5849625),
     with(subset(r, padj<0.05 ),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)


all_res <- list()

for (ct in cts) {
  i <- length(all_res)
  r <- results(ds, contrast=ct, alpha=0.05, pAdjustMethod = "fdr")
  # all_res[[i + 1]] <- rownames(r)[!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) >= 0.5849625] 
  all_res[[i + 1]] <- rownames(r)[!is.na(r$pvalue) & r$pvalue < 0.05 ] 
  names(all_res)[i + 1] <- paste(ct[2:3], collapse="_")
}

gplots::venn(all_res)
