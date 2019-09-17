suppressPackageStartupMessages(library(ImpulseDE2))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(ggplot2))

library(gridExtra)

## simulate data
lsSim <- simulateDataSetImpulseDE2(
  vecTimePointsA=rep(1:6, 3),
  vecTimePointsB=NULL,
  vecBatchesA=rep("B_NULL", 6*3),
  vecBatchesB=NULL,
  scaNConst=100,
  scaNImp=50,
  scaNLin=50,
  scaNSig=50,
  scaNRand=0,
  scaMumax=1000,
  scaSeedInit = 1)
## run DESeq2
dfAnnotationDESeq2 <- lsSim$dfAnnotation

# create a natural cubic spline basis with 4 degrees of freedom
matTimeSplineBasis <- ns(
  dfAnnotationDESeq2$Time, df=4)

#  the "spline" is whatever linear combination of these four curves is
# determined by the regression fitting process.
# Another way to put it: the spline consists of a vector space of curves
# that can be created by taking linear combinations of these four curves.

colnames(matTimeSplineBasis) <- 
  paste0("spline", seq(1, dim(matTimeSplineBasis)[2]))


# add spline basis into annotation data frame
dfAnnotationDESeq2 <- cbind(dfAnnotationDESeq2, matTimeSplineBasis)
# create DESeq2 object with spline basis vectors as individual linear predictors

dds <- DESeqDataSetFromMatrix(
  countData = lsSim$matObservedCounts,
  colData = dfAnnotationDESeq2,
  design = ~spline1 + spline2 + spline3 + spline4)


# estimate dispersion with less stringent dispersion outlier handling
# the outlier calling is controlled via outlierSD of estimateDispersionsMAP()
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds)

## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.

dds <- estimateDispersionsMAP(dds, outlierSD = 10)
# perform a log-likelihood ratio test
dds <- nbinomLRT(dds, full = ~spline1 + spline2 + spline3 + spline4,
                 reduced = ~1)
# extract results
print(head( results(dds) ))


count_data <- data.frame(lsSim$matObservedCounts)
count_data <- dplyr::mutate(count_data, Gene=rownames(count_data))
count_data <- dplyr::select(count_data, Gene, colnames(count_data)[1:18])

long_data <- tidyr::gather(count_data, Sample, Count, colnames(count_data[2:ncol(count_data)])) 


long_data$Time <- dfAnnotationDESeq2[match(long_data$Sample, dfAnnotationDESeq2$Sample), 3]
long_data$ID <- substr(long_data$Sample, nchar(long_data$Sample)-3, nchar(long_data$Sample))  


long_data$p_value <- dfAnnotationDESeq2[match(long_data$Sample, dfAnnotationDESeq2$Sample), 3]


all_genes <- ggplot(long_data, aes(x=Time, y=Count+1, group=Gene, color=Gene))+
  geom_point()+
  geom_path()+
  theme_minimal()+
  scale_y_log10()+
  theme(legend.position = "none")+
  ggtitle("All Genes")




sig_data <- data.frame(results(dds)[c(2, 6)], Gene=rownames(results(dds)[1]))

sig_data <- sig_data[order(sig_data$padj, decreasing = F),]

top68 <- sig_data[1:61,]
bot <- sig_data[200:250,]
siggy <- subset(sig_data, sig_data$padj<0.01)
#df 4 = 61
#df 5 = 68



plotly::ggplotly(ggplot(subset(long_data, long_data$Gene %in% top20$Gene), aes(x=Time, y=Count+1, group=Gene, color=ID))+
                   geom_point()+
                   geom_path()+
                   theme_minimal()+
                   facet_grid(~Gene, scales = "free")+
                   #scale_y_log10()+
                   theme(legend.position = "none")+
                   ggtitle("Top 20")
)









#gene.labs <- paste(top68$Gene, "\n", paste(substr(log10(top68$padj), 1,4)), sep=' ')

#gene.labs <- paste(top68$Gene, "\n", top68$padj, sep = ' ')

ddd <- subset(long_data, long_data$Gene %in% top68$Gene)

ddd$padj <- top68[match(ddd$Gene, top68$Gene), 2]

ddd$gg <- factor(ddd$Gene, levels = unique(ddd[order(ddd$padj),1]))




gene.labs <- paste(factor(top68$gg, levels = top68[order(top68$padj),]$gg), "\n", paste(substr(format(top68$padj, scientific = T), 1, 3),
                                                                                        substr(format(top68$padj, scientific = T), nchar(format(top68$padj, scientific = T))-3, nchar(format(top68$padj, scientific = T))), sep=' '),
                   sep=' ')

names(gene.labs) <- factor(top68$gg, levels = top68[order(top68$padj),]$gg)



top68_graph <- ggplot(ddd, aes(x=Time, y=Count, group=gg, color=ID))+
  geom_point()+
  geom_path()+
  theme_minimal()+
  facet_wrap(~gg, labeller=labeller(gg=gene.labs), scale="free")+
  theme(legend.position = "none")+
  ggtitle("Top 61")



ggsave("top68_graph.pdf", top68_graph, width=25, height = 19)









ddd <- subset(long_data, long_data$Gene %in% bot$Gene)

ddd$padj <- bot[match(ddd$Gene, bot$Gene), 2]

ddd$gg <- factor(ddd$Gene, levels = unique(ddd[order(ddd$padj),1]))




gene.labs <- paste(factor(bot$Gene, levels = bot[order(bot$padj),]$Gene), "\n", paste(substr(format(bot$padj, scientific = T), 1, 3),
                                                                                        substr(format(bot$padj, scientific = T), nchar(format(bot$padj, scientific = T))-3, nchar(format(bot$padj, scientific = T))), sep=' '),
                   sep=' ')

names(gene.labs) <- factor(bot$Gene, levels = bot[order(bot$padj),]$Gene)




bot_graph <- ggplot(ddd, aes(x=Time, y=Count, group=gg, color=ID))+
  geom_point()+
  geom_path()+
  theme_minimal()+
  facet_wrap(~gg, labeller=labeller(gg=gene.labs), scale="free")+
  theme(legend.position = "none")+
  ggtitle("Bottom 50")



ggsave("bottom50_graph.pdf", bot_graph, width=25, height = 19)














play <- data.frame(lsSim$matObservedCounts[1,], dfAnnotationDESeq2$Time, matTimeSplineBasis)


colnames(play)[1:2] <- c("Gene", "Time")


lm0 <- lm(Gene~1, data=play)
lm1 <- lm(Gene~spline1, data=play)
lm2 <- lm(Gene~spline1 + spline2, data=play )
lm3 <- lm(Gene~spline1 + spline2 + spline3, data=play)
lm4 <- lm(Gene~spline1 + spline2 + spline3 + spline4, data=play)
lm5 <- lm(Gene~Time+spline1 + spline2 + spline3 + spline4, data=play)


lm4

#### function that efficiently plots lm results in ggplot

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}







long_splines <- tidyr::gather(dfAnnotationDESeq2, Spline, Y, colnames(dfAnnotationDESeq2)[5:8])

ggplot(long_splines, aes(x=Time, y=Y))+
  geom_point()+
  geom_line()+
  facet_grid(~Spline)
theme_minimal()


grid.arrange(ggplotRegression(lm1), ggplotRegression(lm2), ggplotRegression(lm3), ggplotRegression(lm4))


grid.arrange(plot(x=dfAnnotationDESeq2$Time, y=dfAnnotationDESeq2$spline1), plot(x=dfAnnotationDESeq2$Time, y=dfAnnotationDESeq2$spline2), plot(x=dfAnnotationDESeq2$Time, y=dfAnnotationDESeq2$spline3), plot(x=dfAnnotationDESeq2$Time, y=dfAnnotationDESeq2$spline4))


ggplotRegression(lm0)



