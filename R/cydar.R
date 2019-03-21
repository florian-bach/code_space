library(ncdfFlow)
library(edgeR)
library(cydar)
library(flowCore)


######         make some mock data


ncells <- 20000
nda <- 200
nmarkers <- 31
down.pos <- 1.8
up.pos <- 1.2
conditions <- rep(c("A", "B"), each=3)
combined <- rbind(matrix(rnorm(ncells*nmarkers, 1.5, 0.6), ncol=nmarkers),
                  matrix(rnorm(nda*nmarkers, down.pos, 0.3), ncol=nmarkers),
                  matrix(rnorm(nda*nmarkers, up.pos, 0.3), ncol=nmarkers))
combined[,31] <- rnorm(nrow(combined), 1, 0.5) # last marker is a QC marker. 
combined <- 10^combined # raw intensity values                  
sample.id <- c(sample(length(conditions), ncells, replace=TRUE), 
               sample(which(conditions=="A"), nda, replace=TRUE), 
               sample(which(conditions=="B"), nda, replace=TRUE))
colnames(combined) <- paste0("Marker", seq_len(nmarkers))


#We use this to construct a ncdfFlowSet for our downstream analysis.

collected.exprs <- list()
for (i in seq_along(conditions)) {
  stuff <- list(combined[sample.id==i,,drop=FALSE])
  names(stuff) <- paste0("Sample", i)
  collected.exprs[[i]] <- poolCells(stuff)
}
names(collected.exprs) <- paste0("Sample", seq_along(conditions))
collected.exprs <- ncdfFlowSet(as(collected.exprs, "flowSet"))


# The intensities need to be transformed and gated prior to further analysis. We first pool all cells together into a single flowFrame, which will be used
# for construction of the transformation and gating functions for all samples. This avoids spurious differences from using sample-specific functions.


pool.ff <- poolCells(collected.exprs)

# We use the estimateLogicle method from the flowCore package to obtain a transformation function, and apply it to pool.ff. This performs a biexponential
# transformation with parameters estimated for optimal display.

trans <- estimateLogicle(pool.ff, colnames(pool.ff))
proc.ff <- transform(pool.ff, trans)


gate.31 <- outlierGate(proc.ff, "Marker31", type="upper")
filter.31 <- filter(proc.ff, gate.31)
summary(filter.31@subSet)




##########                  hyperspheres

processed.exprs <- transform(collected.exprs, trans)
processed.exprs <- Subset(processed.exprs, gate.31)
processed.exprs <- processed.exprs[,1:30]

cd <- prepareCellData(processed.exprs)

cd <- countCells(cd, tol=0.5)


###############            The output is another CyData object with extra information added to various fields. In particular, the reported count matrix
###############            contains the set of counts for each hypersphere (row) from each sample (column).

head(assay(cd))

head(intensities(cd))

#######           EdgeR           ###########

y <- DGEList(assay(cd), lib.size=cd$totals)

# assay(CD) = matrix with set of counts for each hypersphere (row) from each sample (column)
#cd$totals = string of numbers indicating number of cells in sample

# translated, the assay(CD) object could be a matrix of cluster percentages (rows) per person (columns)

# to construct the assay(CD object), one would concatenate the dataframes of cluster abundances from all volunteers

keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
cd <- cd[keep,]
y <- y[keep,]

design <- model.matrix(~factor(conditions))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)

######        FDR        ####

qvals <- spatialFDR(intensities(cd), res$table$PValue)

is.sig <- qvals <= 0.05
summary(is.sig)


sig.coords <- intensities(cd)[is.sig,]
sig.res <- res$table[is.sig,]
coords <- prcomp


plotSphereLogFC(coords$x[,1], coords$x[,2], sig.res$logFC)
