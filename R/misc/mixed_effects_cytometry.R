# vac63c ####
library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
#setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))
fcs <- subset(fcs, !grepl(pattern = "307", fcs))
fcs <- subset(fcs, !grepl(pattern = "302", fcs))


vac63_flowset <- flowCore::read.flowSet(fcs)

md <- read.csv("vac63c_metadata.csv", header = T)

md <- subset(md, md$file_name %in% fcs)
md <- md[order(md$timepoint),]


#  
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

sce <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))



# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 7                   
# ggplot2::ggsave("marker_expression.png", p, width=12, height=12)                                      


refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     #"CXCR5",
                     "CD38",
                     #"CD69",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     #"TCRgd",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")

# all_markers <- c("CD45", "CD3", "CD3", "CD14", "CD16", refined_markers)

set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)




meta45_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_50_merge.csv")

d_input <- CATALYST::mergeClusters(sce, k = "meta50", table = meta45_table, id = "flo_merge", overwrite = TRUE)


code_id <- colData(d_input)$cluster_id
cluster_id <- metadata(d_input)$cluster_codes[, 
                                              "flo_merge"][code_id]

colData(d_input)$code_id <- code_id
colData(d_input)$cluster_id <- cluster_id


# vac69a ####


setwd("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

d_input <- vac69a.cytof::read_full(path_to_directory = "/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

code_id <- colData(d_input)$cluster_id
cluster_id <- metadata(d_input)$cluster_codes[, 
                                           "flo_merge"][code_id]

colData(d_input)$code_id <- code_id
colData(d_input)$cluster_id <- cluster_id


# glmmm ####



ei <- metadata(d_input)$experiment_info
#ei$timepoint <- factor(ei$timepoint, levels=c("C10", "Baseline", "DoD", "T6"))

# alt <- ifelse(ei$timepoint=="T6"&ei$volunteer!="V05", TRUE, FALSE)
# ei$alt <- alt
#ei$alt <- ei$alt>40
design <- createDesignMatrix(ei, c("timepoint", "volunteer", "n_infection"))
# 

#design <- model.matrix(~ei$time+ei$volunteer:ei$alt)
# batch_design <- createDesignMatrix(ei, c("timepoint", "t"))

FDR_cutoff <- 0.05


metadata <- as.list(c(metadata(d_input), ei))

cs_by_s <- split(seq_len(ncol(d_input)), colData(d_input)$sample_id)
cs <- unlist(cs_by_s[as.character(ei$sample_id)])
es <- t(assays(d_input)[["exprs"]])[cs, , drop = FALSE]
d_se <- SummarizedExperiment(assays = list(exprs = es), 
                             rowData = colData(d_input)[cs, ], colData = rowData(d_input), 
                             metadata = metadata(d_input))

d_counts <- calcCounts(d_se)

# test da glmm


  
  min_cells = 3
  min_samples = 4
  
  counts <- assays(d_counts)[["counts"]]
  cluster_id <- rowData(d_counts)$cluster_id
  tf <- counts >= min_cells
  ix_keep <- apply(tf, 1, function(r) sum(r) >= min_samples)
  counts <- counts[ix_keep, , drop = FALSE]
  cluster_id <- cluster_id[ix_keep]
  n_cells_smp <- colSums(counts)
  
  
  formula <- createFormula(ei, cols_fixed = c("timepoint"), cols_random = "volunteer")
  
  #no n_infection
  #all_t6_contrast <- t(diffcyt::createContrast(c(0,0,0,1)))
  
  #with n_infection
  all_t6_contrast <- t(diffcyt::createContrast(c(0,0,0,1)))
  
  
  p_vals <- rep(NA, length(cluster_id))
 
  
   for (i in seq_along(cluster_id)) {
    p_vals[i] <- tryCatch({
      y <- counts[i, ]/n_cells_smp
      data_i <- cbind(y, n_cells_smp, formula$data)
      #fit <- lme4::glmer(y~timepoint+(1|volunteer), data = data_i, family = binomial, weights = NULL)
      # fit <- survey::svyglm(y~timepoint+volunteer, design=survey::svydesign(ids=~1,
      #                                                                       weights=~n_cells_smp,
      #                                                                       data=data_i), family=binomial)
      fit <- MASS::glmmPQL(y~timepoint, random= ~ 1 | volunteer, data = data_i, family = binomial, weights=n_cells_smp)
      test <- multcomp::glht(fit, all_t6_contrast)
      summary(test)$test$pvalues
    }, error = function(e) NA)
   }
  
  p_adj <- p.adjust(p_vals, method = "fdr")
  stopifnot(length(p_vals) == length(p_adj))
  out <- data.frame(p_val = p_vals, p_adj = p_adj, stringsAsFactors = FALSE)
  row_data <- as.data.frame(matrix(as.numeric(NA), nrow = nlevels(cluster_id), 
                                   ncol = ncol(out)))
  colnames(row_data) <- colnames(out)
  cluster_id_nm <- as.numeric(cluster_id)
  row_data[cluster_id_nm, ] <- out
  row_data <- cbind(cluster_id = rowData(d_counts)$cluster_id, 
                    row_data)
  res <- d_counts
  rowData(res) <- row_data
 
  first <- row_data
  
  third <- row_data
  
  
  
  
  prim_fold_changes_t6 <- counts[,c(7:9, 16:18, 25:27, 34:36)]
  ter_fold_changes_t6 <- counts[,-c(7:9, 16:18, 25:27, 34:36)]
  
  
  #apply(counts[,grepl("T6", colnames(prim_fold_changes_t6))], 1, mean)
  
  #select each timepoint seperately, calc cluster-wise mean, make fold change, log2 transform
  
  prim_fold_changes_t6 <- log2(apply(prim_fold_changes_t6[,grepl("T6", colnames(prim_fold_changes_t6))], 1, mean) / apply(prim_fold_changes_t6[,grepl("Baseline", colnames(prim_fold_changes_t6))], 1, mean))
  ter_fold_changes_t6 <-  log2(apply(ter_fold_changes_t6[,grepl("T6", colnames(ter_fold_changes_t6))], 1, mean) / apply(ter_fold_changes_t6[,grepl("Baseline", colnames(ter_fold_changes_t6))], 1, mean))
  
  first$logFC <- prim_fold_changes_t6[match(first$cluster_id, names(prim_fold_changes_t6))]
  third$logFC <- ter_fold_changes_t6[match(third$cluster_id, names(ter_fold_changes_t6))]
  
  sig_first <- subset(first, first$p_adj<0.05 & abs(first$logFC)>1)
  #20
  write.csv(x = first, "glmmPQL_first.csv")
  sig_third <- subset(third, third$p_adj<0.05 & abs(third$logFC)>1)
  #6
  write.csv(x = third, "glmmPQL_third.csv")
  
  
  # all(sig_third$cluster_id %in% sig_third$cluster_id)
  # [1] TRUE
  
  table(row_data$p_adj<0.05)
  
  # first
  # FALSE  TRUE 
  # 16    34
  
  # third
  # FALSE  TRUE 
  # 41     9 
  
  gg_counts <- data.frame(counts)
  
  gg_percs <- lapply(1:ncol(counts), function(x) counts[,x] / n_cells_smp[x])
  
  gg_percs <- do.call(cbind, gg_percs)
  
  colnames(gg_percs) <- colnames(counts)
  
  gg_counts <- data.frame(gg_percs)
  
  gg_counts$cluster <- rownames(gg_counts)
  
  gg_counts <- gather(gg_counts, sample_id, count, colnames(gg_counts)[1:36])
  
  gg_counts$volunteer <- substr(gg_counts$sample_id, 1,4)
  
  gg_counts$timepoint <- rep("hello", nrow(gg_counts))
  gg_counts$timepoint <- ifelse(grepl("Baseline", gg_counts$sample_id), "Baseline",   gg_counts$timepoint)
  gg_counts$timepoint <- ifelse(grepl("C45", gg_counts$sample_id), "C45",   gg_counts$timepoint)
  gg_counts$timepoint <- ifelse(grepl("DoD", gg_counts$sample_id), "Diagnosis",   gg_counts$timepoint)
  gg_counts$timepoint <- ifelse(grepl("T6", gg_counts$sample_id), "T6",   gg_counts$timepoint)
  
  gg_counts$n_infection <- "third"
  
  gg_counts$n_infection <- ifelse(gg_counts$volunteer=="v313", "first", gg_counts$n_infection)
  gg_counts$n_infection <- ifelse(gg_counts$volunteer=="v315", "first", gg_counts$n_infection)
  gg_counts$n_infection <- ifelse(gg_counts$volunteer=="v320", "first", gg_counts$n_infection)
  
  
  time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")
  
  
  first_gg_counts <- subset(gg_counts, gg_counts$cluster %in% sig_first$cluster_id)
  third_gg_counts <- subset(gg_counts, gg_counts$cluster %in% sig_third$cluster_id)
  
sig_first_plot <-  ggplot(first_gg_counts, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45")), y=count*100, fill=n_infection))+
    geom_boxplot(aes(fill=n_infection))+
    geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
    facet_wrap(~cluster, scales = "free", ncol = 5, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
    theme_minimal()+
    ylab("% of all CD3+ cells")+
    xlab("Timepoint")+
  
    scale_y_continuous()+
    scale_fill_manual(values = c("first" = time_col[4],
                                 "third" = time_col[1]))
   
ggsave("sig_first_plot.png", sig_first_plot, width=15, height=12) 
  
  
sig_third_plot <-  ggplot(third_gg_counts, aes(x=factor(timepoint, levels=c("Baseline", "Diagnosis", "T6", "C45")), y=count*100, fill=n_infection))+
    geom_boxplot(aes(fill=n_infection))+
    geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
    facet_wrap(~cluster, scales = "free", ncol = 6, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
    theme_minimal()+
    ylab("% of all CD3+ cells")+
    xlab("Timepoint")+
    scale_y_continuous()+
    scale_fill_manual(values = c("first" = time_col[4],
                                 "third" = time_col[1]))

ggsave("sig_third_plot.png", sig_third_plot, width=14, height=7)  
