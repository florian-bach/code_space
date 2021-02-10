











d_input <- sce

code_id <- colData(d_input)$cluster_id
cluster_id <- metadata(d_input)$cluster_codes[, 
                                           "meta45"][code_id]

colData(d_input)$code_id <- code_id
colData(d_input)$cluster_id <- cluster_id

metadata <- as.list(c(metadata(d_input), all_ei))

cs_by_s <- split(seq_len(ncol(d_input)), colData(d_input)$sample_id)
cs <- unlist(cs_by_s[as.character(all_ei$sample_id)])
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
  
  
  formula <- createFormula(all_ei, cols_fixed = c("timepoint", "n_infection"), cols_random = "volunteer")
  
  all_t6_contrast <- t(diffcyt::createContrast(c(0,0,0,1,0)))
  
  
  p_vals <- rep(NA, length(cluster_id))
 
  
   for (i in seq_along(cluster_id)) {
    p_vals[i] <- tryCatch({
      y <- counts[i, ]/n_cells_smp
      data_i <- cbind(y, n_cells_smp, formula$data)
      # fit <- lme4::glmer(y~timepoint+(1|volunteer), data = data_i, family = binomial,
      #              weights = round(n_cells_smp/100))
      fit <- MASS::glmmPQL(y~timepoint+n_infection, random= ~ 1 | volunteer, data = data_i, family = binomial, weights=n_cells_smp/1000)
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

