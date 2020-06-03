catalyst_z_normalize <- function(es, th=2.5) {
  es_n <- apply(es, 1, function(x) {
    sd <- stats::sd(x, na.rm=TRUE)
    x <- x - mean(x, na.rm=TRUE)
    if (sd != 0) x <- x / sd
    x[x >  th] <-  th
    x[x < -th] <- -th
    return(x)
  })
  return(t(es_n))
}




##sandbox
catalyst_z_normalize <- function(es, th=2) {
  es_n <- apply(matrix(t6_map_data), 1, function(x) {
    sd <- stats::sd(x, na.rm=TRUE)
    x <- x - mean(x, na.rm=TRUE)
    if (sd != 0) x <- x / sd
    x[x >  th] <-  th
    x[x < -th] <- -th
    return(x)
  })
  return(t(es_n))
}

##sandbox



catalyst_anno_factors <- function(x, ids, which, type = c("row", "column")) {
  type <- match.arg(type)
  # get non-numeric cell metadata variables
  cd <- colData(x)
  df <- data.frame(cd, check.names = FALSE)
  df <- select_if(df, ~!is.numeric(.))
  df <- mutate_all(df, ~droplevels(factor(.x)))
  
  # store sample matching
  m <- match(ids, df$sample_id)
  
  # get number of matches per variable
  ns <- split(df, df$sample_id) %>% 
    lapply(mutate_all, droplevels) %>% 
    lapply(summarize_all, nlevels) %>% 
    do.call(what = "rbind")
  
  # keep only uniquely mapable factors included in 'which'
  keep <- names(which(colMeans(ns) == 1))
  keep <- setdiff(keep, c("sample_id", "cluster_id"))
  if (is.character(which))
    keep <- intersect(keep, which)
  if (length(keep) == 0) return(NULL)
  df <- df[m, keep, drop = FALSE]
  
  # get list of colors for each annotation
  lvls <- lapply(as.list(df), levels)
  nlvls <- vapply(lvls, length, numeric(1))
  pal <- brewer.pal(8, "Set3")[-2]
  if (any(nlvls > length(pal)))
    pal <- colorRampPalette(pal)(max(nlvls))
  names(is) <- is <- colnames(df)
  cols <- lapply(is, function(i) {
    u <- pal[seq_len(nlvls[i])]
    names(u) <- lvls[[i]]; u})
  
  
  HeatmapAnnotation(which = type, df = df, 
                    col = cols, gp = gpar(col = "white"))
}


cls <- list(c())
names(cls) <- names(is)
# x = sce
# y = diffcyt_result


top_anno <- catalyst_anno_factors(merged_daf, levels(merged_daf$sample_id), TRUE, "column")

flo_berlin <- berlin[c(4:15)]


ComplexHeatmap::Heatmap(
  matrix = y, 
  #'name = paste0("normalized\n"[normalize], "frequency"),
  col = berlin[c(4:15)],
  na_col = "lightgrey", 
  rect_gp = gpar(col = "white")#',  
  #cluster_rows = FALSE,
  #cluster_columns = FALSE,
  #row_names_side = "left",
  #top_annotation = top_anno,
  #right_annotation = right_anno,
)
