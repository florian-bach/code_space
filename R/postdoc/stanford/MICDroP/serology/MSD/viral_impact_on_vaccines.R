analysis_pairs <- list(
  ## Q1: does seroconversion 8→24 affect vaccine IgG at 24 weeks?
  list(
    label      = "sc_8_24_on_vax_24",
    tp         = "24 weeks",
    sc_levels  = c("nothing", "converts 8 to 24")
  ),
  
  ## Q2: does seroconversion 8→24 affect vaccine IgG at 52 weeks?
  list(
    label      = "sc_8_24_on_vax_52",
    tp         = "52 weeks",
    sc_levels  = c("nothing", "converts 8 to 24")
  ),
  
  # Q3: does seroconversion 24→52 affect vaccine IgG at 52 weeks?
  list(
    label      = "sc_24_52_on_vax_52",
    tp         = "52 weeks",
    sc_levels  = c("nothing", "converts 24 to 52")
  )
)
safe_lm <- function(dat) {
  tryCatch(
    lm(log2_titer ~ virus_seroconversion, data = dat),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
}

extract_contrasts <- function(mod) {
  if (!inherits(mod, "lm")) return(NULL)
  tryCatch(
    emmeans(mod, ~ virus_seroconversion) %>%
      contrast("revpairwise") %>%
      summary(infer = TRUE),
    error = function(e) NULL
  )
}


# linear regression ####
regression_results_raw <- map_dfr(analysis_pairs, function(ap) {
    
    wide_long_msd_cross %>%
      filter(
        timepoint            == ap$tp,
        virus_seroconversion %in% ap$sc_levels
      ) %>%
      mutate(
        virus_seroconversion = factor(
          virus_seroconversion,
          levels = ap$sc_levels   # "nothing" is always reference
        )
      ) %>%
      group_by(virus_antigen, vaccine_antigen) %>%
      nest() %>%
      mutate(
        n_ids     = map_int(data, ~ n_distinct(.x$id)),
        n_sc      = map_int(data, ~ n_distinct(.x$virus_seroconversion)),
        model     = map(data, safe_lm),
        converged = map_lgl(model, ~ inherits(.x, "lm")),
        cont      = map(model, extract_contrasts)
      ) %>%
      filter(converged, n_sc == 2, !map_lgl(cont, is.null)) %>%
      select(virus_antigen, vaccine_antigen, n_ids, cont) %>%
      unnest(cont) %>%
      mutate(analysis = ap$label)
  })
  


## ---- Multiple testing correction ----------------------------

regression_results <- regression_results_raw %>%
  group_by(analysis, virus_antigen, vaccine_antigen) %>%
  ungroup() %>%
  group_by(analysis, vaccine_antigen) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

regression_results%>%
  filter(padj<0.25)

regression_sig_pairs <- regression_results %>%
  filter(padj < 0.1)

## ============================================================
## Visualise significant seroconversion → vaccine IgG pairs
## One panel per significant virus × vaccine × analysis combo
## ============================================================




make_conversion_plot <- function(analyses, tp) {
  wide_long_msd_cross %>%
    filter(timepoint == tp) %>%
    inner_join(
      regression_sig_pairs %>%
        filter(analysis %in% analyses) %>%
        distinct(virus_antigen, vaccine_antigen, analysis),
      by = c("virus_antigen", "vaccine_antigen")
    ) %>%
    filter(virus_seroconversion %in% c("converts 8 to 24", "converts 24 to 52", "nothing")) %>%
    mutate(wrapped_virus_serconversion=str_wrap(virus_seroconversion, width = 9))%>%
    mutate(
      wrapped_virus_serconversion = factor(
        wrapped_virus_serconversion,
        levels = c("converts\n8 to 24", "converts\n24 to 52", "nothing"))
    ) %>%
    ggplot(aes(x = wrapped_virus_serconversion, y = vaccine_titer, fill = virus_seroconversion)) +
    geom_violin() +
    geom_boxplot(outliers = FALSE, width = 0.3) +
    ggpubr::stat_compare_means(
      vjust       = 0.2,
      comparisons = list(
        c("converts\n8 to 24", "converts\n24 to 52"),
        c("nothing",          "converts\n24 to 52"),
        c("nothing",          "converts\n8 to 24")
      )
    ) +
    scale_y_log10()+
    facet_wrap(~ vaccine_antigen + virus_antigen, scales = "free_y") +
    ylab(paste("log2 IgG titer at ", tp, sep="")) +
    scale_fill_manual(values = seroconversion_cols) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), legend.position = "none")
}

## Plot 1: significant pairs where outcome is 24-week titer
plot_24wk <- make_conversion_plot(
  analyses = "sc_8_24_on_vax_24",
  tp       = "24 weeks"
)

## Plot 2: significant pairs where outcome is 52-week titer
plot_52wk <- make_conversion_plot(
  analyses = c("sc_8_24_on_vax_52", "sc_24_52_on_vax_52"),
  tp       = "52 weeks"
)

plot_52wk2 <- make_conversion_plot(
  analyses = "sc_8_24_on_vax_24",
  tp       = "52 weeks"
)

combo_plot <- plot_24wk / plot_52wk2

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/figures_for_paper/revised_viral_seroconversion_figure.png", width=8, height=7, dpi=444)



personal_seroconverions <- wide_long_msd_cross %>%
  distinct(id, virus_antigen, virus_seroconversion) %>%
  group_by(id)%>%
  summarise(overall_seroconverions = sum(virus_seroconversion != "nothing"),
            early_seroconverions = sum(virus_seroconversion == "converts 8 to 24"),
            late_seroconverions = sum(virus_seroconversion == "converts 24 to 52"))%>%
  arrange(desc(overall_seroconverions))


long_msd%>%
  left_join(., personal_seroconverions, by="id")%>%
  filter(antigen %in% c(vaccines, vaccines_iga[3:5]),
         timepoint=="52 weeks")%>%
  ggplot(., aes(x=overall_seroconverions, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman")+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~antigen)

long_msd%>%
  left_join(., personal_seroconverions, by="id")%>%
  filter(antigen %in% c(vaccines, vaccines_iga[3:5]),
         timepoint=="52 weeks")%>%
  ggplot(., aes(x=early_seroconverions, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman")+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~antigen)

long_msd%>%
  left_join(., personal_seroconverions, by="id")%>%
  filter(antigen %in% c(vaccines, vaccines_iga[3:5]),
         timepoint=="52 weeks")%>%
  ggplot(., aes(x=late_seroconverions, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman")+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~antigen)


long_msd%>%
  left_join(., personal_seroconverions, by="id")%>%
  filter(antigen %in% c(vaccines, vaccines_iga[3:5]),
         timepoint=="24 weeks")%>%
  ggplot(., aes(x=early_seroconverions, y=titer))+
  geom_smooth(method="lm")+
  geom_point()+
  ggpubr::stat_cor(method="spearman")+
  scale_y_log10()+
  theme_minimal()+
  facet_wrap(~antigen)

 
long_msd%>%
  left_join(., personal_seroconverions, by="id")%>%
  filter(antigen %in% c('Tetanus', "Polio", "Diphtheria"),
         timepoint=="24 weeks")%>%
  ggplot(., aes(x=early_seroconverions, y=titer))+
  geom_smooth(method="lm")+
  geom_point(position = position_jitter(width = 0.2))+
  ggpubr::stat_cor(method="spearman", color="red", label.y = 5.1)+
  scale_y_log10()+
  theme_minimal(base_size = 15)+
  facet_wrap(~antigen)



# sandbox ####


## wilcoxon alternative ####

extract_wilcox <- function(dat) {
  tryCatch({
    groups <- split(dat$log2_titer, dat$virus_seroconversion)
    if (length(groups) != 2) return(NULL)
    
    test <- wilcox.test(groups[[1]], groups[[2]], exact = FALSE)
    
    data.frame(
      estimate = median(groups[[2]]) - median(groups[[1]]),  # median difference: seroconverted - never
      p.value  = test$p.value,
      contrast = paste(names(groups)[2], "-", names(groups)[1])
    )
  }, error = function(e) NULL)
}

wilcox_results_raw <- map_dfr(analysis_pairs, function(ap) {
  
  wide_long_msd_cross %>%
    filter(
      timepoint            == ap$tp,
      virus_seroconversion %in% ap$sc_levels
    ) %>%
    mutate(
      virus_seroconversion = factor(
        virus_seroconversion,
        levels = ap$sc_levels
      )
    ) %>%
    group_by(virus_antigen, vaccine_antigen) %>%
    nest() %>%
    mutate(
      n_ids  = map_int(data, ~ n_distinct(.x$id)),
      n_sc   = map_int(data, ~ n_distinct(.x$virus_seroconversion)),
      cont   = map(data, extract_wilcox)
    ) %>%
    filter(n_sc == 2, !map_lgl(cont, is.null)) %>%
    select(virus_antigen, vaccine_antigen, n_ids, cont) %>%
    unnest(cont) %>%
    mutate(analysis = ap$label)
})


wilcox_results <- wilcox_results_raw %>%
  group_by(analysis, virus_antigen, vaccine_antigen) %>%
  ungroup() %>%
  group_by(analysis, virus_antigen) %>%
  mutate(padj1 = p.adjust(p.value, method = "BH")) %>%
  group_by(analysis, vaccine_antigen) %>%
  mutate(padj2 = p.adjust(p.value, method = "BH")) %>%
  ungroup()

wilcox_sig_pairs <- wilcox_results%>%
  filter(padj2<0.1)


## wilcoxon visualisation ####

## ============================================================
## Visualise significant seroconversion → vaccine IgG pairs
## Wilcoxon version — estimate is median difference
## ============================================================


results <- map_dfr(analysis_pairs, function(ap) {
  
  wide_long_msd_cross %>%
    filter(
      timepoint            == ap$tp,
      virus_seroconversion %in% ap$sc_levels
    ) %>%
    mutate(virus_seroconversion = factor(
      virus_seroconversion,
      levels = ap$sc_levels   # first level = reference ("nothing")
    )) %>%
    group_by(virus_antigen, vaccine_antigen) %>%
    nest() %>%
    mutate(
      n_ids     = map_int(data, ~ n_distinct(.x$id)),
      model     = map(data, ~ tryCatch(
        lm(log2_titer ~ virus_seroconversion, data = .x),
        error   = function(e) NULL,
        warning = function(w) NULL
      )),
      converged = map_lgl(model, ~ inherits(.x, "lm")),
      emm       = map(model, ~ if (inherits(.x, "lm"))
        emmeans(.x, ~ virus_seroconversion) else NULL),
      cont      = map(emm, ~ if (!is.null(.x))
        summary(contrast(.x, "revpairwise"), infer = TRUE) else NULL)
    )%>%
    filter(converged, !map_lgl(cont, is.null)) %>%
    select(virus_antigen, vaccine_antigen, n_ids, cont) %>%
    unnest(cont) %>%
    mutate(analysis = ap$label)
})

## Two-stage FDR: Holm within pair, BH across pairs per analysis
results <- results %>%
  group_by(analysis, virus_antigen, vaccine_antigen) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
  ungroup() %>%
  group_by(analysis, contrast) %>%
  mutate(padj = p.adjust(p_holm, method = "BH")) %>%
  ungroup()
