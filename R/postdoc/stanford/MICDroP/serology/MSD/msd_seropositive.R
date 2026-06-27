
nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")
maternal_treatment_arms <- haven::read_dta("~/Library/CloudStorage/Box-Box/DP+SP study/Databases and preliminary findings/Final database used for analyses/DPSP treatment allocation_FINAL.dta")

vaccines = c("Diptheria",     "Measles" ,      "Mumps",         "Pertussis",     "Polio",
             "Rotavirus" ,    "Rubella",       "Tetanus", "Pneumo.1.4.14")


nulisa_data <- nulisa_data%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"))

metadata_columns <- c("id", "anyDP", "treatmentarm",  "dob", "date", "ageinwks", "gender_categorical", "mstatus", "qPCRparsdens", "visittype", "fever", "febrile", "rogerson", "GAcomputed", "gi", "SGA", "qPCRdich", "mqPCRparsdens")

epi_data <- nulisa_data%>%
  distinct(sample, total_n_para_12, total_n_malaria_12, total_n_malaria_6, total_n_para_6,
           id, dob, date, ageinwks, gender_categorical, mstatus, qPCRparsdens, visittype, fever, febrile, rogerson, GAcomputed, gi, SGA, qPCRdich, mqPCRparsdens, anyHP)%>%
  mutate(treatmentarm=mic_drop_key$treatmentarm[match(as.numeric(id), mic_drop_key$id)],
         anyDP=if_else(treatmentarm==1, "no", "yes"),
         treatmentarm=case_match(treatmentarm,
                                 1~"Placebo",
                                 2~"DP 1 year",
                                 3~"DP 2 years"),
         mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))

msd_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/batch_one.csv")

long_msd <- msd_data%>%
  mutate(sample=paste(SubjectID, "_", "tp", TimePt, sep=""))%>%
  mutate(id=SubjectID, timepoint=paste(TimePt, "weeks"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  select(-SubjectID, -TimePt)%>%
  pivot_longer(cols=-c(sample, id, timepoint), names_to = "antigen", values_to = "titer")

antibodies_and_epi <- epi_data%>%
  inner_join(., long_msd, by="sample")

antibodies_and_nulisa <-  nulisa_data%>%
  select(-id, -timepoint)%>%
  inner_join(., long_msd, by="sample")

# defining seroconversion ####
#  between 8 and 24: calculate average fold change from 8 to 24; assume 2 STDvs above mean = serconversion even if decrease; else if increase at all;
#  between 24 and 52: simple 4fold increase
final_cutoff_frame <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/final_cutoff_frame.csv")
colnames(final_cutoff_frame)[2]="cut_off"

wide_long_msd <-  long_msd %>%
  filter(!is.na(timepoint))%>%
  pivot_wider(names_from = "timepoint", values_from = "titer", id_cols = c("id", "antigen"))%>%
  mutate(log2fc_8_24=log2(`24 weeks`/`8 weeks`),
         log2fc_24_52=log2(`52 weeks`/`24 weeks`))%>%
  group_by(antigen)%>%
  mutate("z_fc_8_24"=scale(log2fc_8_24, center = T))%>%
  pivot_longer(cols = c("8 weeks", "24 weeks","52 weeks"), names_to = "timepoint", values_to = "titer")%>%
  pivot_longer(cols = c("log2fc_8_24", "log2fc_24_52"), names_to = "fc_flavour", values_to = "log2fc")%>%
  mutate(conversion=case_when(timepoint=="52 weeks" & fc_flavour=="log2fc_24_52" & log2fc > 2 ~ "converts 24 to 52",
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & z_fc_8_24 > 2 | 
                              timepoint=="24 weeks" & fc_flavour=="log2fc_8_24" & log2fc > 1 ~ "converts 8 to 24",
                              .default="nothing"))%>%
  left_join(., final_cutoff_frame, by="antigen")%>%
  mutate(sero_positive=if_else(titer>=cut_off, 1, 0))

positivity_df <- wide_long_msd%>%
  distinct(id, timepoint, antigen, sero_positive)%>%
  mutate(mom_rx=maternal_treatment_arms$treatmentarm[match(id-10000, maternal_treatment_arms$id)],
         mom_rx=case_match(mom_rx,
                           1~"SP",
                           2~"DP",
                           3~"DPSP"))


antigens <- unique(positivity_df$antigen)
timepoints <- c("8 weeks", "24 weeks", "52 weeks")


nulisa_and_seropos <-  nulisa_data %>%
  left_join(., positivity_df, by=c("id", "timepoint"))%>%
  filter(timepoint %in% timepoints)


sero_pos_purf <- nulisa_and_seropos%>%
  filter(!is.na(antigen))%>%
  group_by(targetName, antigen, timepoint)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(sero_positive~conc,  family = "binomial",  data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(timepoint, targetName)%>%
  mutate(padj=p.adjust(p, method="fdr"))

sig_purf_conversion_nulisa <- sero_pos_purf%>%
  filter(padj<0.1)


for(i in 1:nrow(sig_purf_conversion_nulisa)){
  
  anti=sig_purf_conversion_nulisa$antigen[i]
  targi=sig_purf_conversion_nulisa$targetName[i]
  timi=sig_purf_conversion_nulisa$timepoint[i]
  
  plt = nulisa_and_seroconversion%>%
    filter(antigen==anti,
           targetName==targi,
           timepoint==timi)%>%
    mutate(sero_positive=case_when(sero_positive==0~"negative",
                                   sero_positive==1~"positive"))%>%
    ggplot(., aes(x=timepoint, y=conc, fill=factor(sero_positive)))+
    geom_boxplot(outliers = F)+
    facet_wrap(~targetName, scales="free")+
    ggtitle(paste0(anti, " at ", timi)) +
    # scale_fill_manual(values = c("orange", "red", "black"))+
    theme_minimal()+
    theme(legend.title = element_blank())
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/seropositivity_nulisa_cor_",anti,"_",targi, ".png"),plt, width = 4, height = 4, dpi = 444, bg="white", limitsize = F)
  
}


# simple table ####
wide_long_msd%>%
  distinct(id, timepoint, antigen, sero_positive)%>%
  group_by(antigen, timepoint)%>%
  summarise(prevalence=mean(sero_positive))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")), y=prevalence, color=antigen,  group=antigen))+
  geom_line()+
  geom_point()+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(~antigen)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank())

# neurocognitive_data

neuro_cog <- read.csv("~/postdoc/stanford/clinical_data/MICDROP/neurocognitive/NCT_infants_share032025.csv")
neuro_cog_edit <- neuro_cog%>%
  mutate(id=subjid, date=as.Date(vdate))%>%
  select(-gender, -GAcomputed, -SGA, -date, -dob)%>%
  pivot_longer(cols = ends_with("composite"), names_to = "composite_kind", values_to = "composite_score")

seropos_plus_neuro <- left_join(positivity_df, neuro_cog_edit, by="id")


seropos_neuro_cog_purrr <- seropos_plus_neuro%>%
  filter(!is.na(composite_kind))%>%
  group_by(antigen, timepoint, composite_kind)%>%
  nest()%>%
  mutate(model = purrr::map(data, ~glm(sero_positive~composite_score,  family = "binomial",  data=.)))%>%
  mutate(model_summary=purrr::map(model, ~summary(.)))%>%
  mutate(p=purrr::map(model_summary, ~coef(.)[8]))%>%
  group_by(timepoint, antigen)%>%
  mutate(padj=p.adjust(p, method="fdr"))

seropos_neuro_cog_purrr_sigs <- seropos_neuro_cog_purrr %>%
  filter(padj<0.05)


for(i in 1:nrow(seropos_neuro_cog_purrr_sigs)){
  
  anti=seropos_neuro_cog_purrr_sigs$antigen[i]
  compo=seropos_neuro_cog_purrr_sigs$composite_kind[i]
  timi=seropos_neuro_cog_purrr_sigs$timepoint[i]
  
  plt = seropos_plus_neuro%>%
    filter(antigen==anti,
           composite_kind==compo,
           timepoint==timi)%>%
    mutate(sero_positive=case_when(sero_positive==0~"negative",
                                   sero_positive==1~"positive"))%>%
    ggplot(., aes(x=timepoint, y=composite_score, fill=factor(sero_positive)))+
    geom_boxplot(outliers = F)+
    facet_wrap(~composite_kind, scales="free")+
    ggtitle(paste0(anti, " at ", timi)) +
    ggpubr::stat_compare_means()+
    # scale_fill_manual(values = c("orange", "red", "black"))+
    theme_minimal()+
    theme(legend.title = element_blank())
  
  ggsave(paste0("~/postdoc/stanford/plasma_analytes/MICDROP/MSD/figures/seropositivity_neurocog_",anti,"_",anti, ".png"),plt, width = 4, height = 4, dpi = 444, bg="white", limitsize = F)
  
}





# safe coerce to character (works if column is factor, character, or a list-column)
safe_as_char <- function(x) {
  if (is.factor(x)) return(as.character(x))
  if (is.list(x)) {
    # if an element has length != 1 we return NA (marks malformed rows)
    vapply(x, function(el) {
      if (length(el) == 1) as.character(el) else NA_character_
    }, FUN.VALUE = character(1))
  } else {
    as.character(x)
  }
}

# full set of levels from the whole dataset (ensures consistent table dims)


conv_levels <- levels(factor(positivity_df$sero_positive))
mom_levels  <- levels(factor(positivity_df$mom_rx))

# process groups robustly using group_split + map_dfr
seropos_momrx_chisq <- positivity_df %>%
  group_by(antigen, timepoint) %>%
  group_split() %>%
  map_dfr(function(df) {
    antigen_name <- df$antigen[1]
    timepoint_name <- df$timepoint[1]
    
    # coerce safely to plain character vectors
    conv <- safe_as_char(df$sero_positive)
    mom  <- safe_as_char(df$mom_rx)
    
    # drop rows where coercion produced NA (malformed list-elements)
    keep <- !is.na(conv) & !is.na(mom)
    conv <- conv[keep]; mom <- mom[keep]
    
    # build contingency table with the global levels (zeros if absent)
    tbl <- table(factor(conv, levels = conv_levels),
                 factor(mom,  levels = mom_levels))
    total <- sum(tbl)
    
    # invalid / trivial tables
    if (total == 0 || nrow(tbl) < 2 || ncol(tbl) < 2) {
      return(tibble(
        antigen = antigen_name,
        statistic = NA_real_,
        parameter = NA_real_,
        p.value = NA_real_,
        timepoint=timepoint_name,
        method = "not_enough_levels_or_data",
        expected_min = NA_real_,
        total = total,
        contingency = list(tbl)
      ))
    }
    
    # expected counts
    exp_tbl <- outer(rowSums(tbl), colSums(tbl), "*") / total
    min_exp <- min(exp_tbl, na.rm = TRUE)
    
    # choose test: standard chi-square if expected >=5, else fisher for 2x2, else simulate
    test_res <- tryCatch({
      if (!is.na(min_exp) && min_exp >= 5) {
        chisq.test(tbl)
      } else if (nrow(tbl) == 2 && ncol(tbl) == 2) {
        fisher.test(tbl)
      } else {
        chisq.test(tbl, simulate.p.value = TRUE, B = 2000)
      }
    }, error = function(e) e)
    
    if (inherits(test_res, "error")) {
      tibble(
        antigen = antigen_name,
        statistic = NA_real_,
        parameter = NA_real_,
        p.value = NA_real_,
        timepoint=timepoint_name,
        method = paste0("test_error: ", test_res$message),
        expected_min = min_exp,
        total = total,
        contingency = list(tbl)
      )
    } else {
      td <- broom::tidy(test_res)
      tibble(
        antigen = antigen_name,
        timepoint=timepoint_name,
        statistic = td$statistic,
        parameter = td$parameter,
        p.value = td$p.value,
        method = td$method,
        expected_min = min_exp,
        total = total,
        contingency = list(tbl)
      )
    }
  }) %>%
  group_by(antigen)%>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

seropos_momrx_chisq%>%
  arrange(p.value)


positivity_df%>%
  filter(antigen %in% c("PIV.4", "PIV.1"), timepoint=="8 weeks")%>%
  group_by(antigen, mom_rx)%>%
  summarise(seroprevalence=mean(sero_positive))%>%
  ggplot(., aes(x="", y=seroprevalence+0.001, fill = mom_rx))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values=colorspace::sequential_hcl(3, palette = "Purple Orange"))+
  facet_wrap(~antigen)+
  xlab("")+
  ylab('seroprevalence')+
  theme_minimal()
