library(immunarch)
library(tidyverse)

vdj <- read.csv("/Volumes/lab_prasj/BIG_Flo/whitney_harrington_microchimerism/data/merged_tcr_contigs.csv")

imm_tbl <- vdj %>%
  transmute(
    Sample = sample,
    Clonotype = raw_clonotype_id,
    CDR3.aa = cdr3,
    CDR3.nt = cdr3_nt,
    V.name = v_gene,
    D.name = d_gene,
    J.name = j_gene,
    Chain = chain
  )

immdata <- list()

immdata$data <- imm_tbl %>%
  split(.$Sample) %>%
  lapply(function(df) df %>% select(-Sample))

immdata$meta <- tibble(
  Sample = names(immdata$data)
)



immdata$data <- lapply(immdata$data, function(df) {
  df %>%
    group_by(Clonotype, CDR3.aa, CDR3.nt, V.name, D.name, J.name, Chain) %>%
    summarise(Clones = n(), .groups = "drop")
})

immdata$data <- lapply(immdata$data, function(df) {
  
  df %>%
    mutate(
      Proportion = Clones / sum(Clones)
    )
})
saveRDS(immdata, "/Volumes/lab_prasj/BIG_Flo/whitney_harrington_microchimerism/data/immunarch_object.rds")

exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Sample"), .meta = immdata$meta)



exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

imm_pr <- repClonality(immdata$data, .method = "clonal.prop")

imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))

imm_rare <- repClonality(immdata$data, .method = "rare")

imm_hom <- repClonality(immdata$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)

vis(imm_top) + vis(imm_top, .by = "Sample", .meta = immdata$meta)

vis(imm_rare) + vis(imm_rare, .by = "Sample", .meta = immdata$meta)

vis(imm_hom) + vis(imm_hom, .by = c("Sample"), .meta = immdata$meta)
