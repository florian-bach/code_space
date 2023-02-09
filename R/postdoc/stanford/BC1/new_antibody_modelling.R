#combine this 
long_raw_df

# id    timepoint antigen              conc flag             flag_type  
# <fct> <fct>     <chr>               <dbl> <dbl+lbl>        <chr>      
#   1 11001 1         CSP GENOVA       0.0151   NA               No Flag    
# 2 11001 1         EBA140 RIII V JB 0.00175   4 [BelowMinStd] BelowMinStd
# 3 11001 1         EBA175 RIII Vkt  0.00271  NA               No Flag    
# 4 11001 1         EBA181 VIII V    0.000102  4 [BelowMinStd] BelowMinStd
# 5 11001 1         Etramp4 Ag2      0.00217  NA               No Flag    
# 6 11001 1         Etramp5 Ag1 kt   0.000143  4 [BelowMinStd] BelowMinStd

wide_raw_df <- long_raw_df %>%
  select(-flag_type)%>%
  pivot_wider(names_from = c(timepoint, antigen), values_from = c(conc, flag))

# with this
clin_data <- clin_data %>%
  filter(id %in% wide_raw_df$id)


infs <- clin_data %>%
  group_by(id) %>%
  select(id, age, anyinfection, sxifinfected) %>%
  mutate(inf_0_6   = sum(if_else(age<6, anyinfection, 0), na.rm = TRUE),
         inf_6_12  = sum(if_else(age>6  & age<12, anyinfection, 0), na.rm = TRUE),
         inf_12_18 = sum(if_else(age>12 & age<18, anyinfection, 0), na.rm = TRUE),
         inf_12_24 = sum(if_else(age>12 & age<24, anyinfection, 0), na.rm = TRUE),
         inf_0_12  = sum(if_else(age<12, anyinfection, 0), na.rm = TRUE),
         inf_0_24  = sum(if_else(age<24, anyinfection, 0), na.rm = TRUE),
         symp_0_6   = sum(if_else(age<6, sxifinfected, 0), na.rm = TRUE),
         symp_6_12  = sum(if_else(age>6  & age<12, sxifinfected, 0), na.rm = TRUE),
         symp_12_18 = sum(if_else(age>12 & age<18, sxifinfected, 0), na.rm = TRUE),
         symp_12_24 = sum(if_else(age>12 & age<24, sxifinfected, 0), na.rm = TRUE),
         symp_0_12  = sum(if_else(age<12, sxifinfected, 0), na.rm = TRUE),
         symp_0_24  = sum(if_else(age<24, sxifinfected, 0), na.rm = TRUE)
         ) %>%
  select(-anyinfection, -sxifinfected, -age) %>%
  distinct()
  
         
combo_data <- cbind(long_raw_df, infs[match(long_raw_df$id, infs$id),])

combo_data <- combo_data[,-7]


combo_data %>%
  filter(antigen %in% modelable_antigens)%>%
  ggplot(aes(x=inf_0_12, y=conc, fill=factor(inf_0_12)))+
  geom_boxplot()+
  facet_grid(timepoint~antigen, labeller = labeller(antigen = label_wrap_gen(width = 6)))+
  scale_y_log10()+
  xlab("Number of Infections in First Year of Life")+
  ylab("Concentration")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  viridis::scale_fill_viridis(option="B", direction = -1, discrete = TRUE)
  
  