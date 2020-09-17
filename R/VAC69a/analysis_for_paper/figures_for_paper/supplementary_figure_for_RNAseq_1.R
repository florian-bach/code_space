# B vs. T cell counts

      data <- data.table::fread("/home/flobuntu/PhD/clinical_data/vac69a/lymph_counts.csv")
      
      long_data <- tidyr::gather(data, Cell_Type, Million_Cells_Per_mL, colnames(data)[2:3])
      
      ggplot(long_data, aes(x=factor(Timepoint), y=Million_Cells_Per_mL))+
        geom_boxplot(aes(fill=Timepoint))+
        geom_point(aes(shape=Volunteer))+
        facet_wrap(~Cell_Type, scales="free")+
        theme_minimal()
              