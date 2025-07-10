library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(plotly)
library(viridis)

data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv", header=T, stringsAsFactors = F)

data <- read.csv("/home/flobuntu/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv", header=T, stringsAsFactors = F)

# data<- na.omit(data)
# data$timepoint <- gsub("D", "DoD", data$timepoint, fixed=T)
# data$timepoint <- gsub("DoDoDoD+6", "T+6", data$timepoint, fixed=T)
# data$timepoint <- gsub("DoD+6", "T+6", data$timepoint, fixed=T)
# 
# write.csv(data, "/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv")

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:ncol(data)])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)

long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")


ggplot(long_data, aes(x=factor(long_data$timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point()+
  geom_line(aes(group=long_data$Volunteer))+
  facet_wrap(~ Analyte, scales = "free")+
  theme_bw()+
  scale_color_manual(values=my_paired_palette)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"))

#viva_data <- filter(long_data, Analyte %in% c("CXCL10", "IL12p70", "IL10", "TNFRII", "IL1RA", "ICAM1"))


pca_theme <- theme(strip.background = element_blank(),
                   legend.position = "none",
                   strip.text = element_text(size=20, face = "bold"),
                   axis.text.x = element_text(angle = 60, hjust = 1, size=14),
                   axis.text.y = element_text(size=16),
                   axis.title.y = element_text(size=20),
                   axis.title.x = element_blank())

viva_data <- long_data
(viva_plot <- ggplot(viva_data, aes(x=factor(viva_data$timepoint, levels=c("C-1", "DoD", "T+6", "C+45")), y=Concentration, color=Volunteer), group=Volunteer)+
  geom_point(size=2.5)+
  geom_line(aes(group=viva_data$Volunteer), size=2)+
  scale_y_log10()+
  facet_wrap(~ Analyte, scales = "free")+
  theme_bw()+
  ylab("Plasma Concentration in pg / mL")+
  scale_color_manual(values=my_paired_palette)+
  theme_bw()+
  pca_theme)
)

ggsave("/Users/s1249052/PhD/plasma/vac69a/log_all_hail_legendplex.png", viva_plot, width=35.5, height = 20)
####     pca plot




data2 <- spread(long_data, Analyte, Concentration)
split_data <- split(data2, data2$Volunteer)

my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
names(my_paired_palette) <- names(list_of_pcas)


transformed_list_of_analytes <- lapply(
  split_data, function(x){scales::rescale(as.matrix(x[,3:ncol(x)], to=c(1,100)))
    })

list_of_pcas <- lapply(transformed_list_of_analytes, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
  cbind(x[,1:2], pca$x)
})


list_of_perc <- lapply(transformed_list_of_analytes, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
})



for(i in 1:length(list_of_pcas)){
  
  p <- NULL
  p <- list_of_pcas[[i]]
  q <-list_of_perc[[i]]
  
  assign(paste(names(list_of_pcas[i]), "_pca_plot", sep=''),
  ggplot(p, aes(x=PC1, y=PC2, colour=Volunteer))+
  scale_color_manual(values=my_paired_palette)+
  geom_text(label=p$timepoint, size=6, fontface="bold")+
  ggtitle(paste(names(list_of_pcas[i])))+
  xlab(paste("PC1 ", data.frame(summary(q)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(q)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=25),
        #axis.text = element_blank(),
        plot.title = element_text(size=16, hjust=0.5)
        )
  )
  }


gridExtra::grid.arrange(v002_pca_plot, v003_pca_plot, v005_pca_plot, v006_pca_plot, v007_pca_plot, v009_pca_plot)




list_of_pcas2 <- lapply(split_data, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
  cbind(x[,1:2], pca$x)
})


list_of_perc2 <- lapply(split_data, function(x){
  pca <- prcomp(x[,3:ncol(x)], center = T)
})

top_hits2 <- lapply(list_of_perc2, function(x){
  head(
    x$rotation[order((x$rotation[,1]), decreasing=T),],
    n=10)
})





top_hits <- lapply(list_of_perc, function(x){
  head(
    x$rotation[order(x$rotation[,1], decreasing=T),],
    n=15)
})

top_hits

  head(
    list_of_perc[1]$v002$rotation[order
      (abs(list_of_perc[1]$v002$rotation[,1]), decreasing=T),],
    n=10)

#####     handmade

pca_viva <-ggplot()+
  # geom_point(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer), size=4)+
  # ggrepel::geom_label_repel(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), point.padding = 0.2, size=4)+
  geom_text(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=7, fontface="bold")+
  # geom_label(aes(x=plasma_pca$x[,1], y=plasma_pca$x[,2], color=data$Volunteer, label=data$timepoint), size=4)+
  xlab(paste("PC1 ", summary(plasma_pca)$importance[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", summary(plasma_pca)$importance[2,2]*100, "%", sep = ""))+
  scale_shape_manual(values=c(16:19))+
  scale_color_manual(values=my_paired_palette)+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=25))
   
ggsave("/Users/s1249052/PhD/plasma/vac69a/pca_viva.png", pca_viva, width=12, height = 10)

#####      correlation matrices      ######

# read in data
data <- read.csv("/Users/s1249052/PhD/plasma/vac69a/Vivax_plasma_analytes2_no_inequalities.csv", stringsAsFactors = F)

long_data <- gather(data, Analyte, Concentration, colnames(data)[3:ncol(data)])

long_data$Analyte <- substr(long_data$Analyte, 1, nchar(long_data$Analyte)-7)
long_data$Analyte <- gsub(".", "", long_data$Analyte, fixed = T)
long_data$Concentration <- as.numeric(long_data$Concentration)
long_data <- subset(long_data, Volunteer != "v009")

broad_data <- spread(long_data, Analyte, Concentration)

list_of_tps <- split(broad_data, broad_data$timepoint)


try <- lapply(list_of_tps, function(x){split(x, x$Volunteer)})
try_list <- purrr::flatten(lapply(try, function(x){x[]}))
names(try_list) <- lapply(try_list, function(x){paste(x$Volunteer, x$timepoint, sep='_')})


# apply pearson correlation to analyte part of matrix
list_of_pcorr <- lapply(try, function(x){cor(x[,3:ncol(x)], method = "p")})
#convert to distance matrix
list_of_pdist <- lapply(list_of_pcorr, function(x){dist(x, method="e", diag=F, upper=F, p=2)})
#do hierarchical clustering 
list_of_pclust <- lapply(list_of_pdist, function(x){hclust(x)})
#add column of rownames to create full matrix with gather later
list_of_pcorr <- lapply(list_of_pcorr, function(x){cbind(x, Analyte_X=rownames(x))})
#make data frame, make long format
list_of_long_p <- lapply(list_of_pcorr, function(x){data.frame(x[], stringsAsFactors = F)})
list_of_long_p <- lapply(list_of_pcorr, function(x){tidyr::gather(as.data.frame(x), Analyte_Y, Ro, colnames(x)[1:ncol(x)-1])})


# matrix theme
corr_matrix_theme <-
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust=0.5))

# make a plot for each timepoint; change option of variable(specific order) to choose
for(i in 1:length(list_of_long_p)){
  
  plot_data <- data.frame(list_of_long_p[i])
  
  #run this if you want the order to change per timepoint
  #specific_order <- list_of_pclust[[i]]$order
  
  #run this if you want the order to be fixed at timepoint x;
  # > names(list_of_pclust)
  # [1] "C-1"  "C+45" "DoD"  "T+6" 
  
  specific_order <- list_of_pclust[[1]]$order
  
  colnames(plot_data) <- c("Analyte_X", "Analyte_Y", "Ro")
  plot_data$Ro <- as.numeric(plot_data$Ro)
  
  plt <- ggplot(plot_data, aes_(x=factor(plot_data$Analyte_X, levels = colnames(list_of_pcorr[[i]])[specific_order]) , y=factor(plot_data$Analyte_Y, levels=colnames(list_of_pcorr[[i]])[specific_order])))+
    geom_tile(aes(fill=plot_data$Ro))+
    scale_fill_viridis(option="A")+
    ggtitle(names(list_of_long_p)[i])+
    labs(fill = expression(paste("Pearson ", rho)))+
    corr_matrix_theme
  
  ggsave(paste("/Users/s1249052/PhD/cytof/vac69a/figures_for_paper/plasma_matrices/", names(list_of_long_p)[i], ".png", sep=''), plt, width=7, height=5)
}






library(lme4)

long_data$timepoint <- factor(long_data$timepoint, levels=c("Baseline", "DoD", "T+6", "C+45"))

list_of_analytes <- split(long_data, long_data$Analyte)
plasma_time_models <- lapply(list_of_analytes, function(x) lmer(Concentration~timepoint+(1|Volunteer), data=x))


plasma_time_results <- lapply(plasma_time_models, FUN=function(x)(summary(x)))

pvals <- lapply(plasma_time_results, function(x) x$tTable)
pvals <- lapply(pvals, function(x){data.frame(x)})
pvals <- mapply(cbind, pvals, "model"= names(pvals), SIMPLIFY = F)
pvals <- lapply(pvals, function(x) data.frame(x, "timepoint"=rownames(x)))

results_df <- data.table::rbindlist(pvals)
results_df <- subset(results_df, results_df$timepoint!="(Intercept)")
results_df$p_adjust <-  p.adjust(results_df$p.value)

stable_through_time <- subset(results_df, p_adjust>0.05)
varies_through_time <- subset(results_df, p_adjust<0.05)

paste(varies_through_time$model, varies_through_time$timepoint)





