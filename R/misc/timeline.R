library(timevis)


groups <- data.frame(
  id=c("Key Deliverables", "CHMI Trial")
  content= c("reVAC63#2", "Dose#1", "Dose#2", "Vac3#1", "Dose#3", "Vac3#2", "Vac12#1", "Vac3#3", "Vac12#2", "Vac12#3")

)

timevis(data = data.frame(
  id      = 1:10,
  content = c("reVAC63#2", "Dose#1", "Dose#2", "Vivax3#1", "Dose#3", "Vivax3#2", "Vivax12#1", "Vivax3#3", "Vivax12#2", "Vivax12#3"),
  start   = c("2018-11-10", "2019-01-10", "2019-06-01", "2019-06-01", "2019-11-01", "2019-11-01", "2019-11-01", "2020-04-14", "2020-04-14", "2020-09-01"),
  #end     = c("2019-01-10","2019-03-10", "2019-08-01", "2019-08-01", "2020-01-01", "2020-01-01", "2020-01-01", "2020-06-14", "2020-06-14", "2020-11-01", "2019-01-10", "2019-06-01"),
  end     = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  group   = c(4,1,1,2,1,2,3,2,3,3)),
  groups  = data.frame(id=1:4, content=c("Dose Escalation (n=6)", "Vivax3 (n=3)", "Vivax12 (n=12)", "P. falciparum (n=16)")),
  showZoom = F, fit=T)


timevis(data, showZoom = F, fit=T)






library(reshape2)
library(ggplot2)
library(directlabels)

tasks <- c("reVAC63#2", "Dose#1", "Dose#2", "Vac3#1", "Dose#3", "Vac3#2", "Vac12#1", "Vac3#3", "Vac12#2", "Vac12#3")

dfr <- data.frame(
  name        = factor(tasks, levels = tasks),
  start.date  = as.Date(c("2018-11-10", "2019-01-10", "2019-06-01", "2019-06-01", "2019-11-01", "2019-11-01", "2019-11-01", "2020-04-14", "2020-04-14", "2020-09-01")),
  end.date    = as.Date(c("2019-01-10","2019-03-10", "2019-08-01", "2019-08-01", "2020-01-01", "2020-01-01", "2020-01-01", "2020-06-14", "2020-06-14", "2020-11-01")),
  color.ful   = factor(c(4,1,1,2,1,2,3,2,3,3))
  )
mdfr <- melt(dfr, measure.vars = c("start.date", "end.date"))

ggplot(mdfr, aes(value, name, color=color.ful)) + 
  geom_line(size = 6) +
  xlab(NULL) + 
  ylab(NULL)+
  geom_dl(aes(label = mdfr$name), method = list(dl.combine("last.points"), cex = 0.8))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  scale_x_date(date_labels = "%b %Y", limits=as.Date(c("2018-09-01", "2021-01-01")))+
  ggtitle("Clinical Trial Timeline")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("timeline.pdf", width = 6, height = 4)



