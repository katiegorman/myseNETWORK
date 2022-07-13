save.image("emergence_counts.RDS")

## code to plot emergence counts pretty


library(ggplot2)
library(RColorBrewer)
library(egg) ## to plot two side by side
library(patchwork) ## to plot 2 and add common legend

load("emergence_counts.RDS")

## bring in csv
emer <- read.csv("emergence counts.csv", header = TRUE)
head(emer)

## assign obs numbers so things are in the right order but we don't have big date gaps
emer$obs <- 1:nrow(emer)

## separate out by year
emer18 <- emer[emer$Year == "2018",]
emer19 <- emer[emer$Year == "2019",]


cols <- c("Bat house" = "#66c2a5", "Telephone pole" = "#8da0cb", "Tree" = "#fc8d62
")

## 2018 plot
p18 <- ggplot(emer18, aes(obs, Bats, fill = Type)) + geom_bar(stat = "identity") +
  theme_minimal() + scale_fill_manual(values=c("#66c2a5", "#8da0cb", "#000000")) +
  geom_text(aes(label=Bats), vjust=1.6, color="white", size=3.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("2018") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of bats emerged")# +
  #theme(legend.position="none")
p18



## 2019 plot
p19 <- ggplot(emer19, aes(obs, Bats, fill = Type)) + geom_bar(stat = "identity") +
  theme_minimal() + scale_fill_manual(values=c("#8da0cb","#000000")) +
  geom_text(aes(label=Bats), vjust=1.6, color="white", size=3.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("2019") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
 theme(legend.position="none") +
  ylim(0,10)
  
p19



## throw them both into one thing with combined legend
combined <- p18 + (p19 + theme(legend.position = "none")) + 
  plot_layout(guides = "collect") #& theme(legend.position = "bottom")
combined


## export plots to tiff
tiff(file = c("G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/Figures/emergence.tiff"), 
     units = "in", width = 6, height = 3, res = 300, compression = 'lzw')

combined

dev.off()
