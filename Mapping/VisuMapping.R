library(ggplot2)
library(gridExtra)
setwd("~/Documents/Combinatoire/")

data <- read.csv("LogsData.csv", h=T, sep = ',')

um <- ggplot(data, aes(x=UM))+
  geom_density(color="darkblue", fill="lightblue", alpha=0.4) +
  ggtitle("Uniquely mapped reads") + coord_cartesian(xlim=c(80, 100))

mm <- ggplot(data, aes(x=MM))+
  geom_density(color="darkred", fill="red", alpha=0.4) +
  ggtitle("Multi-mapped reads") +coord_cartesian(xlim=c(0.5, 2))

reads <- ggplot(data, aes(x=MappedReads))+
  geom_density(color="darkgreen", fill="green", alpha=0.4) +
  ggtitle("Sequencing depth") 

grid.arrange(um, mm, reads, ncol = 3)


library(reshape2)
melted_res <- melt(data[c("MM", "UM")])


p<-ggplot(melted_res, aes(x=value,  fill=variable)) +
  geom_density(alpha=0.4) 
p

p<-ggplot(melted_res, aes(x=variable, y = value, fill=variable)) +
  geom_violin(alpha=0.4) + labs(title="Mapped reads for the 24 At RNASeq samples",x="", y = "% of reads")
p