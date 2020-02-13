library(ggplot2)
source("DEFunctions.R")
library(gridExtra)
source("Visu.R")
data <- read.table("AR1174qPCRControlesCinetiqueFaibleNitrate.txt", h = T, sep = '\t')
#data <- data[!grepl(">", data$Status),]
df <- data[c("Pos", "Cp")]

genes <- list("A" = "UBQ10", "B" = "ACT2", "C"= "NRT2.1", "D"= "NRT1.1", "E" = "BT1", "F"="BT2", "G"="HRS1")
AGIToNames <- list("AT4G05320" = "UBQ10", "AT3G18780" = "ACT2", "AT1G08090"= "NRT2.1", "AT1G12110"= "NRT1.1", "AT5G63160" = "BT1", "AT3G48360"="BT2", "AT1G13300"="HRS1")
TargetGenes <- c("AT3G18780","AT5G63160","AT3G48360", "AT1G13300", "AT1G12110", "AT1G08090", "AT4G05320")
names <- as.vector(unlist(AGIToNames[TargetGenes]))


getGene <- function(letter){return(genes[[letter]])}
df$gene <- as.vector(sapply(substr(df$Pos, 1, 1), getGene))

getCondition <- function(Pos){
  number <- as.numeric(substr(Pos,2, nchar(as.vector(Pos))))
  if(number <9) return("aCO2")
  else return("eCO2")
}
df$co2 <- sapply(df$Pos, getCondition)

ggplot(data = df, aes(x=interaction(co2, gene), y = Cp, col = co2)) + geom_point(alpha = 0.7,size = 3 ) + ggtitle("qPCR Cp values") +
  theme(axis.text.x = element_text(angle = 320,
                                   hjust = 0, colour = "grey50", size = 14), plot.title = element_text( size = 14, face = "bold"))

act2 <- getExpression("AT3G18780", c("cnF", "CnF"))
bt1 <- getExpression("AT5G63160", c("cnF", "CnF"))
bt2 <- getExpression("AT3G48360", c("cnF", "CnF"))
hrs1 <- getExpression("AT1G13300", c("cnF", "CnF"))
nrt1.1 <- getExpression("AT1G12110", c("cnF", "CnF"))
nrt2.1 <- getExpression("AT1G08090", c("cnF", "CnF"))
ubq10 <- getExpression("AT4G05320", c("cnF", "CnF"))



grid.arrange(bt1, bt2, hrs1, nrt1.1, nrt2.1, nrow = 1)

heatmapPerso(normalized.count, genes=c("AT5G63160","AT3G48360", "AT1G13300", "AT1G12110", "AT1G08090"))
heatmapPerso(normalized.count, genes=TargetGenes, conds = c("cnF", "CnF"))

data <- read.table("AR1174qPCRControlesCinetiqueFaibleNitrate.txt", h = T, sep = '\t')
data <- data[data$Cp < 26,c("Cp", "Pos")]
data <- na.omit(data)
diffCp <- data$Cp[seq(from = 1, to = dim(data)[1]-1, by = 2)] - data$Cp[seq(from = 2, to = dim(data)[1], by = 2)]
labels <- as.vector(data$Pos[seq(from = 1, to = dim(data)[1]-1, by = 2)])
meanCp <- (data$Cp[seq(from = 1, to = dim(data)[1]-1, by = 2)] + data$Cp[seq(from = 2, to = dim(data)[1], by = 2)])/2

combinedData <- data.frame(Pos = labels, diffCp = diffCp, meanCp = meanCp)

combinedData$gene <- as.vector(sapply(substr(combinedData$Pos, 1, 1), getGene))
combinedData$co2 <- sapply(combinedData$Pos, getCondition)
combinedData

                       
                       