library(ggplot2)
library("tidyr")
library(stringr)

setwd("D:/These/CombinatoireRNASeqFeNCO2/")
load("normalized.count_At.RData")
load("GenesNitrate_At.RData")

genes <- sharedBy3[1:20]

df <- data.frame(t(normalized.count[genes,]))
df$condition <- str_split_fixed(rownames(df), "_", 2)[,1]
df$exactCondition <- rownames(df)

data <- gather(data = df,
                   key = gene, 
                   value = expression, -condition, -exactCondition)


exp.heatmap <- ggplot(data = data, mapping = aes(x = exactCondition,
                                                     y = gene,
                                                     fill = log(expression+0.1))) +
  geom_tile() +
  xlab(label = "Condition") +
  # facet_grid makes two panels, one for control, one for flu:
  facet_grid(~ condition, switch = "x", scales = "free_x", space = "free_x") + labs(fill = "Normalized expression") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank()) + scale_fill_distiller(palette = "YlGnBu")

exp.heatmap
