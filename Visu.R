suppressMessages(library(stringr, warn.conflicts = F, quietly = T))
suppressMessages(library(tidyr, warn.conflicts = F, quietly = T))
suppressMessages(library(ggplot2, warn.conflicts = F, quietly = T))


#setwd("D:/These/CombinatoireRNASeqFeNCO2/")

heatmapPerso <- function(normalized.count, genes, conds="all", specie="At", geneNames=NA){
  if(specie == "At") load("normalized.count_At.RData")
  if(specie == "Sl") load("normalized.count_Sl.RData")
  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  df <- data.frame(t(normalized.count[genes,conds]))
  
  df$condition <- str_split_fixed(rownames(df), "_", 2)[,1]
  df$exactCondition <- rownames(df)
  data <- gather(data = df,key = gene,value = expression, -condition, -exactCondition)
  print(head(data))
  exp.heatmap <- ggplot(data = data, mapping = aes(x = exactCondition, y = gene,
                                                   fill = log(expression+0.1))) +
    geom_tile() + xlab(label = "Condition") +
    facet_grid(~ condition, switch = "x", scales = "free_x", space = "free_x") + labs(fill = "Normalized expression") +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank()) + scale_fill_distiller(palette = "YlGnBu") 
  if(length(geneNames) > 1) {
    exp.heatmap = exp.heatmap + scale_y_discrete(labels=geneNames)
  }
  print(exp.heatmap)
}



getExpression <- function(gene, conds = "all", specie = "At"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n
  if(specie == "At") load("normalized.count_At.RData")
  if(specie == "Sl") load("normalized.count_Sl.RData")
  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  print(conds)
  df <- normalized.count[gene, conds]
  library(reshape2)
  d<- melt(df, silent=T)
  d$group = str_split_fixed(rownames(d), "_", 2)[,1]
  
  p <- ggplot(data = d, aes(x=group, y=value, fill=group)) + geom_dotplot(binaxis = "y", stackdir = "center") +
    theme(axis.text.x = element_text(angle = 320,
                                     hjust = 0, colour = "grey50"), plot.title = element_text( size = 14, face = "bold")) +
    ggtitle(paste("Normalized expression for ", gene))
  print(df)
  return(p)
}
