#             Those functions are made to be used during Differential Expression Analysis

suppressMessages(library(TCC, warn.conflicts = F, quietly = T))
suppressMessages(library(biomartr))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db, warn.conflicts = F, quietly = T))
suppressMessages(library(clusterProfiler, warn.conflicts = F, quietly = T))
suppressMessages(library(enrichplot, warn.conflicts = F, quietly = T))
suppressMessages(library(DESeq2, warn.conflicts = F, quietly = T))




# So scripts for analysis on a particular dataset are lighter and less redundant


getExpression <- function(normalized.count, gene, conds = "all"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n
  
  if (length(conds) ==1){
    conds = colnames(normalized.count)
  }
  
  else{conds = grepl(conds[1], colnames(normalized.count)) | grepl(conds[2], colnames(normalized.count))}
  df <- normalized.count[gene, conds]
  library(reshape2)
  d<- melt(df)
  d$group = str_split_fixed(rownames(d), "_", 2)[,1]
  
  p <- ggplot(data = d, aes(x=group, y=value, fill=group)) + geom_dotplot(binaxis = "y", stackdir = "center") +
    theme(axis.text.x = element_text(angle = 320,
                                     hjust = 0, colour = "grey50"), plot.title = element_text( size = 14, face = "bold")) +
    ggtitle(paste("Normalized expression for ", gene))
  p
  return(df)
}



suppressMessages(library(TCC, warn.conflicts = F, quietly = T))
suppressMessages(library(biomartr))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db, warn.conflicts = F, quietly = T))
suppressMessages(library(clusterProfiler, warn.conflicts = F, quietly = T))
suppressMessages(library(enrichplot, warn.conflicts = F, quietly = T))
suppressMessages(library(DESeq2, warn.conflicts = F, quietly = T))
OntologyProfileAt <- function(r){
  #Plot ontology enrichment stats of a given list of genes for Arabidopsis thaliana
  ego <- enrichGO(gene = r$entrezgene_id,
                  OrgDb = org.At.tair.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)
  simpOnt <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  simpOnt@result$Description
  print(barplot(simpOnt, showCategory = 40, font.size = 5))
  print(dotplot(simpOnt, showCategory = 40, font.size = 5))
  print(emapplot(simpOnt, layout = "kk"))
}