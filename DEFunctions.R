# Those generic functions are made to be used for Differential Expression Analysis
# So scripts for the analysis on a particular dataset are lighter and less redundant

suppressMessages(library(TCC, warn.conflicts = F, quietly = T))
#suppressMessages(library(biomartr))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(org.At.tair.db, warn.conflicts = F, quietly = T))
suppressMessages(library(clusterProfiler, warn.conflicts = F, quietly = T))
suppressMessages(library(enrichplot, warn.conflicts = F, quietly = T))
suppressMessages(library(DESeq2, warn.conflicts = F, quietly = T))
load("./normalized.count_At.RData")

#shiny::runGitHub("TCC-GUI", "swsoyee", subdir = "TCC-GUI", launch.browser = TRUE)

########################################################## sample matching

annot <- read.csv("Code_for_RNAseq_CO2_N_Fr.csv", h = T, sep = ';')
conditions <- as.vector(unique(annot$Sample))
annot$ID <- paste0('R', annot$Code)

annot$condition <- substr(conditions, 1, nchar(conditions)-1)
cond <- unique(substr(conditions, 1, nchar(conditions)-1))
cond <- cond[grepl("At", cond)]

getCondition <- function(id){
  # get condition without replicates
  return(annot[annot$ID == id, "condition"])
}
getExactCondition <- function(id){
  # get condition with sample
  return(annot[annot$ID == id, "Sample"])
}

getLabel <- function(id, with.rep = T){
  # get condition with sample in a simplified notation
  text <- as.vector(annot[annot$ID == id, "Sample"])
  res <- ''
  nb <- substr(text, nchar(text), nchar(text))
  if(grepl("Ambient", text)){res = paste0(res, "c")}
  else{res = paste0(res, "C")}
  if(grepl("High", text)){res = paste0(res, "N")}
  else{res = paste0(res, "n")}
  if(grepl("Starv", text)){res = paste0(res, "f")}
  else{res = paste0(res, "F")}
  if(with.rep) res = paste0(res, "_", nb)
  return(res)
}

########################################################## show expression

getExpression <- function(gene, conds = "all"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n
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

########################################################## Ontology

OntologyProfileAt <- function(r){
  #Plot ontology enrichment stats of a given data frame of genes for Arabidopsis thaliana
  # the column of the gene ids must be named entrezgene_id
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

########################################################## sample matching

dualDE <- function(data, labels, pval=0.01, method="edger", flc_filter = 0){
  # selecting the right labels for pairwise comparison
  data <- data[,grepl(labels[1], colnames(data)) | grepl(labels[2], colnames(data))]
  group <- str_split_fixed(colnames(data), "_", 2)[,1]
  
  # tcc object
  tcc <- new("TCC", data, group)
  
  #2 steps normalisation
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1, FDR = pval, floorPDEG = 0.05)
  print(tcc$norm.factors)
  tcc$DEGES$execution.time
  s <- sample(rownames(tcc$count), size = 200)
  heatmap(as.matrix(tcc$count[s,]), main = "Before normalisation")
  normalized.count <- getNormalizedData(tcc)
  heatmap(as.matrix(normalized.count[s,]), main = "After normalisation")
  
  #DEtest
  tcc <- estimateDE(tcc, test.method = method, FDR = pval, design = model.matrix(~group))
  result <- getResult(tcc, sort = TRUE)
  DEgenes <- subset(result,estimatedDEG==1 & abs(m.value) > flc_filter)
  DEgenes$upreg = ifelse(DEgenes$m.value > 1, 1, 0)
  print(paste(dim(DEgenes)[1], " genes DE"))
  head(result)
  plotMDS(normalized.count, main="Multidimensional scaling plot of distances between gene expression profiles")
  plot(tcc)
  heatmap(normalized.count[DEgenes$gene_id,])
  return(DEgenes)
}

getOneFactorComparisons <- function(factor = 'C'){
  if(factor == 'C') pos = 1
  else if(factor == 'n') pos = 2
  else if(factor =='f') pos = 3
  poss = seq(1:3)
  others = poss[poss != pos]
  res = c()
  for(comp in comparables){
    labels = str_split_fixed(comp, '-', 2)
    diffFactor = substr(labels[1], pos, pos) != substr(labels[2], pos, pos)
    other1 = substr(labels[1], others[1], others[1]) == substr(labels[2], others[1], others[1])
    other2 = substr(labels[1], others[2], others[2]) == substr(labels[2], others[2], others[2])
    if (diffFactor & other1 & other2 & grepl(factor, labels[1])){
      res = c(res, comp)
    }
  }
  return(res)
}
