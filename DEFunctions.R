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
suppressMessages(library(coseq, warn.conflicts = F, quietly = T))
suppressMessages(library(ade4, warn.conflicts = F, quietly = T))
suppressMessages(library(adegraphics, warn.conflicts = F, quietly = T))
suppressMessages(library(RColorBrewer, warn.conflicts = F, quietly = T))



#shiny::runGitHub("TCC-GUI", "swsoyee", subdir = "TCC-GUI", launch.browser = TRUE)
#GOMicroTom <- read.table("~/Documents/MicroTom/Shared_Genes_Solyc_Sly1.1.csv", sep = '\t', h = T)
mart = useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "athaliana_eg_gene")

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

getExpression <- function(gene, conds = "all", specie = "At"){
  # Plots the expression levels of a given gene, using the normized.count data provoded.
  # conditions are all the columns of the data by default, or can be specified
  # biological replicated should be identified by _n
  if(specie == "At") load("./normalized.count_At.RData")
  if(specie == "Sl") load("./normalized.count_Sl.RData")
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

OntologyProfile <- function(ids, specie="At"){
  #Plot ontology enrichment stats of a given a set of entrezgene IDs
  # only for Arabidopsis
  if(specie == "At"){
    results <- getBM( filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "description", "external_gene_name", "entrezgene_id"),
                      values = ids, mart = mart)
    results <- results[!rownames(results) %in% which(duplicated(results$ensembl_gene_id)), ]
    kable(results)
    
    ego <- enrichGO(gene = results$entrezgene_id,
                    OrgDb = org.At.tair.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable = TRUE)
    
    #simpOnt <- simplify(ego)
    # , cutoff=0.7, by="p.adjust", select_fun=min
    #simpOnt@result$Description
    print(barplot(ego, showCategory = 40, font.size = 5))
    #print(dotplot(ego, showCategory = 40, font.size = 5))
    print(emapplot(ego, layout = "kk"))
    return(results)
  }
  if(specie=="Sl"){
    ont <- GOMicroTom[GOMicroTom$Sly %in% genes1$gene_id,]
    ont <- ont[match(genes1$gene_id[which(genes1$gene_id %in% ont$Sly)], ont$Sly),]
    head(ont, )
    print(paste0(sum(!genes1$gene_id %in% GOMicroTom$Sly), "not in the Micro-Tom annotation..."))
    return(ont)
  }
}

########################################################## sample matching

dualDE <- function(data, labels, pval=0.01, method="edger", flc_filter = 0, plot=T){
  # selecting the right labels for pairwise comparison
  headers <- c(colnames(data)[(grepl(labels[1], colnames(data)))] , colnames(data)[grepl(labels[2], colnames(data))])
  data <- data[headers]
  group <- str_split_fixed(colnames(data), "_", 2)[,1]
  group <- factor(group)
  group <- relevel(group, labels[1])
  # tcc object
  tcc <- new("TCC", data, group)
  print(model.matrix(~group))
  print(colnames(data))
  #2 steps normalisation
  tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1, FDR = pval, floorPDEG = 0.05)
  print(tcc$norm.factors)
  tcc$DEGES$execution.time
  s <- sample(rownames(tcc$count), size = 200)
  normalized.count <- getNormalizedData(tcc)
  if(plot){
    heatmap(as.matrix(tcc$count[s,]), main = "Before normalisation")
    heatmap(as.matrix(normalized.count[s,]), main = "After normalisation")
  }
  #DEtest
  tcc <- estimateDE(tcc, test.method = method, FDR = pval, design = model.matrix(~group))
  result <- getResult(tcc, sort = TRUE)
  DEgenes <- subset(result,estimatedDEG==1 & abs(m.value) > flc_filter)
  DEgenes$upreg = ifelse(DEgenes$m.value > 0, 1, 0)
  print(paste(dim(DEgenes)[1], " genes DE"))
  head(result)
  if(plot){
    plotMDS(normalized.count, main="Multidimensional scaling plot of distances between gene expression profiles")
    plot(tcc)
    heatmap(normalized.count[DEgenes$gene_id,])
  }
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

###################### Summary  graph of many dualDE from saved objects

plotDEGNumber <- function(){
  load("genesco2Sl.RData")
  load("genesIronSl.RData")
  load("genesNitrateSl.RData")
  
  load("genesco2At.RData")
  load("genesIronAt.RData")
  load("genesNitrateAt.RData")
  
  genesco2Sl$factor = "CO2"
  genesNitrateSl$factor = "Nitrate"
  genesIronSl$factor = "Iron"
  
  
  dfSl <- rbind(genesco2Sl, genesIronSl, genesNitrateSl)
  dfSl$Specie <- "Micro-Tom"
  
  genesco2At$factor = "CO2"
  genesNitrateAt$factor = "Nitrate"
  genesIronAt$factor = "Iron"
  
  
  dfAt <- rbind(genesco2At, genesIronAt, genesNitrateAt)
  dfAt$Specie <- "A. thaliana"
  ggplot(dfAt, aes(fill=reg, y=value, x=variable)) +
    geom_bar(position="stack", stat="identity", size=1, alpha=0.5) + 
    ggtitle("CO2, Nitrate and Iron effet on gene regulation - A") + scale_fill_discrete("Set2") +
    xlab("") + ylab("Number of differentially expressed genes") + coord_flip() + facet_grid(~ factor) 
  
  df <- rbind(dfSl, dfAt)
  
  ggplot(df, aes(fill=reg, y=value, x=variable)) +
    geom_bar(position="stack", stat="identity", size=0.3, alpha=0.5, color = "black") + 
    ggtitle("CO2, Nitrate and Iron effet on gene regulation") +
    xlab("") + ylab("Number of differentially expressed genes") + coord_flip() + facet_grid(~ Specie) 
}

######################## Poisson mixture model for gene clustering on expression

clustering <- function(DEgenes, data, conds="all", nb_clusters = 2:12){
  if (length(conds) ==1){
    conds = colnames(data)
    groups <- str_split_fixed(conds, '_', 2)[,1]
  }
  dataC <- data[DEgenes,conds]
  conds <- str_split_fixed(colnames(dataC), '_', 2)[,1]
  run_pois <- coseq(dataC, conds=groups, K=nb_clusters, model="Poisson",iter = 5, transformation = "none")
  print(coseq::plot(run_pois, conds = groups, collapse_reps="average", graphs = c("ICL", "boxplots", "profiles", "probapost_barplots")))
  summary(run_pois)
  clusters_per_genes <- coseq::clusters(run_pois)
  return(clusters_per_genes)
}


###################### ACP on gene expression, and visualisation of clusters

ACP <- function(data, clusteredGenes){
  acp <- dudi.pca(log(data+0.1), center = TRUE, scale = TRUE, scannf = FALSE, nf = 4)
  print(summary(acp))
  data$cluster = clusteredGenes[as.vector(rownames(data))]
  
  s.corcircle(acp$co, xax=1, yax=2, fullcircle = FALSE, pback.col = "lightgrey")
  adegraphics::s.class(acp$li, xax =  1, yax = 2, as.factor(data$cluster), labels = as.character(levels(as.factor(data$cluster))), col = brewer.pal(n = 10, name = "Paired"), chullSize = 1, ellipseSize = 0, plabels.cex = 0.7, pbackground.col = "grey85", main = "Clusters dans le plan principal", ylim = c(-9, 9))
  
  s.corcircle(acp$co, xax=2, yax=3, fullcircle = FALSE, pback.col = "lightgrey")
  
  adegraphics::s.class(acp$li, xax =  2, yax = 3, as.factor(data$cluster), labels = as.character(levels(as.factor(data$cluster))), col = brewer.pal(n = 10, name = "Paired"), chullSize = 1, ellipseSize = 0, plabels.cex = 0.7, pbackground.col = "grey85", main = "Clusters dans le plan principal")
  
  s.corcircle(acp$co, xax=4, yax=2, fullcircle = FALSE, pback.col = "lightgrey")
  
  adegraphics::s.class(acp$li, xax =  4, yax = 2, as.factor(data$cluster), labels = as.character(levels(as.factor(data$cluster))), col = brewer.pal(n = 10, name = "Paired"), chullSize = 1, ellipseSize = 0, plabels.cex = 0.7, pbackground.col = "grey85", main = "Clusters dans le plan principal")
}
