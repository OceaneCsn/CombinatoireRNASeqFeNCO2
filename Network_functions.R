library(PLNmodels)
library(ggplot2)
library(psych)


PLN_network <- function(data, DEGenes, plot_path=F){
  # covariables
  groups <- str_split_fixed(colnames(data), "_", 2)[,1]
  co2 <- str_split_fixed(groups, "", 3)[,1]
  nitrate <- factor(str_split_fixed(groups, "", 3)[,2])
  nitrate <- relevel(nitrate, "N")
  fer <- factor(str_split_fixed(groups, "", 3)[,3])
  fer = relevel(fer, "F")
  covariates <- data.frame(row.names =colnames(data), co2,nitrate, fer)
  
  # preparation des données
  counts <- round(t(data[DEGenes,]), 0)
  plnData <- prepare_data(counts = counts, covariates = covariates)
  network_models <- PLNnetwork(Abundance ~ nitrate + fer + co2 +offset(log(Offset)), data = plnData)
  network_models
  network_models$criteria %>% head() %>% knitr::kable()
  plot(network_models, "diagnostic")
  plot(network_models)
  if(plot_path==T){
    coefficient_path(network_models, corr = TRUE) %>% 
      ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) + 
      geom_line(show.legend = FALSE) +  coord_trans(x="log10") + theme_bw()
  }
  
  model_StARS <- getBestModel(network_models, "StARS")
  
  net <- plot(model_StARS)
  plot.igraph(net)
  #plot.igraph(net, vertex.size = 10, vertex.label.cex = 0.5) 
  plot(model_StARS, type = "support", output = "corrplot")
  
  # Verification des predictions du modele
  data.frame(
    fitted   = as.vector(fitted(model_StARS)),
    observed = as.vector(counts)
  ) %>% 
    ggplot(aes(x = observed, y = fitted)) + 
    geom_point(size = .5, alpha =.25 ) + 
    scale_x_log10() + 
    scale_y_log10() + 
    theme_bw() + annotation_logticks()
  return(net)
}




netStats <- function(g){
  degree<- degree(g)
  betweenness<- betweenness(g, weights=NA)
  Node_nw_st<- data.frame(degree, betweenness)
  print(ggplot( data = Node_nw_st, aes(x=degree)) +geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.7) +
    ggtitle("Degree distribution") +
    theme(plot.title = element_text(size=15)))
  print(ggplot( data = Node_nw_st, aes(x=betweenness)) +geom_histogram( fill="#E69F00", color="#e9ecef", alpha=0.7) +
    ggtitle("Betweeness distribution") +
    theme(plot.title = element_text(size=15)))
  Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
  Node_nw_st <- cbind(Node_nw_st, Rank_stat)
  return(Node_nw_st[order(-Node_nw_st$Rank_stat),])
}
