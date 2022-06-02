tidy_propagation_metrics <- function(A, random.removal.list = NULL, cell.distances = NULL){
  
  communicability.matrices <- communicability(A,binary = FALSE)
  
  # tidy functions only work with dataframes, but this still works, and is fast
  comm.df <- reshape2::melt(communicability.matrices[[2]],
                            value.name = "weighted.communicability")
  
  # extract cell of sp1 and cell of sp2
  comm.df$sp1 <- sub("\\-.*", "", comm.df$Var1)
  comm.df$cell1 <- sub(".*-", "", comm.df$Var1)
  comm.df$sp2 <- sub("\\-.*", "", comm.df$Var2)
  comm.df$cell2 <- sub(".*-", "", comm.df$Var2)
  
  # comm.df$scaled.binary.communicability <- scales::rescale(comm.df$binary.communicability, 
  #                                                          to = c(0,1))
  # 
  # # this should be valid because the two matrices have the same dimensions and names
  # dfw <- reshape2::melt(communicability.matrices[[2]],
  #                       value.name = "weighted.communicability")
  # comm.df$weighted.communicability <- dfw$weighted.communicability
  
  # -------------------------------------------------------------------------
  # get path lengths as well
  graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                                 mode = "undirected",
                                                 weighted = "1",diag = FALSE)
  # this takes ~10min
  path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
  
  # turning it to df is quick
  df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
  
  df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
  df.path.lengths$cell1 <- sub(".*-", "", df.path.lengths$Var1)
  df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
  df.path.lengths$cell2 <- sub(".*-", "", df.path.lengths$Var2)
  
  comm.df$Var1 <- NULL
  comm.df$Var2 <- NULL
  
  df.path.lengths$Var1 <- NULL
  df.path.lengths$Var2 <- NULL
  
  # sp1,cell1,sp2,cell2 should be common
  comm.df <- left_join(comm.df,df.path.lengths)
  
  # -------------------------------------------------------------------------
  # get spatial distance between cells
  # comm.df <- left_join(comm.df,cell.distances,by = c("cell1" = "cell_from",
  #                                                    "cell2" = "cell_to"))
  
  comm.df <- comm.df[,c("sp1","cell1",
                        "sp2","cell2",
                        # "binary.communicability",
                        # "scaled.binary.communicability",
                        "weighted.communicability",
                        "shortest.path.length")]
                        # "distance")]
  
  comm.df.aggregated <- comm.df %>%
    mutate(shortest.path.length = na_if(shortest.path.length, Inf)) %>%
    summarise(#binary_communicability_mean = mean(binary.communicability,na.rm = T),
      #binary_communicability_sd = sd(binary.communicability,na.rm = T),
      weighted_communicability_mean = mean(weighted.communicability,na.rm = T),
      weighted_communicability_sd = sd(weighted.communicability,na.rm = T),
      shortest_path_length_mean = mean(shortest.path.length,na.rm = T),
      shortest_path_length_sd = sd(shortest.path.length,na.rm = T))

    
  # -------------------------------------------------------------------------
  # global communication efficiency - a network level property
  
  landscape.graph <- igraph::graph_from_adjacency_matrix(A,
                                                         weighted = T)
  landscape.gce <- GCE_weighted(g = landscape.graph,normalised = T,
                                directed = T)
  
  gce.df <- data.frame(normalised_gce = landscape.gce$normalised)
  
  metrics.df <- cbind(comm.df.aggregated,gce.df) %>%
    pivot_longer(everything(),names_to = "metric",values_to = "value")
  
# -------------------------------------------------------------------------
# if a list of matrices with random removals are supplied, I can also obtain the relative 
  # contribution of the removed sp as compared to these random removals
  
  if(!is.null(random.removal.list)){
    
    rand.metrics <- list()
    
    for(i in 1:length(random.removal.list)){
      Ar <- random.removal.list[[i]]
      
      num.interspecific.links <- sum(Ar != 0 & Ar != 1)
      
      if(num.interspecific.links != 0){
      
      communicability.matrices <- communicability(Ar,binary = F)
      
      # tidy functions only work with dataframes, but this still works, and is fast
      comm.rand.df <- reshape2::melt(communicability.matrices[[2]],
                                value.name = "weighted.communicability")
      
      # extract cell of sp1 and cell of sp2
      comm.rand.df$sp1 <- sub("\\-.*", "", comm.rand.df$Var1)
      comm.rand.df$cell1 <- sub(".*-", "", comm.rand.df$Var1)
      comm.rand.df$sp2 <- sub("\\-.*", "", comm.rand.df$Var2)
      comm.rand.df$cell2 <- sub(".*-", "", comm.rand.df$Var2)
      
      # comm.rand.df$scaled.binary.communicability <- scales::rescale(comm.rand.df$binary.communicability, 
      #                                                          to = c(0,1))
      # 
      # # this should be valid because the two matrices have the same dimensions and names
      # dfw <- reshape2::melt(communicability.matrices[[2]],
      #                       value.name = "weighted.communicability")
      # comm.rand.df$weighted.communicability <- dfw$weighted.communicability
      
      # -------------------------------------------------------------------------
      # get path lengths as well
      graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = Ar,
                                                     mode = "undirected",
                                                     weighted = "1",diag = FALSE)
      # this takes ~10min
      path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
      
      # turning it to df is quick
      df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
      
      df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
      df.path.lengths$cell1 <- sub(".*-", "", df.path.lengths$Var1)
      df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
      df.path.lengths$cell2 <- sub(".*-", "", df.path.lengths$Var2)
      
      comm.rand.df$Var1 <- NULL
      comm.rand.df$Var2 <- NULL
      
      df.path.lengths$Var1 <- NULL
      df.path.lengths$Var2 <- NULL
      
      # sp1,cell1,sp2,cell2 should be common
      comm.rand.df <- left_join(comm.rand.df,df.path.lengths)
      
      # -------------------------------------------------------------------------
      # get spatial distance between cells
      # comm.rand.df <- left_join(comm.rand.df,cell.distances,by = c("cell1" = "cell_from",
      #                                                    "cell2" = "cell_to"))
      
      comm.rand.df <- comm.rand.df[,c("sp1","cell1",
                            "sp2","cell2",
                            # "binary.communicability",
                            # "scaled.binary.communicability",
                            "weighted.communicability",
                            "shortest.path.length")]
      # "distance")]
      
      comm.rand.df.aggregated <- comm.rand.df %>%
        mutate(shortest.path.length = na_if(shortest.path.length, Inf)) %>%
        summarise(#binary_communicability_mean = mean(binary.communicability,na.rm = T),
          #binary_communicability_sd = sd(binary.communicability,na.rm = T),
          weighted_communicability_mean = mean(weighted.communicability,na.rm = T),
          weighted_communicability_sd = sd(weighted.communicability,na.rm = T),
          shortest_path_length_mean = mean(shortest.path.length,na.rm = T),
          shortest_path_length_sd = sd(shortest.path.length,na.rm = T))
      
      # -------------------------------------------------------------------------
      # global communication efficiency - a network level property
      
      landscape.graph <- igraph::graph_from_adjacency_matrix(Ar,
                                                             weighted = T)
      landscape.gce <- GCE_weighted(g = landscape.graph,normalised = T,
                                    directed = T)
      gce.rand.df <- data.frame(normalised_gce = landscape.gce$normalised)
      
      rand.metrics[[length(rand.metrics)+1]] <- cbind(comm.rand.df.aggregated,gce.rand.df)
      
      }# if interspecific links (i.e. valid matrix)
    }# for each random removal
    
    rand.df <- bind_rows(rand.metrics)
    
    rand.long <- rand.df %>% pivot_longer(everything(),
                              names_to = "metric",
                              values_to = "value") %>%
      group_by(metric) %>%
      summarise(rand.avg.value = mean(value))
    
    metrics.df <- left_join(metrics.df,rand.long)
    metrics.df$sp.contribution <- metrics.df$value - metrics.df$rand.avg.value
    metrics.df$rand.avg.value <- NULL
  }# if random removals
  
# -------------------------------------------------------------------------
# return the long dataframe including all aggregated metrics
  
  return(metrics.df)
  
}