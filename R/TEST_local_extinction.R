
num.nodes <- 50
deg.lambda <- 10

my.dist <- extraDistr::rtpois(n = num.nodes,lambda = deg.lambda,a = 0)

# make sum even
if (sum(my.dist) %% 2 != 0) { my.dist[1] <- my.dist[1] + 1 }

my.net <- igraph::sample_degseq(my.dist,method = "vl")
my.metaweb <- as.matrix(igraph::as_adjacency_matrix(my.net,type = "both"))
weights <- abs(gamlss.dist::rSHASHo(sum(my.dist), mu = 0, 
                                    sigma = .1, nu = 0, tau = 3))
my.metaweb[my.metaweb == 1] <- weights
diag(my.metaweb) <- 1

rem <- nrow(my.metaweb)

met.list <- list()
rem.matrix <- my.metaweb
all.diag <- F
while(rem >2 & !all.diag){
  
  sp.pos <- sample(nrow(rem.matrix),1)
  
  rem.matrix <- rem.matrix[-sp.pos,-sp.pos]
  if(sum(rem.matrix) == nrow(rem.matrix)){
    all.diag <- TRUE
  }else{
  
  rem <- rem - 1
  com <- communicability(A = rem.matrix)
  
# -------------------------------------------------------------------------
  comm.df <- reshape2::melt(com[[2]],
                            value.name = "com")
  comm2 <- comm.df %>% summarise(com.mean = mean(com,na.rm = T),
    com.sd = sd(com,na.rm = T))
# -------------------------------------------------------------------------
  rem.graph <- igraph::graph_from_adjacency_matrix(rem.matrix,
                                                         weighted = T)
  gce <- GCE_weighted(g = rem.graph,normalised = T,
                                directed = T)
# -------------------------------------------------------------------------
  met.list[[length(met.list)+1]] <- data.frame(removed = nrow(my.metaweb) - rem,
                                               comm.mean = comm2$com.mean,
                                               comm.sd = comm2$com.sd,
                                               gce = gce$normalised)
  }# if only diagonal elements
}# while

met.df <- bind_rows(met.list)

# -------------------------------------------------------------------------
met.long <- met.df %>% pivot_longer(cols = c(comm.mean,comm.sd,gce),names_to = "metric")

ggplot(met.long, aes(x = removed, y = value)) + 
  geom_point() + 
  facet_grid(metric~., scales = "free_y") +
  NULL






