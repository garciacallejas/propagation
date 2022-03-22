
library(extraDistr)
library(igraph)
library(gamlss.dist)

# -------------------------------------------------------------------------
richness <- 50
degree.gradient <- seq(1,20,length.out = 10)
network.replicates <- 10

# interaction strength parameters
int.mean <- 0
int.sd <- 1
tau <- 1.5

for(i.dist in 1:length(degree.gradient)){
  
  for(i.net in 1:network.replicates){
    
    my.dist <- extraDistr::rtpois(n = richness,lambda = degree.gradient[i.dist],a = 0)
    my.net <- igraph::sample_degseq(my.dist,method = "vl")
    my.matrix <- as.matrix(igraph::as_adjacency_matrix(my.net,type = "both"))
    # cat("lambda:",degree.gradient[i.dist],"- connectance:",(sum(my.dist)/richness^2),"\n")
    
    # assign interaction strengths according to an "extended" normal dist
    weights <- abs(gamlss.dist::rSHASHo(sum(my.dist), mu = int.mean, 
                             sigma = int.sd, nu = 0, tau = tau))
    my.matrix[my.matrix == 1] <- weights
    diag(my.matrix) <- 1
    
  }# for i.net
}# for i.dist



