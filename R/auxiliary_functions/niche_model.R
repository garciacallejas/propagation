#' Niche Model Food Web
#'
#' adapted from https://github.com/jjborrelli/trophic/blob/master/R/topological.R
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return An adjacency matrix for a niche model food web.
#' @export
#'
#' @section Reference:
#' Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. Nature 404:180–183.
#'
#' @examples
#' niche(20, .1)
niche_model <- function(S, C){
  cond <- FALSE
  while(!cond){
    n.i <- sort(runif(S), decreasing = F)
    r.i <- rbeta(S,1,((1/(2*C))-1))*n.i
    c.i <- runif(S, r.i/2, n.i)
    
    a <- matrix(0, nrow = S, ncol = S)
    
    for(i in 2:S){
      for(j in 1:S){
        if(n.i[j] > (c.i[i] - (.5 * r.i[i])) & n.i[j] < (c.i[i] + .5 * r.i[i])){
          a[j, i] <- 1
        }
      }
    }
    
    cond <- igraph::is_connected(igraph::graph_from_adjacency_matrix(a))
  }
  
  return(a)
}