
#' Communicability at the network level
#'
#' needs the intsegration package (https://github.com/gbertagnolli/intsegration)
#'
#' @param A interaction matrix, with values > 0
#' @param weighted whether to calculate weighted communicability or, otherwise, binary
#' @param structural.zeros optional matrix of the same dimensions as A, 
#' with zeroes in links that are structurally forbidden
#' @param return.pairwise.comm whether to return pairwise communicability matrices
#'
#' @return list with two components: communicability at the network level (raw.comm) 
#' and communicability of the maximally connected network (max.comm). Optionally,
#' the pairwise communicability matrices can be returned as two more components.
#' 
#' @export
#'
#' @examples
communicability_network <- function(A,
                                    weighted = TRUE,
                                    structural.zeros = NULL, 
                                    return.pairwise.comm = FALSE){
  
  N <- nrow(A)
  
  # ensure diagonals are set to zero
  diag(A) <- 0
  
  com.list <- communicability(A,weighted = weighted,return.scaled = F)
  bin.comm <- com.list[[1]]
  w.comm <- com.list[[2]]
  
  # aggregated metric
  if(weighted){
    comm_obs <- 1. / N / (N - 1) * sum(w.comm, na.rm = T)
  }else{
    comm_obs <- 1. / N / (N - 1) * sum(bin.comm, na.rm = T)
  }
  
  inv.A <- 1/A
  #this works
  inv.A[is.infinite(inv.A)] <- .Machine$double.xmax
  all.shortest.paths <- intsegration::rcpp_floyd_flow(inv.A)
  #flows
  phi.A <- all.shortest.paths$F
  
  # this is wrong
  # phi.A <- maotai::shortestpath(A)
  # for numerical convergence
  # phi.A[is.infinite(phi.A)] <- .Machine$double.xmax
  
  # structural zeros
  if(!is.null(structural.zeros)){
    phi.A[structural.zeros == 0] <- 0
  }
  
  # diagonal must be zero as well, check just in case
  diag(phi.A) <- 0
  
  # and for numerical rounding, assign zeros to extremely low values
  # since these are spurious
  phi.A[phi.A < 1e-10] <- 0
  
  # communicability of the "ideal" matrix
  comm.phi <- communicability(phi.A,return.scaled = F)
  bin.comm.phi <- comm.phi[[1]]
  w.comm.phi <- comm.phi[[2]]
  
  # aggregated metric
  if(weighted){
    comm_ideal <- 1. / N / (N - 1) * sum(w.comm.phi, na.rm = T)
  }else{
    comm_ideal <- 1. / N / (N - 1) * sum(bin.comm.phi, na.rm = T)
  }
  
  if(return.pairwise.comm){
    return(list(raw.comm = comm_obs,
                max.comm = comm_ideal,
                bin.pairwise.comm = bin.comm,
                weighted.pairwise.comm = w.comm))
  }else{
    # return the ratio, or the raw value
    return(list(raw.comm = comm_obs,
                max.comm = comm_ideal))
  }
  
  
}
