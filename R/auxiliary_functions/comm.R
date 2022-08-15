
#' communicability at the network level
#'
#' needs the intsegration package (https://github.com/gbertagnolli/intsegration)
#'
#' @param A 
#' @param weighted 
#' @param structural.zeros 
#'
#' @return
#' @export
#'
#' @examples
comm <- function(A,
                 normalised = TRUE,
                 weighted = FALSE,
                 structural.zeros = NULL, 
                 return.pairwise.comm = FALSE){
  
    N <- nrow(A)
  
  # ensure diagonals are set to zero
  diag(A) <- 0
  
  com.list <- communicability(A)
  bin.comm <- com.list[[1]]
  w.comm <- com.list[[2]]
  
  # aggregated metric
  if(weighted){
    comm_obs <- 1. / N / (N - 1) * sum(w.comm, na.rm = T)
  }else{
    comm_obs <- 1. / N / (N - 1) * sum(bin.comm, na.rm = T)
  }
  
  
  # calculate the "ideal" network as a baseline
  if(normalised){
    
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
    
    # and for numerical rounding, assing zeros to extremely low values
    # since there are spurious
    phi.A[phi.A < 1e-10] <- 0
    
    # communicability of the "ideal" matrix
    comm.phi <- communicability(phi.A)
    bin.comm.phi <- comm.phi[[1]]
    w.comm.phi <- comm.phi[[2]]
    
    # aggregated metric
    if(weighted){
      comm_ideal <- 1. / N / (N - 1) * sum(w.comm.phi, na.rm = T)
    }else{
      comm_ideal <- 1. / N / (N - 1) * sum(bin.comm.phi, na.rm = T)
    }
  }else{
    comm_ideal <- NA
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

# this is an alternative implementation using tailored functions
# from bertagnolli (package "intsegration")
# it is a bit more difficult to follow, I think
# inv.A <- 1/A
# #this works
# inv.A[is.infinite(inv.A)] <- .Machine$double.xmax
# all.shortest.pahts <- intsegration::rcpp_floyd_flow(inv.A)
# #flows
# phi.A2 <- all.shortest.pahts$F
