
#' communicability of an adjacency matrix
#'
#' Calculates communicability (10.1103/PhysRevE.77.036111) for a matrix
#' of arbitrary size. It returns binary and, if the matrix is weighted,
#' weighted communicability. The latter is calculated simply as the binary 
#' counterpart, by exponentiating the scaled weighted matrix.
#'
#' @param A numeric matrix
#'
#' @return list with two components: binary.comm.matrix and weighted.comm.matrix
#' @export
#'
#' @examples
communicability <- function(A){
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  # is my matrix binary?
  is.binary <- sum(A > 1) == 0
  if(is.binary){
    # binary
    # beware, it takes a while for a 200M elements matrix - RAM is a limiting factor
    binary.comm.matrix <- expm::expm(A)
    weighted.comm.matrix <- 0
  }else{
    # binary and weighted
    binary.matrix <- ifelse(A != 0, 1, 0)
    binary.comm.matrix <- expm::expm(binary.matrix)
    # TODO consider weighting these elements by community size
    
    # weighted
    
    scaled.comm.matrix <- range01(A)
    weighted.comm.matrix <- expm::expm(scaled.comm.matrix)
    
    # below, a previous approach, from Crofts and Higham 2009
    # not used because simply exponentiating the weighted matrix gives 
    # good results, coherent with binary communicability values.
    
    # abs.int.matrix <- abs(A)
    # degree.diag.matrix <- matrix(0,nrow = nrow(binary.matrix),
    #                              ncol = ncol(binary.matrix))
    # diag(degree.diag.matrix) <- rowSums(abs.int.matrix) # strength - suma de los pesos
    # 
    # d.inverse <- solve(degree.diag.matrix)
    # d.pow <- sqrt(d.inverse)
    # 
    # multiplied.matrix <- d.pow %*% abs.int.matrix %*% d.pow # eq. 2.2 of Crofts and Higham 2009
    # weighted.comm.matrix <- expm::expm(multiplied.matrix)
  }
  
  return(list(binary.matrix = binary.comm.matrix, 
              weighted.matrix = weighted.comm.matrix))
  
}