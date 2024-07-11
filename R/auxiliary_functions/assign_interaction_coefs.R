
#' Assign interaction coefficients to a sign matrix
#'
#' @param sign.matrix sign matrix (coefficientes 0,1,-1)
#' @param mean.a mean of the non-diagonal coefficients. Positive mean for a>0, negative for a<0. 
#' Signs are already specified, so the normal distribution from where coefficients are drawn is truncated at 0.
#' @param sd.a standard deviation of the non-diagonal coefficients
#' @param mean.d mean of the diagonal coefficients. No sign is specified here, so in theory it's possible to 
#' specify positive diagonal coefficients. Default is -1
#' @param sd.d standard deviation of the diagonal coefficients. Default is 0.
#' @param symmetric are predator-prey coefficients symmetric?
#'
#' @return coefficients matrix
#' @export
#'
#' @examples
# sign.matrix <- matrix(c(-1,0,1,0,-1,0,-1,0,-1),nrow = 3)
# mean.a <- 0
# sd.a <- 0.25

assign_interaction_coefs <- function(sign.matrix,mean.a,sd.a,mean.d = -1,sd.d = 0,symmetric = TRUE){
  S <- nrow(sign.matrix)
  # non-diagonal coefficients
  coef.matrix <- abs(matrix(truncnorm::rtruncnorm(S^2,a = 0,mean = mean.a,sd = sd.a),nr=S,nc=S))
  coef.matrix <- coef.matrix * sign.matrix
  # diagonal coefficients
  diag(coef.matrix) <- rnorm(S,mean.d,sd.d)
  
  # if antagonistic coefficients are symmetric, update it
  # note that other interactions will not be symmetric
  # check only the upper diagonal
  if(symmetric){
    for(i.row in 1:(S-1)){
      for(i.col in (1+i.row):S){
        if((coef.matrix[i.row,i.col] < 0 & coef.matrix[i.col,i.row] > 0) | 
           (coef.matrix[i.row,i.col] > 0 & coef.matrix[i.col,i.row] < 0)){
          coef.matrix[i.col,i.row] <- -coef.matrix[i.row,i.col]
        }# if antagonistic interaction
      }# for i.col
    }# for i.row
  }# if symmetric antagonistic coefs
  
  return(coef.matrix)
}
