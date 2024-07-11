
#' generate a random horizontal interaction matrix
#'
#' @param S number of nodes
#' @param c connectance
#' @param tau parameter of the rSHASHO distribution
#' @param min.diag.dom minimum ratio of diagonal dominance
#' @param restricted.positive only positive values?
#' @param int.mean mean of interaction strengths
#' @param int.sd sd of interaction strenghts
#' @param scale_to_1 wheter to scale each row so that the max value is 1 and all others
#' are scaled in reference to it
#'
#' @return
#' @export
#'
#' @examples
horizontal_community_matrix <- function(S = 5,
                                        c = 0.5,
                                        tau = 1.5,
                                        min.diag.dom = 0,
                                        restricted.positive = TRUE,
                                        int.mean = 0,
                                        int.sd = 1,
                                        scale_to_1 = TRUE){
    
    a.rows <- S
    a.cols <- S
    l <- round(c * (a.rows*a.cols))
    
    A <- matrix(0,nrow = a.rows,ncol = a.cols)
    if(restricted.positive){
        ints <- abs(gamlss.dist::rSHASHo(l, mu = int.mean, 
                                         sigma = int.sd, nu = 0, tau = tau))
    }else{
        ints <- gamlss.dist::rSHASHo(l, mu = int.mean, 
                                     sigma = int.sd, nu = 0, tau = tau)
    }
    
    # randomly assign interaction strengths outside the diagonal
    for(i in 1:l){
        my.sample.row <- sample(1:a.rows,1,replace = T)
        my.sample.col <- sample(1:a.cols,1,replace = T)
        
        while(A[my.sample.row,my.sample.col] != 0 & 
              my.sample.row == my.sample.col){
            my.sample.row <- sample(1:a.rows,1,replace = T)
            my.sample.col <- sample(1:a.cols,1,replace = T)
        }
        A[my.sample.row,my.sample.col] <- ints[i]
    }# for i
    
    # diag values
    for(i.row in 1:a.rows){
        non.diag <- abs(sum(A[i.row,]))
        if(min.diag.dom > 0){
            # values around that needed to achieve dominance in this row
            A[i.row,i.row] <- abs(rnorm(1,mean = (non.diag + min.diag.dom),sd = .1))
        }else{
            # values from the same distribution as the rest
            A[i.row,i.row] <- abs(gamlss.dist::rSHASHo(1, mu = int.mean, 
                                                       sigma = int.sd, nu = 0, tau = tau))
        }
        # cat(i.row,"diag:",A[i.row,i.row],"-non diag:",non.diag,"-dominance:",A[i.row,i.row] - non.diag,"\n")
    }
    
    if(scale_to_1){
      for(i.row in 1:nrow(A)){
        my.row <- A[i.row,]
        A[i.row,] <- scales::rescale(my.row,to = c(0,1))
      }
    }
    
    return(A)
}
