generateBA = function(n = 100, n0 = 2){
  mat = matrix(0, nrow= n, ncol = n)
  for(i in 1:n0){
    for(j in 1:n0){
      if(i != j){
        mat[i,j] = 1
        mat[j,i] = 1
      }
    }
  }
  for(i in n0:n){
    list = c()
    for(k in 1:(i-1)){
      list = c(list, sum(mat[,k]))
    }
    link = sample(c(1:(i-1)), size = 1, prob = list)
    mat[link,i] = 1
    mat[i,link] = 1
  }
  return(mat)
} 