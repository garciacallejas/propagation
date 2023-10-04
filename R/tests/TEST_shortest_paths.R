
# test of different ways to obtain shortest path matrices, as I find it a bit
# confusing still

library(maotai)
library(intsegration)
library(e1071)
library(Rfast)

A <- matrix(c(0,1.5,0,0,1.5,0,4,1,0,4,0,2,0,1,2,0),nrow = 4)
inv.A <- 1/A

# maotai
phi.1 <- maotai::shortestpath(A)
phi.inv.1 <- maotai::shortestpath(inv.A)

# bertagnolli

#this works
inv.A[is.infinite(inv.A)] <- .Machine$double.xmax
all.shortest.pahts <- intsegration::rcpp_floyd_flow(inv.A)
#flows
phi.2 <- all.shortest.pahts$F

# e1071
A.2 <- A
A.2[A == 0] <- Inf
inv.A.2 <- inv.A
inv.A.2[inv.A.2 == 0] <- NA
phi.3 <- e1071::allShortestPaths(A.2)

# Rfast
phi.4 <- Rfast::floyd(A.2)
phi.inv.4 <- Rfast::floyd(inv.A)

# see notebook for a written explanation of why intsegration is the only
# proper way

# basically it keeps track of shortest-path distances and their 
# associated flows, through D and F matrices below
# this is different from floyd on any matrix.

# int n = C.nrow();
# NumericMatrix D(n, n);
# NumericMatrix F(n, n);
# 
# for (int i = 0; i < n; i++) {
#   for (int j = 0; j < n; j++) {
#     D[i + n * j] = C[i + n * j];
#     F[i + n * j] = 1 / C[i + n * j];
#   }
# }
# for (int i = 0; i < n; i++) {
#   D[i + n * i] = 0;              /* no self cycle */
#     F[i + n * i] = 0;              /* no self cycle */
# }
# for (int k = 0; k <n; k++) {
#   for (int i = 0; i < n; i++) {
#     for (int j = 0; j < n; j++) {
#       if (D[i + n * k] + D[k + n * j] < D[i + n * j]) {
#         D[i + n * j] = D[i + n * k] + D[k + n * j];
#         F[i + n * j] = F[i + n * k] + F[k + n * j];
#       }
#     }
#   }
# }


