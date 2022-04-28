
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
library(expm) # for matrix exponentiation
library(tidyverse)

# this is a test
# a.rows <- 50
# a.cols <- 50
# c <- 0.2
# l <- c * (a.rows*a.cols)
# 
# A <- matrix(0,nrow = a.rows,ncol = a.cols)
#   ints <- abs(gamlss.dist::rSHASHo(l, mu = 0, 
#                                    sigma = 1, nu = 0, tau = 1))
#   # hist(ints,breaks = 20)
#   for(i in 1:l){
#     my.sample.row <- sample(1:a.rows,1,replace = T)
#     my.sample.col <- sample(1:a.cols,1,replace = T)
#     
#     while(A[my.sample.row,my.sample.col] != 0 & 
#           my.sample.row == my.sample.col){
#       my.sample.row <- sample(1:a.rows,1,replace = T)
#       my.sample.col <- sample(1:a.cols,1,replace = T)
#     }
#     A[my.sample.row,my.sample.col] <- ints[i]
#   }
# 
# 
# A.scaled <- range01(A)
# diag(A.scaled) <- 1
# A.bin <- A
# A.bin[A.bin>0] <- 1
# 
# commw <- expm(A.scaled)
# commb <- expm(A.bin)
# 
# df <- expand.grid(i = 1:a.rows,j = 1:a.cols,weight = 0,binary.comm = 0,weighted.comm = 0)
# 
# for(i in 1:nrow(df)){
#   df$binary.comm[i] <- commb[df$i[i],df$j[i]]
#   df$weight[i] <- A.scaled[df$i[i],df$j[i]]
#   df$weighted.comm[i] <- commw[df$i[i],df$j[i]]
# }

load("results/communicability/comm_dd02_sa01_dk03_re4_30sp_400cells.RData")
comm.df$id1 <- paste(comm.df$sp1,"-",comm.df$cell1,sep="")
comm.df$id2 <- paste(comm.df$sp2,"-",comm.df$cell2,sep="")

df <- comm.df[,c("id1","id2","binary.communicability","weighted.communicability")]

df.nodiag <- subset(df,id1 != id2)
df.nodiag$scaled.bin.comm <- range01(df.nodiag$binary.communicability)

df.nodes <- df.nodiag %>% 
  group_by(id1) %>%
  summarise(#node.weight = sum(weight),
            node.bin.comm = sum(scaled.bin.comm),
            node.weighted.comm = sum(weighted.communicability)) %>%
  filter(node.bin.comm != 0 & node.weighted.comm != 0)

ggplot(df.nodes,aes(x = node.bin.comm,y = node.weighted.comm)) + 
  # geom_point(aes(size = node.weight), shape = 21, alpha = .5) +
  geom_point(alpha = .5) +
  theme_bw() +
  NULL

