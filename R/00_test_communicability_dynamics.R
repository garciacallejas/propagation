
# compare net effects (dynamics) and communicability (structure)

library(tidyverse)
library(patchwork)
source("R/niche_model.R")
source("R/GLV_functions.R")
source("R/communicability.R")
source("R/horizontal_community_matrix.R")

# -------------------------------------------------------------------------

# how many simulated networks
num.nets <- 10
# richness
S <- seq(10,100,by = 10)
# connectance
C <- .25

# threshold for considering positive abundances
# if a steady-state solution has sum abundances lower than this,
# it will be assumed that the system has no positive solution
abundance.threshold <- 1

net.type <- "horizontal"
# net.type <- "food web"

# -------------------------------------------------------------------------
# signs are the "real" ones in these functions, e.g. competition is negative
# interaction.matrix <- -interaction.matrix
# diag(interaction.matrix) <- -1
# -------------------------------------------------------------------------

net.metrics.list <- list()

for(i.s in 1:length(S)){
  my.s <- S[i.s]
  for(i.net in 1:num.nets){
    if(net.type == "horizontal"){
      A <- horizontal_community_matrix(S = my.s,
                                       c = C,
                                       min.diag.dom = 1)
    }else{
      A <- niche_model(my.s,C)
    }

    # ensure zero diagonal?
    # diag(A) <- 0
    
    # parameters of the dynamical system
    N0 <- rep(1,my.s) #initial conditions for species abundances
    # r0 <- runif(my.s,min = 0,max = 2) # assume al growth rates are > 0
    r0 <- rep(1,my.s)
    # competition must be negative
    parameters <- list(r = r0, A = A * -1)
    
    # solve the dynamical system
    eq.abund <- rootSolve::stode(y = N0,
                                 func=GLV, 
                                 parms = parameters,
                                 positive = TRUE)[[1]]
    
    # if there is a non-negative steady-state solution, it will be eq.abund
    # if there is no positive solution, all eq.abund elements will be zero
    if(sum(eq.abund) > 0){
      
      # calculate jacobian
      J = rootSolve::jacobian.full(y=eq.abund,func=GLV,parms = parameters)
      
      # negative of the inverse jacobian: net effects matrix
      invJ <- - MASS::ginv(J)
      
    }else{
      
      J <- matrix(NA,nrow = nrow(A),
                  ncol = ncol(A))
      
      invJ <- matrix(NA,nrow = nrow(A),
                     ncol = ncol(A))
    }# if-else feasible solution
    
    # calculate communicability
    comm <- communicability(A)
    bin.comm <- comm[[1]]
    w.comm <- comm[[2]]
    
    # -------------------------------------------------------------------------
    # get path lengths as well
    graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                                   mode = "undirected",
                                                   weighted = "1",diag = FALSE)
    # this takes ~10min
    path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
    
    # turning it to df is quick
    df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
    names(df.path.lengths)[1:2] <- c("sp1","sp2")
    
    # generate dataframe
    df1 <- expand.grid(net = i.net,
                       structure = net.type,
                       sp1 = 1:my.s, sp2 = 1:my.s, 
                       direct.effect = 0,
                       net.effect = 0,
                       binary.communicability = 0,
                       weighted.communicability = 0)
    
    for(i.pos in 1:nrow(df1)){
      df1$richness <- my.s
      df1$direct.effect[i.pos] <- J[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$net.effect[i.pos] <- invJ[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$binary.communicability[i.pos] <- bin.comm[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$weighted.communicability[i.pos] <- w.comm[df1$sp1[i.pos],df1$sp2[i.pos]]
      
    }
    
    df2 <- left_join(df1,df.path.lengths)
    
    net.metrics.list[[length(net.metrics.list)+1]] <- df2
    
  }# for i.net
}# for i.s
net.metrics.df <- bind_rows(net.metrics.list)

# write.csv2(net.metrics.df,"results/communicability_indirect_effects_horizontal.csv",row.names = F)

net.clean <- subset(net.metrics.df, !is.na(net.effect))
net.plot.data.subset <- subset(net.metrics.df, sp1 != sp2 & net.effect > -1 & net.effect < 1)
net.plot.data <- subset(net.metrics.df, sp1 != sp2)

# net.dist <- net.plot.data %>% 
#   group_by(richness) %>%
#   summarise(null.95.perc = qnorm(.999,mean(abs(net.effect)),sd(abs(net.effect))))

# net.plot.data.2 <- left_join(net.plot.data,net.dist)
# net.plot.data.2$obs.perc.99 <- ifelse(abs(net.plot.data.2$net.effect) > 
#                                         net.plot.data.2$null.95.perc, 
#                                       "|net effect| > 99th","|net effect| <= 99th")
# 
# net.plot.data$net.value.flag <- ifelse(net.plot.data$net.effect > 1 |
#                                          net.plot.data$net.effect < -1,
#                                        "|net effect| > 1","|net effect| <= 1")

# -------------------------------------------------------------------------
# figures

# p0 <- ggplot(net.plot.data, aes(x = shortest.path.length, y = abs(net.effect))) + 
#   geom_point() + 
#   facet_wrap(richness~., scales = "free") +
# NULL
# p0

p1 <- ggplot(net.plot.data.2, aes(x = abs(net.effect), y = abs(weighted.communicability))) + 
  geom_point(color = "grey10",size = .75, alpha = .7) + 
  facet_wrap(richness~., scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(y = "weighted communicability", x = "net effect (absolute value)") +
  NULL
# p1

p1.2 <- ggplot(net.plot.data.subset, aes(x = abs(net.effect), y = abs(weighted.communicability))) + 
  geom_point(size = .75, color = "grey10", alpha = .7) + 
  facet_wrap(richness~., scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(y = "weighted communicability", x = "net effect (absolute value)") +
  NULL
# p1.2

# -------------------------------------------------------------------------

ggsave(filename = paste("results/images/simulations/communicability_net_effects.pdf",sep=""),
       plot = p1,
       device = cairo_pdf,
       width = 5,height = 3,dpi = 100)

ggsave(filename = paste("results/images/simulations/communicability_net_effects_subset.pdf",sep=""),
       plot = p1.2,
       device = cairo_pdf,
       width = 8,height = 5,dpi = 300)

# -------------------------------------------------------------------------
# cor(df1$net.effect,df1$weighted.communicability)

# metrics.grouped <- net.metrics.df %>%
#   group_by(richness) %>%
#   summarise(cor.net.w = cor(abs(net.effect),abs(weighted.communicability)),
#             cor.net.bin = cor(abs(net.effect),abs(binary.communicability)),
#             cor.dir.w = cor(abs(direct.effect),abs(weighted.communicability)),
#             cor.dir.bin = cor(abs(direct.effect),abs(binary.communicability)),)

