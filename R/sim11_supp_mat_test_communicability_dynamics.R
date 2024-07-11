
# compare net effects (dynamics) and communicability (structure)

library(tidyverse)
library(patchwork)
library(ggpointdensity)
library(colorblindr)
library(gghalves)
library(ggbeeswarm)

# source("R/auxiliary_functions/niche_model.R")
source("R/auxiliary_functions/assign_sign_matrix.R")
source("R/auxiliary_functions/assign_interaction_coefs.R")
source("R/auxiliary_functions/GLV_functions.R")
source("R/auxiliary_functions/communicability.R")
source("R/auxiliary_functions/horizontal_community_matrix.R")

set.seed(1984)

# -------------------------------------------------------------------------
# general parameters

net.types <- c("competitive","food web", "mutualistic")

# how many simulated networks
num.nets <- 100
# richness
S <- seq(20,100,by = 20)
# connectance
C <- .25
# connectance of mutualistic systems
C.mut <- .25

# interaction strengths for the generative models
mean.a <- 0.2
sd.a <- 0.1

# parameters specific of the mutualistic systems
mean.b <- .2 # interaction strength same as other interactions
sd.b <- .1 # interaction strength same as other interactions
c.scaling <- 1e-4
alpha.scaling <- 1e-2
b.scaling <- 1e-5
mean.r.plants <- 0.005
sd.r.plants <- 0.005
mean.r.pols <- 0.05
sd.r.pols <- 0.05
N0.plants.scaling <- 1e5
N0.pols.scaling <- 1e5

# -------------------------------------------------------------------------
# signs are the "real" ones in these functions, e.g. competition is negative

# net.metrics.list <- list()
correlations.list <- list()
# i.type <- i.s <- 3

# loop through each interaction type, and generate network replicates
# of varying richness (S)
for(i.type in 1:length(net.types)){
  my.type <- net.types[i.type]
  for(i.s in 1:length(S)){
    my.s <- S[i.s]
    
    # parameters of the dynamical system
    # this is valid for competitive systems and food webs
    # the parameters for mutualistic systems are different
    
    N0 <- rep(1,my.s) # initial conditions for species abundances
    r0 <- rep(1,my.s) # growth rates are all equal
    
    # now, generate interaction matrices and solve 
    # the associated dynamical system
    # the generative model is different depending on the
    # interaction type (random matrix for competition, niche model for
    # food webs, random bipartite matrix for mutualistic)
    # and the mathematical model is the same for 
    # competitive and food webs (standard lineal LV)
    # but different for mutualistic systems (Garcia-Algarra et al. 2014)
    for(i.net in 1:num.nets){
      if(my.type == "competitive"){
        A <- horizontal_community_matrix(S = my.s,
                                         c = C,
                                         min.diag.dom = 1,
                                         int.mean = mean.a,
                                         int.sd = sd.a,
                                         scale_to_1 = T)
        
        # competition must be negative
        parameters <- list(r = r0, A = A * -1)
        
        # solve the dynamical system
        eq.abund <- rootSolve::stode(y = N0,
                                     func=GLV, 
                                     parms = parameters,
                                     positive = TRUE)[[1]]
        
      }else if(my.type == "food web"){
        # A <- niche_model(my.s,C)
        
        niche <- sort(runif(my.s,0,1))
        
        # second connectance is for antagonistic interactions
        sign.matrix <- assign_sign_matrix(S = my.s,
                                          connectances = c(0,C,0,0,0),
                                          antagonism = "structured",
                                          connectance.type = "undirected",
                                          niche = niche)
        
        A <- assign_interaction_coefs(sign.matrix = sign.matrix,
                                      mean.a = mean.a,
                                      sd.a = sd.a,
                                      symmetric = TRUE)
        
        parameters <- list(r = r0, A = A)
        
        # solve the dynamical system
        eq.abund <- rootSolve::stode(y = N0,
                                     func=GLV, 
                                     parms = parameters,
                                     positive = TRUE)[[1]]
        
      }else if(my.type == "mutualistic"){
        
        # use Garcia-Algarra 2014
        # this bipartite system takes two bipartite interaction 
        # matrices, that I join in a unipartite matrix
        # for calculating communicability
        
        num.plants <- my.s/2
        num.pols <- my.s/2
        
        all.links <- numeric(num.plants*num.pols)
        num.links <- C.mut*num.plants*num.pols
        
        popl.links <- all.links
        popl.links[sample(length(popl.links),num.links,replace = F)] <- abs(rnorm(num.links,
                                                                                  mean = mean.b,
                                                                                  sd = sd.b))
        
        B_popl <- matrix(popl.links*b.scaling,nrow = num.plants,ncol = num.pols)
        
        plpo.links <- all.links
        # same links as the other matrix, to have asymmetric but bidirectional interactions
        plpo.links[which(popl.links != 0)] <- abs(rnorm(num.links,
                                                        mean = mean.b,
                                                        sd = sd.b))
        
        B_plpo <- matrix(plpo.links*b.scaling,nrow = num.pols,ncol = num.plants)
        
        # for calculating communicability
        A <- matrix(0,nrow = num.plants+num.pols,ncol = num.plants+num.pols)
        A[1:num.plants,(num.plants+1):(num.plants+num.pols)] <- B_plpo
        A[(num.plants+1):(num.plants+num.pols),1:num.plants] <- B_popl
        
        # growth rates
        rplants <- abs(rnorm(num.plants,mean.r.plants,sd.r.plants))
        rpols <- abs(rnorm(num.pols,mean.r.pols,sd.r.pols))
        alphaplants <- abs(rnorm(num.plants,1,1)*alpha.scaling)
        alphapols <- abs(rnorm(num.pols,1,1)*alpha.scaling)
        cplants <- rep(1*c.scaling,num.plants)
        cpols <- rep(1*c.scaling,num.pols)
        
        parameters <- list(
          rplants = rplants,
          rpols = rpols,
          alphaplants = alphaplants,
          alphapols = alphapols,
          cplants = cplants,
          cpols = cpols,
          B_plpo,B_popl)
        
        Nplants <- runif(num.plants,0,1)*N0.plants.scaling
        Npols <- runif(num.pols,0,1)*N0.pols.scaling
        
        y <- c(Nplants,Npols)
        
        # solve the dynamical system
        eq.abund <- rootSolve::stode(y = y,
                                     func=mutualistic_model, 
                                     parms = parameters,
                                     positive = TRUE)[[1]]
        # eq.abund
        
      }
      
      # if there is a non-negative steady-state solution, it will be eq.abund
      # if there is no positive solution, all eq.abund elements will be zero
      if(sum(eq.abund) > 0){
        
        # calculate jacobian
        # use different model if it's mutualistic or otherwise
        if(my.type == "mutualistic"){
          # I use different dynamic models if mutualistic or other interactions
          J <- rootSolve::jacobian.full(y=eq.abund,func=mutualistic_model,parms = parameters)
        }else{
          J <- rootSolve::jacobian.full(y=eq.abund,func=GLV,parms = parameters)
        }
        
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
                                                     #mode = "undirected",
                                                     mode = "max",
                                                     weighted = "1",diag = FALSE)
      path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
      
      # turning it to df is quick
      df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
      names(df.path.lengths)[1:2] <- c("sp1","sp2")
      
      # generate dataframe
      df1 <- expand.grid(net = i.net,
                         structure = my.type,
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
      
      # net.metrics.list[[length(net.metrics.list)+1]] <- df2
      # obtain correlations right here, right now
      my.cor.df <- df2 %>%
        filter(sp1 != sp2) %>% 
        mutate(effect_type = ifelse(direct.effect != 0,"direct + indirect","indirect")) %>%
        group_by(net, richness, structure, effect_type) %>%
        summarise(spearman.cor = cor(weighted.communicability, net.effect, method = "spearman"))
      correlations.list[[length(correlations.list)+1]] <- my.cor.df
    }# for i.net
  }# for i.s
}
# net.metrics.df <- bind_rows(net.metrics.list)
all.correlations <- bind_rows(correlations.list)

write.csv2("results/supp_mat_communicability_dynamic_effects_correlations.csv",row.names = F)

# -------------------------------------------------------------------------
cor.plot <- ggplot(all.correlations, aes(x = as.factor(richness), y = abs(spearman.cor))) + 
  # geom_half_point(aes(color = effect_type),
  #                 transformation = position_quasirandom(width = 0.1),
  #                 side = "l", size = 0.5, alpha = 0.5) +
  # geom_half_boxplot(aes(fill = effect_type), side = "r",outlier.size = 0.8) +
  geom_boxplot(aes(fill = effect_type)) + 
  facet_grid(cols = vars(structure)) +
  scale_color_OkabeIto(darken = 0.2) +
  scale_fill_OkabeIto(darken = 0.2) +
  theme_bw() +
  theme(strip.background = element_blank())+
  labs(x = "community richness", y = "spearman correlation")
NULL

# -------------------------------------------------------------------------

# ggsave(filename = paste("results/images/simulations/communicability_net_effects.pdf",sep=""),
#        plot = cor.plot,
#        device = cairo_pdf,
#        width = 8,height = 5,dpi = 300)
