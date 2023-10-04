
# simulate small communities and calculate communicability
# Supplementary section "communicability in small metacommunities"

# -------------------------------------------------------------------------
library(tidyverse)
source("R/communicability_network.R")
source("R/communicability.R")
# source("R/horizontal_community_matrix.R")
# -------------------------------------------------------------------------
# three patches, three configurations. This needs to be handcrafted

S <- 4
n.blocks <- 3

connected.community <- matrix(1, nrow = S, ncol = S)
diag(connected.community) <- 0

full.block <- matrix(0,nrow = S*n.blocks, ncol = S*n.blocks)
full.block[1:S,1:S] <- connected.community
full.block[(S+1):(2*S),(S+1):(2*S)] <- connected.community
full.block[(2*S+1):(3*S),(2*S+1):(3*S)] <- connected.community

# set dispersal links
for(i in 1:nrow(full.block)){
  for(j in 1:ncol(full.block)){
    if(i%%S == j%%S & i!=j){
      full.block[i,j] <- 1
      full.block[j,i] <- 1
    }
  }
}

# -------------------------------------------------------------------------
# modified local structure

local.block <- full.block

local.block[2,1] <- local.block[1,2] <- 0
local.block[3,2] <- local.block[2,3] <- 0
local.block[4,1] <- local.block[1,4] <- 0

local.block[5,6] <- local.block[6,5] <- 0
local.block[8,5] <- local.block[5,8] <- 0
local.block[7,8] <- local.block[8,7] <- 0

local.block[9,11] <- local.block[11,9] <- 0
local.block[11,10] <- local.block[10,11] <- 0
local.block[10,12] <- local.block[12,10] <- 0

# -------------------------------------------------------------------------
# modified dispersal
dispersal.block <- matrix(0,nrow = S*n.blocks, ncol = S*n.blocks)
dispersal.block[1:S,1:S] <- connected.community
dispersal.block[(S+1):(2*S),(S+1):(2*S)] <- connected.community
dispersal.block[(2*S+1):(3*S),(2*S+1):(3*S)] <- connected.community

# set dispersal links
for(i in 1:nrow(dispersal.block)){
  for(j in 1:ncol(dispersal.block)){
    if(i%%S == 2 & j%%S == 2 & i!=j){
      dispersal.block[i,j] <- 1
      dispersal.block[j,i] <- 1
    }
  }
}

# -------------------------------------------------------------------------
# modified sp composition

composition.block <- full.block

composition.block[3,] <- composition.block[,3] <- 0
composition.block[4,] <- composition.block[,4] <- 0

composition.block[5,] <- composition.block[,5] <- 0
composition.block[7,] <- composition.block[,7] <- 0

composition.block[9,] <- composition.block[,9] <- 0
composition.block[10,] <- composition.block[,10] <- 0

# -------------------------------------------------------------------------

fullc <- communicability_network(full.block,weighted = F)
dispc <- communicability_network(dispersal.block,weighted = F)
localc <- communicability_network(local.block,weighted = F)
compc <- communicability_network(composition.block,weighted = F)

com.df <- data.frame(network = c("full","decreased local\ninteractions",
                                  "decreased\ndispersal",
                                  "decreased\ncomposition"),
                     communicability = c(fullc[[1]],dispc[[1]],localc[[1]],compc[[1]]),
                     maximal.communicability = c(fullc[[2]],dispc[[2]],localc[[2]],compc[[2]]))
com.df$normalised.communicability <- com.df$communicability/com.df$maximal.communicability

# -------------------------------------------------------------------------
com.df$log.maximal.communicability <- log(com.df$maximal.communicability)
com.df$maximal.communicability <- NULL
com.df.long <- pivot_longer(com.df,communicability:log.maximal.communicability,names_to = "metric",values_to = "value")
com.df.long$network <- factor(com.df.long$network,levels = c("full","decreased local\ninteractions",
                                                             "decreased\ndispersal",
                                                             "decreased\ncomposition"))

comp <- ggplot(com.df.long, aes(x = network, y = value)) + 
  geom_col() + 
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  labs(x = "") +
  theme(strip.background = element_blank())+
  NULL
# comp

# ggsave(filename = paste("results/images/simulations/example_networks_metrics.pdf",sep=""),plot = comp,
#        device = cairo_pdf,
#        width = 11,height = 5,dpi = 300)

# -------------------------------------------------------------------------

local.matrix <- matrix(0,nrow = 4, ncol = 4, dimnames = list(c("sp1","sp2","sp3","sp4"),
                                                             c("sp1","sp2","sp3","sp4")))

mat.template <- local.matrix %>% 
  as.data.frame() %>%
  rownames_to_column("sp") %>%
  pivot_longer(-c(sp), names_to = "spcol", values_to = "interaction") %>%
  mutate(sp = factor(sp, levels = c("sp4","sp3","sp2","sp1"))) %>%
  ggplot(aes(x=spcol, y=sp)) + 
  geom_tile(color = "grey20",fill = "white") + 
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.line = element_blank()) +
  NULL

# ggsave(filename = paste("results/images/simulations/matrix_template.pdf",sep=""),plot = mat.template,
#        device = cairo_pdf,
#        width = 3,height = 3,dpi = 300)
