
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors


# OUTPUTS


# -------------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(patchwork)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

netcom.df <- read.csv2("results/sim_network_level_communicability.csv")

# to assign species-level metrics: degree and dispersal distances
local.networks <- read.csv2("results/sim_degree_dist_networks.csv")
local.dispersal <- read.csv2("results/dispersal_kernels.csv")

# species presences/absences
load("results/presence_dataframes.RData")

# -------------------------------------------------------------------------
# recover factors

network.categories <- network.categories$network.category
landscape.categories <- landscape.categories$landscape.category
dispersal.categories <- unique(dispersal.categories$dispersal.category)
replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
netcom.df$scaled.communicability <- scales::rescale(netcom.df$raw.communicability)

land.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = landscape.category)) + 
  geom_density_ridges()
net.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = network.category)) + 
  geom_density_ridges()
disp.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = dispersal.category)) + 
  geom_density_ridges()

# land.dist/net.dist/disp.dist

m1 <- lm(scaled.communicability ~ landscape.category + network.category + dispersal.category, 
         data = netcom.df)

# -------------------------------------------------------------------------
# node-level analyses

comm.files <- list.files("/home/david/Work/datasets/NZ/results/communicability",
                         pattern = "RData",full.names = T)
comm.file.names <- list.files("/home/david/Work/datasets/NZ/results/communicability",
                              pattern = "RData",full.names = F)
# load(my.files[1])
# sample.perc <- .05
# n.models <- 100
species.data.list <- list()

# in case we want to sample only a few files instead of the whole set
# for fast testing, for example
n.files <- length(comm.files)

s <- sample(1:length(comm.files),n.files)
my.files <- comm.files[s]
my.file.names <- comm.file.names[s]

for(j in 1:length(my.files)){
  
  # get factor info
  my.net.cat <- substr(my.file.names[j],6,9)
  my.landscape.cat <- substr(my.file.names[j],11,14)
  my.dispersal.cat <- substr(my.file.names[j],16,19)
  
  my.replicate <- substr(my.file.names[j],23,24)
  if(grepl("_",my.replicate)){
    my.replicate <- substr(my.replicate,1,1)
  }
  my.replicate <- as.numeric(my.replicate)
  
  # number of presences in the grid for each species
  my.presence.df <- sp.presence[[my.landscape.cat]][[my.net.cat]][[my.replicate]]
  my.sp.presences <- my.presence.df %>% 
    group_by(sp) %>%
    summarise(presences = sum(presence))
  
  my.sp.dispersal <- local.dispersal %>% 
    filter(replicate == my.replicate &
             dispersal.category == my.dispersal.cat) %>%
    select(sp,dispersal.distance)
  
  my.sp.degree <- local.networks %>%
    filter(generative.model == my.net.cat & replicate == my.replicate) %>%
    group_by(node_from) %>%
    summarise(degree = sum(value != 0))
  names(my.sp.degree)[1] <- "sp"
  
  my.sp.metrics <- left_join(my.sp.degree,left_join(my.sp.dispersal,my.sp.presences))
  
  # get communicability
  load(my.files[j])
  
  # subset
  valid.pairs <- subset(comm.df, sp1 != sp2 & weighted.communicability != 0)
  
  # aggregate per species
  my.sp.data <- valid.pairs %>% group_by(sp1) %>%
    summarise(wc = sum(weighted.communicability),bc = sum(binary.communicability)) %>%
    mutate(network.category = my.net.cat,
           landscape.category = my.landscape.cat,
           dispersal.category = my.dispersal.cat)
  names(my.sp.data)[1] <- "sp"
  
  species.data.list[[length(species.data.list)+1]] <- left_join(my.sp.data,my.sp.metrics)
  
}# for j (each file)

species.data.df <- bind_rows(species.data.list)
write.csv2(species.data.df,"results/sim_species_level_communicability_TEMP.csv")

# -------------------------------------------------------------------------
sp.data.scaled <- species.data.df %>%
  mutate(across(degree:presences,scales::rescale)) %>%
  select(sp,wc,bc,degree:presences)

# TODO improve
m.sp <- lm(log(wc) ~ degree + 
             dispersal.distance + 
             presences, data = sp.data.scaled)
summary(m.sp)

# -------------------------------------------------------------------------
# OLD CODE BELOW
# # -------------------------------------------------------------------------
# # idea: take a small sample of each dataframe (say 5-10%), create a df,
# # and run the model on it. Bootstrap the process 100 or 1000 times.
# # This should be theoretically valid, as these are indpendent instances
# # of the same model, thus the interpretation of the coefficients does not vary.
# 
# # i.model <- 1
# # 
# # sample.perc <- .05
# # n.models <- 100
# # my.models <- list()
# 
# my.files <- list.files("results/communicability",full.names = T)
# 
# # for(i.model in 1:n.models){
#   
#   pairs.data <- list()
#   node.data <- list()
#   my.files <- my.files[sample(1:length(my.files),20)]
#   
#   for(i.file in 1:length(my.files)){
#     
#     # get factor info
#     my.net.cat <- substr(my.files[i.file],30,33)
#     my.landscape.cat <- substr(my.files[i.file],35,38)
#     my.dispersal.cat <- substr(my.files[i.file],40,43)
#     
#     # get communicability
#     load(my.files[i.file])
#     
#     # subset
#     valid.pairs <- subset(comm.df, sp1 != sp2 & weighted.communicability != 0)
#     
#     # aggregate per node
#   my.node.data <- comm.df %>% group_by(sp1,cell1) %>%
#     summarise(wc = sum(weighted.communicability),bc = sum(binary.communicability)) %>%
#     mutate(network.category = my.net.cat,
#            landscape.category = my.landscape.cat,
#            dispersal.category = my.dispersal.cat) %>%
#     arrange(cell1,sp1)
#     
#     # extract sample
#     my.sample <- valid.pairs[sample(nrow(valid.pairs),round(nrow(valid.pairs)*sample.perc),F),]
#     clean.sample <- data.frame(network.category = my.net.cat,
#                                landscape.category = my.landscape.cat,
#                                dispersal.category = my.dispersal.cat,
#                                wc = my.sample[,c("weighted.communicability")])
#     
#     pairs.data[[length(pairs.data)+1]] <- clean.sample
#     node.data[[length(node.data)+1]] <- my.node.data 
#   }# for i.file
#   
#   pairs.df <- bind_rows(pairs.data)
#   nodes.df <- bind_rows(node.data)
#   
# # }# for i.model
# 
#   pairs.means <- pairs.df %>%
#     group_by(network.category,landscape.category,dispersal.category) %>%
#     summarise(mean.wc = mean(wc),sd.wc = sd(wc)) 
#   
#   m1 <- lm(wc~network.category + landscape.category + dispersal.category,data = pairs.df)
#   mm1 <- lm(mean.wc ~ network.category + landscape.category + dispersal.category,data = pairs.means)
#   mn1 <- lm(wc ~ network.category + landscape.category + dispersal.category, data = nodes.df)
#   
#   plot.wc.bc <- ggplot(nodes.df) + 
#     geom_point(aes(x = log(bc), y = log(wc))) +
#     theme_bw() +
#     NULL
#   
#   plot.node.land <- ggplot(nodes.df,aes(x = network.category,y = wc)) + 
#     geom_boxplot(aes(fill = landscape.category)) +
#     # facet_grid(dispersal.category~.)+
#     theme_bw() +
#     # ylim(c(0,1)) +
#     NULL

# plot.means.land <- ggplot(result.means) +
#   geom_point(aes(x = landscape.category, y = mean.wc, fill = network.category), shape = 21) + 
#   geom_errorbar(aes(x = landscape.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
#   facet_grid(network.category~dispersal.category)+
#   theme_bw() + 
#   NULL
# 
# plot.means.net <- ggplot(result.means) +
#   geom_point(aes(x = network.category, y = mean.wc, fill = landscape.category), shape = 21) + 
#   geom_errorbar(aes(x = network.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
#   facet_grid(landscape.category~dispersal.category)+
#   theme_bw() + 
#   NULL

# plot.land <- ggplot(result.df,aes(x = landscape.category,y = wc)) + 
#   geom_boxplot(aes(fill = network.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL
# 
# plot.net <- ggplot(result.df,aes(x = network.category,y = wc)) + 
#   geom_boxplot(aes(fill = landscape.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL






