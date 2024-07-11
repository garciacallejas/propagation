
library(tidyverse)
library(igraph)
library(bipartite)
library(ggbipart)

# -------------------------------------------------------------------------
# this script plots the plant-frugivore metaweb from Peralta et al 2020
# that is used in Fig. 1 of the manuscript

# -------------------------------------------------------------------------
grid.size <- 10

sp.int.orig <- read.csv2("results/plant_bird_interactions_clean.csv")
metaweb <- sp.int.orig %>% mutate(metaweb_interaction = 1) %>%
  dplyr::select(PLANTSPECIES,BIRDSPECIES,metaweb_interaction) %>% unique()
names(metaweb) <- c("plant.sp","bird.sp","metaweb_interaction")
metaweb$group <- "metaweb"
# both bipartite and ggbipart work with adjacency matrices

sp.cells <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",
                            sep=""))
# this is to harmonise datasets
names(sp.cells)[which(names(sp.cells) == "presence")] <- "observations"
# keep presences only
sp.cells <- subset(sp.cells, observations == 1)

# -------------------------------------------------------------------------
sp.int <- subset(sp.int.orig, PLANTSPECIES %in% sp.cells$species & 
                   BIRDSPECIES %in% sp.cells$species)

plant.sp <- sort(unique(sp.int$PLANTSPECIES))
bird.sp <- sort(unique(sp.int$BIRDSPECIES))
all.sp <- c(plant.sp,bird.sp)

# -------------------------------------------------------------------------
metaweb.clean <- subset(metaweb,plant.sp %in% all.sp & bird.sp %in% all.sp)
adjmat <- frame2webs(dframe = metaweb.clean,varnames = c("plant.sp","bird.sp","group"))[[1]]

# -------------------------------------------------------------------------
# other unused options for plotting:

# ggbipart
# g <- bip_railway(adjmat, label=T) + coord_flip()

# bipartite
# plotweb(adjmat)

# -------------------------------------------------------------------------
# igraph
sp.order <- metaweb.clean %>% group_by(bird.sp) %>% summarise(num.int = n()) %>% arrange(desc(num.int))
metaweb.clean$bird.sp <- factor(metaweb.clean$bird.sp,levels = sp.order$bird.sp)
metaweb.clean <- arrange(metaweb.clean,bird.sp)

my.graph <- graph_from_data_frame(metaweb.clean,directed = F)
V(my.graph)$type <- bipartite_mapping(my.graph)$type
V(my.graph)$color <- ifelse(V(my.graph)$type, "darkgoldenrod2", "darkolivegreen4")
# V(my.graph)$shape <- ifelse(V(g)$type, "circle", "square")
E(my.graph)$color <- "gray"
V(my.graph)$size <- 7

my.graph2 <- my.graph
V(my.graph2)$label <- ""
V(my.graph2)$size <- ifelse(V(my.graph2)$type, 2, 1.3)
V(my.graph2)$frame.color <- "grey30"
V(my.graph2)$frame.size <- .5
E(my.graph2)$color <- "gray50"
E(my.graph2)$width <- .8 # does not seem to work with pdf files
# plot(my.graph2, layout=layout.bipartite,asp=0.25)

# png("results/images/tests_figure_1/bipartite_metaweb.png",1600,800)
# pdf("results/images/tests_figure_1/bipartite_metaweb.pdf",6,3)
# plot(my.graph2, layout=layout.bipartite,asp=0.25)
# dev.off()








