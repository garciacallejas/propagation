

library(tidyverse)
library(igraph)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

# this is a huge file (~12GB on memory)
load("results/binary_communicability_matrix_2CELLS.Rdata")

bird.sp.cells.wide <- read.csv2("data/bird_cell_presences_2CELLS.csv")
bird.sp <- sort(unique(bird.sp.cells.wide$species))

# -------------------------------------------------------------------------

# path lengths
# creating the graph is quick
graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = binary.comm.matrix,
                                               mode = "undirected",
                                               weighted = "1",diag = FALSE)
# this takes ~10min
path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")

# turning it to df is quick
df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")

df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
df.path.lengths$grid.id.sp1 <- sub(".*-", "", df.path.lengths$Var1)
df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
df.path.lengths$grid.id.sp2 <- sub(".*-", "", df.path.lengths$Var2)  

df.path.lengths$guild.sp1 <- ifelse(df.path.lengths$sp1 %in% bird.sp, "birds", "plants")
df.path.lengths$guild.sp2 <- ifelse(df.path.lengths$sp2 %in% bird.sp, "birds", "plants")

df.path.lengths <- df.path.lengths[,c("sp1","guild.sp1","grid.id.sp1",
                                      "sp2","guild.sp2","grid.id.sp2",
                                      "shortest.path.length")]

# -------------------------------------------------------------------------

# it is unfeasible to save a csv file of this size
save(df.path.lengths,file = "results/shortest_path_lengths_2CELLS.Rdata")
