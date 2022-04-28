
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors


# OUTPUTS

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

# -------------------------------------------------------------------------
# recover factors

network.categories <- network.categories$network.category
landscape.categories <- landscape.categories$landscape.category
dispersal.categories <- unique(dispersal.categories$dispersal.category)
replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
# idea: take a small sample of each dataframe (say 5-10%), create a df,
# and run the model on it. Bootstrap the process 100 or 1000 times.
# This should be theoretically valid, as these are indpendent instances
# of the same model, thus the interpretation of the coefficients does not vary.

i.model <- i.file <- 1

sample.perc <- .05
n.models <- 100
my.models <- list()

my.files <- list.files("results/communicability",full.names = T)

for(i.model in 1:n.models){
  
  my.data <- list()
  
  for(i.file in 1:length(my.files)){
    
    load(my.files[i.file])
    unique.pairs <- subset(comm.df, sp1 < sp2)
    
    my.sample <- unique.pairs[sample(nrow(unique.pairs),round(nrow(unique.pairs)*sample.perc),F),]
    
  }# for i.file
  
}# for i.model






