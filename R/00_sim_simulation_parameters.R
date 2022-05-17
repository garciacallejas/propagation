# set inputs for the landscape simulations, so that other files do not need 
# to modify them

# INPUTS
# none

# OUTPUTS
# csv with simulation parameters

# -------------------------------------------------------------------------

sim.parameters <- data.frame(ncol = 0,nrow = 0,
                             num.landscape.categories = 0,
                             richness = 0,
                             num.network.categories = 0,
                             min.lambda = 0,
                             max.lambda = 0,
                             int.mean = 0,
                             int.sd = 0,
                             tau = 0,
                             diag.dom = 0,
                             num.dispersal.categories = 0,
                             min.rate = 0,
                             max.rate = 0,
                             num.category.replicates = 0
                             )

# -------------------------------------------------------------------------
# landscape dimensions
sim.parameters$ncol <- 10
sim.parameters$nrow <- 10

# landscape spatial autocorrelation categories
sim.parameters$num.landscape.categories <- 5

# -------------------------------------------------------------------------
# local communities

sim.parameters$richness <-  30
sim.parameters$num.network.categories <- 5

# poisson mean
sim.parameters$min.lambda <- 3
sim.parameters$max.lambda <- 15 # this should vary with richness. for S = 50, 15 gives connectance = 0.3

# some constants for sampling interaction strengths
sim.parameters$int.mean <- 0
sim.parameters$int.sd <- 1
sim.parameters$tau <- 1.5
sim.parameters$diag.dom <- 0

# -------------------------------------------------------------------------
# dispersal

sim.parameters$num.dispersal.categories <- 5

# exponetial rate
sim.parameters$min.rate <- .75
sim.parameters$max.rate <- .25 

# -------------------------------------------------------------------------
# overall number of replicates per category
sim.parameters$num.category.replicates <- 10

# -------------------------------------------------------------------------
write.csv2(sim.parameters,"results/sim_landscape_matrices/parameters_v2.csv")

