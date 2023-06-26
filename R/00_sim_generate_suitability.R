
# script to generate suitability values [0,1] for S species

# INPUTS
# - richness
# - optimum values per species (currently generated randomly)
# - standard deviation around the optimum (currently 0.1 for all sp)

# OUTPUTS
# - dataframe with optimum and deviation values for every sp: "results/sp_suitability.csv"
# -------------------------------------------------------------------------
param <- read.csv2("results/sim_landscape_matrices/parameters_v3.csv")

richness <-  param$richness
sp.names <- paste("sp",sprintf("%02d", 1:richness),sep="")

# -------------------------------------------------------------------------

mean.optimum <- runif(richness)
sd.optimum <- 0.2

suit.df <- data.frame(sp = sp.names, optimum = mean.optimum,sd = sd.optimum)

# -------------------------------------------------------------------------

write.csv2(suit.df, "results/sp_suitability.csv",row.names = F)
