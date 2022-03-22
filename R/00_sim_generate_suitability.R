
# script to generate suitability values [0,1] for S species

# INPUTS
# - richness
# - optimum values per species (currently generated randomly)
# - standard deviation around the optimum (currently 0.1 for all sp)

# OUTPUTS
# - dataframe with optimum and deviation values for every sp: "results/sp_suitability.csv"
# -------------------------------------------------------------------------

richness <- 30
sp.names <- paste("sp",1:richness,sep="")

# -------------------------------------------------------------------------

mean.optimum <- runif(richness)
sd.optimum <- 0.1

suit.df <- data.frame(sp = sp.names, optimum = mean.optimum,sd = sd.optimum)

# -------------------------------------------------------------------------

write.csv2(suit.df, "results/sp_suitability.csv",row.names = F)
