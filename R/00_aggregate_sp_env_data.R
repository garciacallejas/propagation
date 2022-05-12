
# tidy species and environment data

# INPUTS
# - species observations

# OUTPUTS
# - site x observations dataframe

# -------------------------------------------------------------------------

# cell size in meters
grid.size <- 100000

# where are species observations
sp.path <- "results/sp_observations/"

# -------------------------------------------------------------------------

sp.files <- list.files(path = sp.path, 
                       pattern = paste(grid.size/1e3,"km",sep=""),
                       full.names = T)


# -------------------------------------------------------------------------

tidy.obs <- list()

for(i.sp in 1:length(sp.files)){
  
  my.file <- read.csv2(sp.files[i.sp])
  
  tidy.obs[[i.sp]] <- my.file %>%
    group_by(species,cell_id) %>%
    summarise(observations = n())
  
}

obs.df <- bind_rows(tidy.obs)

# -------------------------------------------------------------------------

write.csv2(obs.df, paste("results/sp_observations_",grid.size/1e3,"km.csv",sep=""))


