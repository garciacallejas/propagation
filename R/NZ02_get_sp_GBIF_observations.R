
#-------------------------------------------------------------------------

library(tidyverse)
library(rgbif)
library("rnaturalearth")
library(sf)
library(ggspatial)
library(countrycode)
library(CoordinateCleaner)

# -------------------------------------------------------------------------

# set time interval for observations - those older than this year will be discarded
year.first.obs <- 2000

# path to download observations
gbif.path <- "results/sp_observations/gbif"

# gbif user
gbif.user <- "david.garcia.callejas"
gbif.pwd <- "DIY1337TRa"
gbif.email <- "david.garcia.callejas@gmail.com"

# -------------------------------------------------------------------------

# species names in Peralta dataset (already clean)

sp.list <- read.csv2("data/species_list.csv")
all.sp <- sort(unique(sp.list$species))

# 1 - get gbif taxon ids and species names

presences<- as.vector(NULL)
spkeys <- as.vector(NULL)

sp.presences.list <- list()
for (i in all.sp){
  if(i == "Dacrydium_cupressinum"){
    tax.backbone <- name_backbone(name=i, rank='species',kingdom = "plants")
  }else{
    tax.backbone <- name_backbone(name=i, rank='species')
  }
  # cat(i,"-",tax.backbone$usageKey,"\n")
  key<-tax.backbone$usageKey
  dat.search<- occ_count(taxonKey = key,country = "NZ")
  
  sp.presences.list[[length(sp.presences.list)+1]] <- 
    data.frame(species = i, observations = dat.search, key = key)
  
  # presences <-c(presences,dat.search)
  # spkeys <- c(spkeys,key)
}

# df of taxon ids
sp.presences <- bind_rows(sp.presences.list)
# sp.presences <- data.frame(species = all.sp,
#                            observations = presences,
#                            key = spkeys)

# -------------------------------------------------------------------------
# send data request
res <- occ_download(pred_in("taxonKey", sp.presences$key), 
                    pred("hasCoordinate", TRUE), 
                    pred("country","NZ"),
                    pred_gte("year",year.first.obs),
                    format = "SIMPLE_CSV",
                    user = gbif.user,
                    email = gbif.email,
                    pwd = gbif.pwd) 

# estado de la descarga
occ_download_meta(res)

still_running <- TRUE
status_ping <- 30 #seconds

while (still_running) {
  meta <- occ_download_meta(res)
  status <- meta$status
  print(status)
  still_running <- status %in% c("PREPARING", "RUNNING")
  Sys.sleep(status_ping) # Suspender ejecución expresiones R durante intervalo de tiempo determinado.
}

# descargar
down.key <- occ_download_meta(res)$key 
dat <- occ_download_get(key = down.key, path = gbif.path, overwrite = TRUE)

# load in R
dat.imp <- occ_download_import(dat,select = c("basisOfRecord",
                                              "species","decimalLongitude",
                                              "decimalLatitude",
                                              "coordinateUncertaintyInMeters",
                                              "year","countryCode",
                                              "gbifID","institutionCode"))

# -------------------------------------------------------------------------
# filter, clean up

#convert country code from ISO2c to ISO3c
dat.imp$countryCode <- countrycode(dat.imp$countryCode, 
                                 origin =  'iso2c', 
                                 destination = 'iso3c')

#flag problems
# dat.imp <- data.frame(dat.imp)
flags <- clean_coordinates(x = dat.imp, 
                           lon = "decimalLongitude", lat = "decimalLatitude",
                           countries = "countryCode", 
                           species = "species",
                           tests = c("capitals", "centroids", "equal",
                                     "gbif", "institutions",
                                     "zeros")) # most test are on by default

# summary(flags)

# exclude problematic records
dat_cl <- dat.imp[flags$.summary,]

# uncertainty of spatial coordinates
# hist(dat_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)

dat_cl2 <- dat_cl %>% 
    dplyr::filter(coordinateUncertaintyInMeters / 1000 <= 100 | 
                    is.na(coordinateUncertaintyInMeters))

dat_cl2$species <- str_replace(dat_cl2$species," ","_")

# -------------------------------------------------------------------------

# store in csv
write.csv2(dat_cl2,"results/species_observations_GBIF_2000.csv",row.names = F)

# citation
cit <- gbif_citation(dat)
cit$download #Citación global
# length(cit$datasets) 
saveRDS(cit, paste0(gbif.path,"/full_sp_citation_2000.Rdata"))


