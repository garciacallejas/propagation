
require(rgdal)
require(raster)

NZ <- readOGR("../datasets/spatial_data/coastline/nz-coastlines-and-islands-polygons-topo-150k.shp", stringsAsFactors = F)
# plot(NZ)
# here, one area unit is 10km2
# -------------------------------------------------------------------------
# remove small islands

areas <- lapply(NZ@polygons, 
                function(x) sapply(x@Polygons, function(y) y@area))

# Which are the big polygons (area > 1, e.g. 10km2)?
bigpolys <- lapply(areas, function(x) which(x > 1))
big.positions <- rep(0,length(bigpolys))

# totally sure there is a more elegant way
for(i.pol in 1:length(bigpolys)){
  if(length(bigpolys[[i.pol]]) >= 1 && bigpolys[[i.pol]] >= 1){
    big.positions[i.pol] <- i.pol
  }
}

big.positions <- big.positions[big.positions>0]

# simply subset the object
NZ2 <- NZ[big.positions,]
# plot(NZ2)

writeOGR(NZ2, "../datasets/spatial_data/", "NZ_main_islands", 
         driver = "ESRI Shapefile")
