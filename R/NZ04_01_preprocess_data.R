y = read.csv2("data/sp_observations_dataset_long_10km.csv")
cells = unique(y$cell_id)
ny = length(cells)
species = unique(y$species)
ns = length(species)
datasets = unique(y$dataset)
nd = length(datasets)

Y = NULL
for(i in 1:nd){
  Y1 = matrix(NA,nrow=ny,ncol=ns)
  y1 = subset(y,dataset==datasets[i])
  cells1 = unique(y1$cell_id)
  for(fcell in cells1){
    y2 = subset(y1,cell_id==fcell)
    Y1[which(fcell==cells),]=0
    Y1[which(fcell==cells),match(y2$species,species)] = y2$observations
  }
  colnames(Y1) = species
  rownames(Y1) = cells
  Y1 = Y1[,colSums(Y1,na.rm = T)>0]
  Y[[i]] = Y1
}
names(Y) = datasets

x = read.csv2("data/environmental_factors_long_10km.csv")
covariates = unique(x$environmental_factor)
nc = length(covariates)
XA = matrix(NA, nrow=ny, ncol=nc)
for(fcell in cells){
  x1 = subset(x,cell_id==fcell)
  XA[which(fcell==cells),match(x1$environmental_factor,covariates)] = x1$value
}
colnames(XA) = covariates

x = read.csv2("data/land_use_frequencies_10km.csv")
covariates = unique(x$habitat_type)
nc = length(covariates)
XB = matrix(NA, nrow=ny, ncol=nc)
for(fcell in cells){
  x1 = subset(x,cell.id==fcell)
  XB[which(fcell==cells),match(x1$habitat_type,covariates)] = x1$frequency
}
colnames(XB) = covariates

x = read.csv2("data/years_sampled_dataset_long_10km.csv")
XC = matrix(1, nrow=ny, ncol = nd)
for(i in 1:nd){
  x1 = subset(x,dataset==datasets[[i]])
  XC[match(x1$cell_id,cells),i] = x1$years_sampled
  XC[,i] = log(XC[,i])
}
colnames(XC) = paste0("effort.",datasets)
X=as.data.frame(cbind(XA,XB,XC))

rownames(X) = cells

sel=!is.na(rowSums(X))
X = X[sel,]
for(i in 1:nd){
  Y[[i]] = Y[[i]][sel,]
}

xy = read.csv2("data/NZ_grid_coords_10km.csv")
xy=xy[match(rownames(X),xy$cell_id),]
all(xy$cell_id==rownames(X))
XY=cbind(xy$lon_centroid,xy$lat_centroid)
colnames(XY) = c("x","y")
rownames(XY) = rownames(X)
plot(XY[,1],XY[,2])

tr = read.csv2("data/trait_data.csv")
head(tr)
traits = unique(tr$trait)
nt = length(traits)
Tr = matrix(NA,nrow = ns, ncol = nt+2)
for(sp in species){
  tr1 = subset(tr,species==sp)
  Tr[which(sp==species),1] = unique(tr1$guild)
  Tr[which(sp==species),2] = unique(tr1$status)
  Tr[which(sp==species),match(tr1$trait,traits)+2] = tr1$mean.value
}
colnames(Tr) = c("guild","status",traits)
rownames(Tr) = species
Tr = as.data.frame(Tr)
Tr$guild = as.factor(Tr$guild)
Tr$status = as.factor(Tr$status)
for(j in 3:ncol(Tr)){Tr[,j] = as.numeric(Tr[,j])}

cells = rownames(X)
cells = sprintf("%.5d", as.numeric(cells))
rownames(X) = cells
rownames(XY) = cells
for(i in 1:nd) rownames(Y[[i]]) = cells

save(Y,X,XY,Tr,file="processed data/allData.RData")
