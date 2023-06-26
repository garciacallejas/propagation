
library(tidyverse)

# -------------------------------------------------------------------------

# setwd("U:/all stuff/manuscripts/InPreparation/Jason David NZ networks")
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
  Y[[i]] = Y1
}
names(Y) = datasets

# -------------------------------------------------------------------------

x = read.csv2("data/environmental_factors_long_10km.csv")
covariates = unique(x$environmental_factor)
nc = length(covariates)
XA = matrix(NA, nrow=ny, ncol=nc)
for(fcell in cells){
  x1 = subset(x,cell_id==fcell)
  XA[which(fcell==cells),match(x1$environmental_factor,covariates)] = x1$value
}
colnames(XA) = covariates

# -------------------------------------------------------------------------

x = read.csv2("data/land_use_frequencies_10km.csv")
covariates = unique(x$habitat_type)
nc = length(covariates)
XB = matrix(NA, nrow=ny, ncol=nc)
for(fcell in cells){
  x1 = subset(x,cell.id==fcell)
  XB[which(fcell==cells),match(x1$habitat_type,covariates)] = x1$frequency
}
colnames(XB) = covariates

# -------------------------------------------------------------------------
# correlations between different sets of covariates

land.use.subset <- c("forest","shrubland","anthropic_habitats")
XB.subset <- XB[,land.use.subset]
# XB.subset <- na.omit(XB.subset)
# corland <- cor(XB[,land.use.subset])
# corrplot::corrplot(corland)

env.subset <- c("elevation","precip_ann","precip_seasonality")
XA.subset <- XA[,env.subset]
# XA.subset <- na.omit(XA.subset)

corenv <- cor(XA.subset)
# corrplot::corrplot(corenv)

X.all <- as.data.frame(cbind(XA.subset,XB.subset))
rownames(X.all) = cells

# X.all.nona <- na.omit(X.all)
# cor.all <- cor(X.all.nona)
# corrplot::corrplot(cor.all)

# -------------------------------------------------------------------------

X=as.data.frame(cbind(XA,XB))
rownames(X) = cells

sel=!is.na(rowSums(X))
X = X[sel,]
for(i in 1:nd){
  Y[[i]] = Y[[i]][sel,]
}

tr = x = read.csv2("data/trait_data.csv")
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

save(Y,X,Tr,file="data/allData.RData")
