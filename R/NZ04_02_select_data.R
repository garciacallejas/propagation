set.seed(1)

th = 10
load(file="processed data/allData.RData") #Y,X,XY,Tr
nd = length(Y)
print(names(Y))
ny = nrow(Y[[1]])
for(i in 1:nd){
  Y[[i]] = 1*(Y[[i]]>0)
  Y[[i]] = Y[[i]][,colSums(Y[[i]],na.rm=T)>=th]
}
species = intersect(c(colnames(Y$NVS),colnames(Y$TIER1)),colnames(Y$GBIF))
for(i in 1:nd){
  Y[[i]]=Y[[i]][,is.element(colnames(Y[[i]]),species)]
}
Tr = Tr[match(species,row.names(Tr)),]
ns = length(species)
Y.PO = Y$GBIF
Y.PA = Y.PO
sel = match(colnames(Y$TIER1),colnames(Y.PA))
Y.PA[,sel] = Y$TIER1
sel = match(colnames(Y$NVS),colnames(Y.PA))
Y.PA[,sel] = Y$NVS

cells = rownames(X)
cells = sprintf("%.5d", as.numeric(cells))
rownames(X) = cells
rownames(XY) = cells
rownames(Y.PA) = cells
rownames(Y.PO) = cells

save(Y.PA,Y.PO,Tr,X,XY,file = "processed data/selectedData.RData")
