# setwd("U:/all stuff/manuscripts/InPreparation/Jason David NZ networks")
library(Hmsc)

load(file="data/allData.RData") #Y,X,Tr
nd = length(Y)
ny = nrow(Y[[1]])
ns = ncol(Y[[1]])

nn = matrix(NA,nrow=ns,ncol=nd)
for(i in 1:nd){
  nn[,i] = pmin(colSums(Y[[i]]>0,na.rm = T),colSums(Y[[i]]==0,na.rm = T))
}
colnames(nn) = names(Y)
rownames(nn) = colnames(Y[[1]])
head(nn)
nn.min = pmin(pmin(nn[,1],nn[,2]))
sp = which(nn.min==max(nn.min))[1]

models = list()

for(i in 1:nd){
  models[[i]] = Hmsc(Y=as.matrix(1*(Y[[i]][,sp])>0),
                     XData=X,
                     XFormula = ~elevation+precip_ann+precip_seasonality+forest+shrubland+anthropic_habitats,
                     distr = "probit")
}
names(models) = names(Y)
for(i in 1:nd){
  models[[i]] = sampleMcmc(models[[i]],samples=250,transient = 125,nChains = 4)
}

predY = list()
for(i in 1:nd){
  predY[[i]] = computePredictedValues(models[[i]])
}
MF = list()
MF$AUC = matrix(NA,nrow=nd,ncol=nd)
MF$TjurR2 = matrix(NA,nrow=nd,ncol=nd)
rownames(MF$AUC) = paste0("data = ",names(Y))
colnames(MF$AUC) = paste0("model = ",names(Y))
rownames(MF$TjurR2) = paste0("data = ",names(Y))
colnames(MF$TjurR2) = paste0("model = ",names(Y))
for(i in 1:nd){
  for(j in 1:nd){
    MF$AUC[i,j] = evaluateModelFit(hM=models[[i]], predY=predY[[j]])$AUC
    MF$TjurR2[i,j] = evaluateModelFit(hM=models[[i]], predY=predY[[j]])$TjurR2
  }
}
MF
