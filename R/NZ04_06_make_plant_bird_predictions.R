library(Hmsc)
load(file="models/models_thin_100_samples_250_chains_4.Rdata")
m=models$JSDM
XDataNew=m$XData
for(i in 1:length(XDataNew)){
  print(m$spNames[i])
  print(mean(XDataNew[[i]]$effort))
  print(sqrt(var(XDataNew[[i]]$effort)))  
  XDataNew[[i]]$effort = mean(XDataNew[[i]]$effort)
}
predY = predict(m,XData = XDataNew)

EpredY = Reduce("+",predY)/length(predY)
all(rownames(EpredY)==m$studyDesign$cell)
xy=m$ranLevels$cell$s
all(rownames(EpredY)==rownames(xy))
EpredY.xy = cbind(EpredY,xy)
colnames(EpredY.xy)
plot(EpredY.xy[,"x"],EpredY.xy[,"y"])
save(EpredY.xy,file="plant_bird_predictions.RData")
