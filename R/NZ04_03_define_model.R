library(Hmsc)
load("processed data/selectedData.RData") #X,XY,Tr,Y.PA,Y.PO

ny = nrow(Y.PA)
ns = ncol(Y.PA)

studyDesign = data.frame(cell = as.factor(c(row.names(XY),row.names(XY))), method = as.factor(c(rep("PA",ny),rep("PO",ny))))

YY = rbind(Y.PA,Y.PO)
X.PA = X
X.PO = X
X.PA$effort = (X$effort.NVS+X$effort.TIER1)/2
X.PO$effort = X$effort.GBIF
X.PA$method = "PA"
X.PO$method = "PO"
XX = rbind(X.PA,X.PO)
XX$method = as.factor(XX$method)

sel = rowSums(is.na(YY))==0
YY = YY[sel,]
XX = droplevels(XX[sel,])
studyDesign = droplevels(studyDesign[sel,])

sKnot = constructKnots(sData = XY, nKnots = 13,minKnotDist = 50000)
sKnot[,1] = 10000*round(sKnot[,1]/10000)
sKnot[,2] = 10000*round(sKnot[,2]/10000)+5000
plot(XY)
points(sKnot,col="red")
rL.cell = HmscRandomLevel(sData=XY, sMethod = "GPP", sKnot = sKnot)

min(dist(rbind(XY,as.matrix(sKnot))))

XFormula = ~elevation + precip_ann + precip_seasonality +
  forest + shrubland + anthropic_habitats + method*effort

plot(data.frame(XX$elevation,XX$precip_ann,XX$precip_seasonality,XX$forest,
                XX$shrubland,XX$anthropic_habitats,XX$method,XX$effort))

JSDM = Hmsc(Y=YY, XData=XX,
            XFormula = XFormula,
            studyDesign = studyDesign,
            ranLevels = list(cell=setPriors(rL.cell,nfMax=10)),
            distr = "probit")

models = list(JSDM = JSDM)
save(models, file = "models/unfitted_models.RData")
JSDM
