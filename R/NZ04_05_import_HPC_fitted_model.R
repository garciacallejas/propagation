library(Hmsc)
library(jsonify)

nChains = 4
samples = 250
fp = "64"
for(thin in c(1,10,100)){
  load("models/unfitted_models.RData")
  for(i in 1:length(models)){
    m = models[[i]]
    importFromHPC = from_json(readRDS(file = paste0("models/post_file_",names(models[i]),"_",as.character(thin),"_",fp,".rds"))[[1]])
    postList = importFromHPC[1:nChains]
    cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))
    transient = round(0.5*thin*samples)
    models[[i]] = importPosteriorFromHPC(m, postList, samples, thin, transient)
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   "_fp_",fp,
                   ".Rdata",sep = "")
  save(models,file=filename)
}
