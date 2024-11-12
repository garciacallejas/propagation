library(Hmsc)
library(jsonify)

load(paste0("models/unfitted_models.RData"))

thin = 1
samples = 250
nChains = 4
transient = round(0.5*thin*samples)
verbose = 10

for(i in 1:length(models)){
  m = models[[i]]
  init_obj = sampleMcmc(m, samples=samples, thin=thin,
                        transient=transient, nChains=nChains,
                        verbose=verbose, engine="HPC")
  init_file_path = file.path(getwd(), paste0("models/init_file_",names(models)[i],".rds"))
  saveRDS(to_json(init_obj), file=init_file_path)
}
pre="64" #computational precision
todofile = paste0("models/todo",pre,".txt")
write("#!/bin/bash",file=todofile,append=FALSE)
for(thin in c(1,10,100)){
  for(i in 1:length(models)){
    write(paste0("python3 -m hmsc.run_gibbs_sampler --input ./init_file_",names(models)[i],
                 ".rds --output ./post_file_",names(models)[i],"_",
                 as.character(thin),"_",pre,".rds --samples ",as.character(samples),
                 " --transient ",as.character(round(samples*0.5*thin)),
                 " --thin ",as.character(thin),
                 " --verbose ",as.character(verbose),
                 " --fp ",pre),
          file=todofile,append = TRUE) 
  }
}

