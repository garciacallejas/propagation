
library(pROC)
library(tidyverse)

# -------------------------------------------------------------------------
# this script performs the AUC validation

# -------------------------------------------------------------------------

load("data/local_nz_networks/local_networks_empirical_modelled_NZ.Rdata")

networks <- c(paste0("NZ0",as.character(1:9)),paste0("NZ",as.character(10:16)))
n <- length(networks)

#excluding species that were not modelled
full_interaction_df <- full_interaction_df[!is.na(full_interaction_df$metaweb_interaction),]

Y = full_interaction_df[,paste0(networks,"_empirical_interaction")]

#setting data to 0 rather than NA if an interaction was not observed
Y[is.na(Y)]=0

sum(rowSums(Y)[full_interaction_df$metaweb_interaction==0],na.rm = T)
#17 empirical interactions observed even if metaweb_interaction is zero
sum(rowSums(Y)[full_interaction_df$metaweb_interaction==1],na.rm = T)
#408 empirical interactions observed when metaweb_interaction is one

#selecting only cases where metaweb shows one
full_interaction_df = full_interaction_df[full_interaction_df$metaweb_interaction==1,]
Y = as.data.frame(full_interaction_df[,paste0(networks,"_empirical_interaction")])
# Y[is.na(Y)]=0

pred = full_interaction_df[,paste0(networks,"_bird.sp.prob")] * full_interaction_df[,paste0(networks,"_plant.sp.prob")]
names(pred) <- paste0(networks,"_cooccurrence.prob")

auc.list <- list()
for(i.obs in 1:length(Y)){
  for(i.mod in 1:length(pred)){
    auc.list[[length(auc.list)+1]] <- data.frame(observed_network = substr(names(Y)[i.obs],0,4),
                                                 modelled_network = substr(names(pred)[i.mod],0,4),
                                                 AUC = as.numeric(auc(predictor = pred[,i.mod],response = Y[,i.obs])))
    
  }
}
auc.df <- bind_rows(auc.list)

hist.auc <- ggplot(subset(auc.df, observed_network == modelled_network)) + 
  geom_histogram(aes(x = AUC), fill = "grey40", color = "grey20") + 
  # xlim(c(0,1)) +
  theme_bw() +
  ylab("number of local communities") +
  NULL
# hist.auc

# now, for a second figure:
# are the predictions for the same local network 
# a better match than the predictions for any other local network?

# in particular: relative frequency of AUC of same local network predictions
# is higher than that from other locations

freq.df <- data.frame(network_id = networks, freq = NA)
for(i in 1:n){
  my.net <- networks[i]
  my.AUC <- subset(auc.df, observed_network == my.net)
  freq.df$freq[i] <- mean(my.AUC$AUC[my.AUC$modelled_network == my.net] > my.AUC$AUC[-(my.AUC$modelled_network == my.net)])
}

freq.plot <- ggplot(freq.df) +
  geom_histogram(aes(x = freq), fill = "grey40", color = "grey20") + 
  # xlim(c(0,1)) +
  theme_bw() +
  ylab("number of local communities") +
  xlab("r") +
  NULL
# freq.plot

# -------------------------------------------------------------------------

# ggsave("results/images/network_level_AUC.png",
#        plot = hist.auc, width = 5, height = 3)
# 
# ggsave("results/images/spatial_comparison_AUC.png",
#        plot = freq.plot, width = 5, height = 3)



