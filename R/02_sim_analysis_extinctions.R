
# analyze and plot extinction scenarios


# -------------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(patchwork)

# -------------------------------------------------------------------------

# param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")
# 
# network.categories <- read.csv2("results/network_gradient_categories.csv")
# landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
# dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
# cell.distances <- read.csv2("results/cell_distances.csv")
# cell.distances$cell_from <- as.character(cell.distances$cell_from)
# cell.distances$cell_to <- as.character(cell.distances$cell_to)

# global communication efficency
gce.df <- read.csv2("results/extinction_sequences/gce_extinction_sequences.csv")
# just to be sure
gce.df <- subset(gce.df, extinction_sequence != "random")
gce.rand.df <- read.csv2("results/extinction_sequences/gce_random_extinctions.csv")

names(gce.df)[which(names(gce.df) == "normalised.gce")] <- "normalised_gce_mean"
gce.df$normalised_gce_sd <- 0

gce.rand.ag <- gce.rand.df %>%
  group_by(extinction_sequence,species_removed) %>%
  summarise(normalised_gce_mean = mean(normalised.gce),
            normalised_gce_sd = sd(normalised.gce))

gce.df.all <- bind_rows(gce.df,gce.rand.ag)

# -------------------------------------------------------------------------
# communicability dataframes
comm.files <- list.files("results/extinction_sequences/communicability",
                         pattern = "comm*",
                         full.names = T,)

comm.list <- list()
comm.bin <- list()
for(i.file in 1:length(comm.files)){
  load(comm.files[i.file]) 
  comm.bin[[i.file]] <- data.frame(bin.com = comm.df$binary.communicability, 
                                   extinction_sequence = comm.df$extinction_sequence,
                                   species_removed = comm.df$species_removed)
  comm.list[[i.file]] <- comm.df %>%
    mutate(shortest.path.length = na_if(shortest.path.length, Inf)) %>%
    group_by(extinction_sequence,species_removed) %>%
    summarise(#binary_communicability_mean = mean(binary.communicability,na.rm = T),
              #binary_communicability_sd = sd(binary.communicability,na.rm = T),
              weighted_communicability_mean = mean(weighted.communicability,na.rm = T),
              weighted_communicability_sd = sd(weighted.communicability,na.rm = T),
              shortest_path_length_mean = mean(shortest.path.length,na.rm = T),
              shortest_path_length_sd = sd(shortest.path.length,na.rm = T))
  
}
# this needs to be done for the whole data - I need to scale all raw values 
# before averaging
bin.com.df <- bind_rows(comm.bin) %>% 
  mutate(scaled.bin.com = scales::rescale(bin.com,to = c(0,1))) %>%
  group_by(extinction_sequence,species_removed) %>%
  summarise(binary_communicability_mean = mean(scaled.bin.com,na.rm = T),
    binary_communicability_sd = sd(scaled.bin.com,na.rm = T))
other.metrics <- bind_rows(comm.list)
node.metrics <- left_join(other.metrics,bin.com.df)

node.metrics.mean <- node.metrics %>% select(extinction_sequence, species_removed,
                                             weighted_communicability_mean,
                                             binary_communicability_mean,
                                             shortest_path_length_mean) %>%
  pivot_longer(cols = c(weighted_communicability_mean,
                        binary_communicability_mean,
                        shortest_path_length_mean), names_to = "metric",values_to = "value") %>%
  mutate(measurement = "mean")
node.metrics.sd <- node.metrics %>% select(extinction_sequence, species_removed,
                                           weighted_communicability_sd,
                                           binary_communicability_sd,
                                           shortest_path_length_sd) %>%
  pivot_longer(cols = c(weighted_communicability_sd,
                        binary_communicability_sd,
                        shortest_path_length_sd), names_to = "metric",values_to = "value") %>%
  mutate(measurement = "sd")
node.metrics.long <- bind_rows(node.metrics.mean,node.metrics.sd) %>%
  # node.metrics.long$metric2 <- 
  mutate(metric = recode(metric,weighted_communicability_mean = "weighted_communicability",
                         weighted_communicability_sd = "weighted_communicability",
                         binary_communicability_mean = "binary_communicability",
                         binary_communicability_sd = "binary_communicability",
                         shortest_path_length_mean = "shortest_path_length",
                         shortest_path_length_sd = "shortest_path_length")) %>%
  pivot_wider(names_from = measurement,values_from = value)

# likewise for the random replicates, but after obtaining the mean and sd
# NOTE: first version

comm.files.rand <- list.files("results/extinction_sequences/communicability/random_replicates",
                         full.names = T)

comm.list.rand <- list()
comm.bin.rand <- list()
for(i.file in 1:length(comm.files.rand)){
  load(comm.files.rand[i.file]) 
  comm.bin.rand[[i.file]] <- data.frame(bin.com = comm.df$binary.communicability, 
                                   extinction_sequence = comm.df$extinction_sequence,
                                   replicate = comm.df$replicate,
                                   species_removed = comm.df$species_removed)
  comm.list.rand[[i.file]] <- comm.df %>%
    mutate(shortest.path.length = na_if(shortest.path.length, Inf)) %>%
    group_by(extinction_sequence,replicate,species_removed) %>%
    summarise(#binary_communicability_mean = mean(binary.communicability,na.rm = T),
      #binary_communicability_sd = sd(binary.communicability,na.rm = T),
      weighted_communicability_mean = mean(weighted.communicability,na.rm = T),
      weighted_communicability_sd = sd(weighted.communicability,na.rm = T),
      shortest_path_length_mean = mean(shortest.path.length,na.rm = T),
      shortest_path_length_sd = sd(shortest.path.length,na.rm = T))
  
}
# this needs to be done for the whole data - I need to scale all raw values 
# before averaging
bin.com.df.rand <- bind_rows(comm.bin.rand) %>% 
  mutate(scaled.bin.com = scales::rescale(bin.com,to = c(0,1))) %>%
  group_by(extinction_sequence,replicate,species_removed) %>%
  summarise(binary_communicability_mean = mean(scaled.bin.com,na.rm = T),
            binary_communicability_sd = sd(scaled.bin.com,na.rm = T))
other.metrics.rand <- bind_rows(comm.list.rand)
node.metrics.rand <- left_join(other.metrics.rand,bin.com.df.rand)

# summarise across replicates - here
# NOTE this is a first version, I should think of better ways of
# propagating uncertainty

node.metrics.rand.2 <- node.metrics.rand %>%
  group_by(extinction_sequence,species_removed) %>%
  summarise(weighted_communicability_sd = sd(weighted_communicability_mean),
            weighted_communicability_mean = mean(weighted_communicability_mean),
            binary_communicability_sd = sd(binary_communicability_mean),
            binary_communicability_mean = mean(binary_communicability_mean),
            shortest_path_length_sd = sd(shortest_path_length_mean),
            shortest_path_length_mean = mean(shortest_path_length_mean))
# 
node.metrics.mean.rand <- node.metrics.rand.2 %>% select(extinction_sequence, 
                                             # replicate,
                                             species_removed,
                                             weighted_communicability_mean,
                                             binary_communicability_mean,
                                             shortest_path_length_mean) %>%
  pivot_longer(cols = c(weighted_communicability_mean,
                        binary_communicability_mean,
                        shortest_path_length_mean), names_to = "metric",values_to = "value") %>%
  mutate(measurement = "mean")
node.metrics.sd.rand <- node.metrics.rand.2 %>% select(extinction_sequence, 
                                           # replicate,
                                           species_removed,
                                           weighted_communicability_sd,
                                           binary_communicability_sd,
                                           shortest_path_length_sd) %>%
  pivot_longer(cols = c(weighted_communicability_sd,
                        binary_communicability_sd,
                        shortest_path_length_sd), names_to = "metric",values_to = "value") %>%
  mutate(measurement = "sd")
node.metrics.long.rand <- bind_rows(node.metrics.mean.rand,node.metrics.sd.rand) %>%
  # node.metrics.long$metric2 <- 
  mutate(metric = recode(metric,weighted_communicability_mean = "weighted_communicability",
                         weighted_communicability_sd = "weighted_communicability",
                         binary_communicability_mean = "binary_communicability",
                         binary_communicability_sd = "binary_communicability",
                         shortest_path_length_mean = "shortest_path_length",
                         shortest_path_length_sd = "shortest_path_length")) %>%
  pivot_wider(names_from = measurement,values_from = value)

node.metrics.long.2 <- bind_rows(node.metrics.long,node.metrics.long.rand)
# -------------------------------------------------------------------------
# network-level metric
net.plot <- ggplot(gce.df.all, aes(x = species_removed, y = normalised_gce_mean, 
                               group = extinction_sequence)) + 
  geom_line(aes(color = extinction_sequence)) + 
  geom_errorbar(aes(x = species_removed, y = normalised_gce_mean,
                                     ymin = normalised_gce_mean - normalised_gce_sd,
                                     ymax = normalised_gce_mean + normalised_gce_sd,
                                     color = extinction_sequence), alpha = .6) +
  geom_point(aes(fill = extinction_sequence), shape = 21) +
  theme_bw() +
  NULL
net.plot

# -------------------------------------------------------------------------
# node-level metrics
metrics.included <- c("weighted_communicability","shortest_path_length")
plot.data <- subset(node.metrics.long.2, metric %in% metrics.included)
pd <- position_dodge(.7)
node.metrics.plot <- ggplot(plot.data, aes(group = extinction_sequence)) +
  geom_ribbon(aes(x = species_removed, ymin = mean - sd, ymax = mean + sd, fill = extinction_sequence), alpha = .2) +
  geom_line(aes(x = species_removed, y = mean, color = extinction_sequence)) + 
  # geom_errorbar(aes(x = species_removed, ymin = mean - sd, ymax = mean + sd, color = extinction_sequence), alpha = .5) +
  geom_point(aes(x = species_removed, y = mean, fill = extinction_sequence), shape = 21) +
  facet_grid(metric~., scales = "free_y")+
  theme_bw() +
  NULL
node.metrics.plot
