
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors

# OUTPUTS

# -------------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(patchwork)
library(lmerTest)
library(glmmTMB)
library(ggpointdensity)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v3.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

netcom.df <- read.csv2("results/sim_network_level_communicability.csv")
species.data.df <- read.csv2("results/sim_species_level_communicability.csv")

# to assign species-level metrics: degree and dispersal distances
local.networks <- read.csv2("results/sim_degree_dist_networks.csv")
local.dispersal <- read.csv2("results/dispersal_kernels.csv")

# species presences/absences
load("results/presence_dataframes.RData")

# -------------------------------------------------------------------------
# recover factors

# TODO assign parameter values to netcom.df
netcom.df <- left_join(netcom.df,network.categories)
netcom.df <- left_join(netcom.df,landscape.categories)

dispersal.exp.rates <- unique(dispersal.categories[,c(2:3)])
netcom.df <- left_join(netcom.df,dispersal.exp.rates)

netcom.df$landscape.category <- as.factor(netcom.df$spatial.autocorr.value)
netcom.df$network.category <- as.factor(netcom.df$poisson.lambda)
netcom.df$dispersal.category <- factor(netcom.df$exponential.rate,
                                       levels = sort(unique(netcom.df$exponential.rate),decreasing = T))

# network.categories <- network.categories$network.category
# landscape.categories <- landscape.categories$landscape.category
# dispersal.categories <- unique(dispersal.categories$dispersal.category)
# network.categories <- network.categories$network.category
# landscape.categories <- landscape.categories$landscape.category
# dispersal.categories <- unique(dispersal.categories$dispersal.category)

replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
netcom.df$scaled.communicability <- scales::rescale(netcom.df$raw.communicability)

# test
netcom.df$scaled.communicability[netcom.df$scaled.communicability == 0] <- 1e-6

# land.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = landscape.category)) + 
#   geom_density_ridges()
# net.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = network.category)) + 
#   geom_density_ridges()
# disp.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = dispersal.category)) + 
#   geom_density_ridges()

netcom.df.long <- pivot_longer(netcom.df, dispersal.category:landscape.category,
                               names_to = "category.type",
                               values_to = "category.level")

netcom.df.long$category.type <- as.factor(netcom.df.long$category.type)
netcom.df.long$category.level <- as.factor(netcom.df.long$category.level)

# for some reason this does not work
# comm.boxplot <- ggplot(netcom.df.long, aes(x = category.level, y = scaled.communicability)) + 
#   geom_boxplot() + 
#   facet_grid(.~category.type, drop = TRUE) + 
#   NULL
# comm.boxplot

disp.plot <- ggplot(subset(netcom.df.long, category.type == "dispersal.category"), 
                    aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "dispersal rate (d)", y = "landscape communicability") +
  ggtitle("a)") +
  NULL
# disp.plot

landscape.plot <- ggplot(subset(netcom.df.long, category.type == "landscape.category"), 
                    aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "landscape spatial\n autocorrelation (c)", y = "") +
  ggtitle("b)") +
  NULL
# landscape.plot
net.plot <- ggplot(subset(netcom.df.long, category.type == "network.category"), 
                         aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "network degree\n distribution (lambda)", 
       y = "") +
  ggtitle("c)") +
  NULL
# net.plot

combined.plot <- disp.plot + landscape.plot + net.plot

# -------------------------------------------------------------------------
# ggsave("results/images/simulation_categories.pdf",
#        plot = combined.plot,
#        device = cairo_pdf,
#        width = 8, height = 4,dpi = 300)

# -------------------------------------------------------------------------

# simulation results should not be analysed with formal statistical tests, 
# as they can be tricked by increasing power. Correlations are appropriate.

# In addition, I was not able to find a well-specified model for the species-level
# analysis, after trying several ideas. But the qualitative trends seem robust.

# land.dist/net.dist/disp.dist

# m1 <- lm(log(scaled.communicability) ~ landscape.category + network.category + dispersal.category, 
#          data = netcom.df)

# kolmogorov-smirnov test is significant, but the visual patterns are ok
# DHARMa::testResiduals(m1)
# summary(m1)

# -------------------------------------------------------------------------

sp.data.scaled <- species.data.df %>%
  mutate(across(degree:presences,scales::rescale)) %>%
  mutate(id = paste(network.category,landscape.category,dispersal.category,sep="")) %>%
  select(id,sp,wc,bc,degree:presences)

sp.data.scaled.2 <- sp.data.scaled %>% 
  group_by(id) %>%
  mutate(scaled.wc = scales::rescale(wc))


# -------------------------------------------------------------------------

correlations.df <- data.frame(id = unique(sp.data.scaled.2$id), 
                               dispersal.distance = NA, 
                               degree = NA,
                               presences = NA)

for(i in 1:nrow(correlations.df)){
  my.data <- subset(sp.data.scaled.2,id == correlations.df$id[i])
  
  correlations.df$dispersal.distance[i] <- cor.test(my.data$scaled.wc,my.data$dispersal.distance,method = "spearman")$estimate
  correlations.df$degree[i] <- cor.test(my.data$scaled.wc,my.data$degree,method = "spearman")$estimate
  correlations.df$presences[i] <- cor.test(my.data$scaled.wc,my.data$presences,method = "spearman")$estimate
  
}
cor.df.long <- pivot_longer(correlations.df,dispersal.distance:presences,names_to = "factor",values_to = "rho")

cor.df.long$factor[cor.df.long$factor == "dispersal.distance"] <- "dispersal distance"
cor.df.long$factor[cor.df.long$factor == "presences"] <- "number of\npresences"

corr.dist.plots <- ggplot(cor.df.long, aes(x = factor, y = rho)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = expression(rho)) +
  NULL
corr.dist.plots

sp.data.long <- sp.data.scaled.2 %>% pivot_longer(degree:presences,names_to = "factor",values_to = "value")

sp.data.long$factor[sp.data.long$factor == "dispersal.distance"] <- "dispersal distance"
sp.data.long$factor[sp.data.long$factor == "presences"] <- "number of presences"

corr.raw.plots <- ggplot(sp.data.long, aes(x = value, y = scaled.wc)) +
  geom_pointdensity() +
  scale_color_continuous(type = "viridis",name = "number of\nobservations") +
  # geom_hex(bins = 1000) +
  # geom_point(alpha = .4) +
  facet_grid(factor~.,scales = "free_y") +
  # scale_fill_continuous(type = "viridis") +
  theme_bw() +
  theme(strip.background = element_blank())+
  NULL
corr.raw.plots

# -------------------------------------------------------------------------

cor.complete <- corr.dist.plots + corr.raw.plots + plot_layout(widths = c(1, 2))

ggsave("results/images/simulations/species_level_communicability_simulations.pdf",
       plot = cor.complete,
       device = cairo_pdf,
       width = 12, height = 6,dpi = 300)

# sp.data.scaled.2$scaled.wc[sp.data.scaled.2$scaled.wc == 0] <- 1e-7
# sp.data.scaled.2$presences[sp.data.scaled.2$presences == 0] <- 1e-7
# 
# t1 <- subset(sp.data.scaled.2,id == sp.data.scaled.2$id[86])

# -------------------------------------------------------------------------

# m.glm <- glm(scaled.wc ~ degree + dispersal.distance + log(presences),
#              data = t1,family = Gamma(link = "log"))
# 
# m.glmm <- glmmTMB(scaled.wc ~ log(presences) + degree + dispersal.distance + (1|id),
#                   data = sp.data.scaled.2,family = Gamma(link = "log"))
# 
# # mrand <- lmer(wc ~ presences + degree + dispersal.distance + (1|id),
# #               data = sp.data.scaled)
# 
# simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.glmm)
# plot(simulationOutput)
# DHARMa::testResiduals(m.glmm)

# sp2 <- subset(sp.data.scaled,degree < 0.6 & presences < 0.75)

# ggplot(t1,aes(y = scaled.wc, x = presences)) +
#   geom_point() +
#   geom_smooth() +
#   NULL
# 
# mg <- glm(wc ~ presences,
#   data = sp2, family = Gamma(link = "log"))
# 
# simulationOutput <- DHARMa::simulateResiduals(fittedModel = mg)
# plot(simulationOutput)
# DHARMa::plotResiduals(simulationOutput, species.data.df$degree)
# 
# DHARMa::testResiduals(mg)
# summary(mg)

# -------------------------------------------------------------------------
# OLD CODE BELOW
# # -------------------------------------------------------------------------
# # idea: take a small sample of each dataframe (say 5-10%), create a df,
# # and run the model on it. Bootstrap the process 100 or 1000 times.
# # This should be theoretically valid, as these are indpendent instances
# # of the same model, thus the interpretation of the coefficients does not vary.
# 
# # i.model <- 1
# # 
# # sample.perc <- .05
# # n.models <- 100
# # my.models <- list()
# 
# my.files <- list.files("results/communicability",full.names = T)
# 
# # for(i.model in 1:n.models){
#   
#   pairs.data <- list()
#   node.data <- list()
#   my.files <- my.files[sample(1:length(my.files),20)]
#   
#   for(i.file in 1:length(my.files)){
#     
#     # get factor info
#     my.net.cat <- substr(my.files[i.file],30,33)
#     my.landscape.cat <- substr(my.files[i.file],35,38)
#     my.dispersal.cat <- substr(my.files[i.file],40,43)
#     
#     # get communicability
#     load(my.files[i.file])
#     
#     # subset
#     valid.pairs <- subset(comm.df, sp1 != sp2 & weighted.communicability != 0)
#     
#     # aggregate per node
#   my.node.data <- comm.df %>% group_by(sp1,cell1) %>%
#     summarise(wc = sum(weighted.communicability),bc = sum(binary.communicability)) %>%
#     mutate(network.category = my.net.cat,
#            landscape.category = my.landscape.cat,
#            dispersal.category = my.dispersal.cat) %>%
#     arrange(cell1,sp1)
#     
#     # extract sample
#     my.sample <- valid.pairs[sample(nrow(valid.pairs),round(nrow(valid.pairs)*sample.perc),F),]
#     clean.sample <- data.frame(network.category = my.net.cat,
#                                landscape.category = my.landscape.cat,
#                                dispersal.category = my.dispersal.cat,
#                                wc = my.sample[,c("weighted.communicability")])
#     
#     pairs.data[[length(pairs.data)+1]] <- clean.sample
#     node.data[[length(node.data)+1]] <- my.node.data 
#   }# for i.file
#   
#   pairs.df <- bind_rows(pairs.data)
#   nodes.df <- bind_rows(node.data)
#   
# # }# for i.model
# 
#   pairs.means <- pairs.df %>%
#     group_by(network.category,landscape.category,dispersal.category) %>%
#     summarise(mean.wc = mean(wc),sd.wc = sd(wc)) 
#   
#   m1 <- lm(wc~network.category + landscape.category + dispersal.category,data = pairs.df)
#   mm1 <- lm(mean.wc ~ network.category + landscape.category + dispersal.category,data = pairs.means)
#   mn1 <- lm(wc ~ network.category + landscape.category + dispersal.category, data = nodes.df)
#   
#   plot.wc.bc <- ggplot(nodes.df) + 
#     geom_point(aes(x = log(bc), y = log(wc))) +
#     theme_bw() +
#     NULL
#   
#   plot.node.land <- ggplot(nodes.df,aes(x = network.category,y = wc)) + 
#     geom_boxplot(aes(fill = landscape.category)) +
#     # facet_grid(dispersal.category~.)+
#     theme_bw() +
#     # ylim(c(0,1)) +
#     NULL

# plot.means.land <- ggplot(result.means) +
#   geom_point(aes(x = landscape.category, y = mean.wc, fill = network.category), shape = 21) + 
#   geom_errorbar(aes(x = landscape.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
#   facet_grid(network.category~dispersal.category)+
#   theme_bw() + 
#   NULL
# 
# plot.means.net <- ggplot(result.means) +
#   geom_point(aes(x = network.category, y = mean.wc, fill = landscape.category), shape = 21) + 
#   geom_errorbar(aes(x = network.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
#   facet_grid(landscape.category~dispersal.category)+
#   theme_bw() + 
#   NULL

# plot.land <- ggplot(result.df,aes(x = landscape.category,y = wc)) + 
#   geom_boxplot(aes(fill = network.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL
# 
# plot.net <- ggplot(result.df,aes(x = network.category,y = wc)) + 
#   geom_boxplot(aes(fill = landscape.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL






