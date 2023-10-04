
# Script to plot communicability-path length relationships

# INPUTS: 
# - pairwise communicability and path lengths: "{external_path}/results/pairwise_metrics_combined.Rdata"

# OUTPUTS:
# - several plots: "results/images/*.pdf-png"

# NOTE:
# some inputs are too big for git/github. load them externally and keep 
# the path {external_path} always the same

external_path <- "/home/david/Work/datasets/NZ/"

# -------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(patchwork)

# -------------------------------------------------------------------------

# load(paste(external_path,"results/communicability_pairwise.Rdata",sep=""))
# load(paste(external_path,"results/shortest_path_lengths.Rdata",sep=""))
# 
# df.combined <- left_join(df.path.lengths,df1[,c("sp1","grid.id.sp1",
#                                                 "sp2","grid.id.sp2",
#                                                 "scaled.binary.communicability")])
# 
# save(df.combined,file = paste(external_path,"results/pairwise_metrics_combined.Rdata",sep=""))

load(paste(external_path,"results/pairwise_metrics_combined.Rdata",sep=""))

# -------------------------------------------------------------------------
# plots

metrics.clean <- subset(df.combined, !(is.infinite(shortest.path.length)))

# many pairs have comm = 0 but paths of length > 0. Check those
# tt <- subset(metrics.clean,scaled.binary.communicability == 0)
# table(tt$shortest.path.length)
metrics.clean$shortest.path.length <- as.factor(metrics.clean$shortest.path.length)
path.length.comm.plot <- ggplot(metrics.clean, aes(x = shortest.path.length, 
                                                   y = scaled.binary.communicability)) + 
  geom_hex() +
  NULL

path.length.comm.plot

# # richness
# bird.richness.plot <- ggplot(data = grid.metrics.data) + 
#   geom_sf(aes(fill = birds)) + 
#   scale_fill_viridis_c() +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "") +
#   theme_bw() +
#   NULL
# 
# plant.richness.plot <- ggplot(data = grid.metrics.data) + 
#   geom_sf(aes(fill = plants)) + 
#   scale_fill_viridis_c() +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "") +
#   theme_bw() +
#   NULL
# 
# # propagation
# bird.plot <- ggplot(data = grid.plot.data) + 
#   geom_sf(aes(fill = birds)) + 
#   scale_fill_viridis_c() +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "") +
#   theme_bw() +
#   NULL
# 
# plant.plot <- ggplot(data = grid.plot.data) + 
#   geom_sf(aes(fill = plants)) + 
#   scale_fill_viridis_c() +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "") +
#   theme_bw() +
#   NULL
# 
# metrics.plot <- plant.richness.plot + bird.richness.plot + plot_annotation(tag_levels = 'A')
# prop.plot <- plant.plot + bird.plot + plot_annotation(tag_levels = 'A')
# 
# # -------------------------------------------------------------------------
# 
# ggsave(filename = "results/images/richness_plot.pdf",
#        plot = metrics.plot,
#        device = cairo_pdf,
#        width = 10, height = 6,dpi = 300)
# ggsave(filename = "results/images/richness_plot.png",
#        plot = metrics.plot,
#        # device = cairo_pdf,
#        width = 10, height = 6,dpi = 300)
# 
# ggsave(filename = paste("results/images/propagation_plot_",my.source.cell,"_",my.source.guild,".pdf",sep=""),
#        plot = prop.plot,
#        device = cairo_pdf,
#        width = 10, height = 6,dpi = 300)
# ggsave(filename = paste("results/images/propagation_plot_",my.source.cell,"_",my.source.guild,".png",sep=""),
#        plot = prop.plot,
#        # device = cairo_pdf,
#        width = 10, height = 6,dpi = 300)


