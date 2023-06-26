

# -------------------------------------------------------------------------

library(tidyverse)
library(patchwork)

# -------------------------------------------------------------------------

load("results/landscape_matrices.RData")
autocorr.values <- read.csv2("results/spatial_autocorrelation_categories.csv")

# -------------------------------------------------------------------------
n.rows <- nrow(landscape.list[[1]][[1]])

rep <- 1

plot.list <- list()

for(i in 1:length(landscape.list)){
  my.matrix <- landscape.list[[i]][[rep]]
  
  my.df <- reshape2::melt(my.matrix)
  names(my.df) <- c("y","x","landscape.attribute")
  mean.att <- mean(my.df$landscape.attribute)
  sd.att <- sd(my.df$landscape.attribute)
  
  my.plot <- ggplot(my.df) +
    geom_tile(aes(x = x, y = y, fill = landscape.attribute)) + 
    scale_fill_viridis_c() +
    scale_x_discrete(expand = c(0, 0), breaks = NULL) + 
    scale_y_discrete(expand = c(0,0), breaks = NULL) +
    # scale_y_discrete(breaks=NULL) +
    # scale_x_discrete(breaks=NULL) +
    labs(x="",y = "") +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle(paste(n.rows, "rows\nspatial autocorrelation:",autocorr.values$spatial.autocorr.value[i],
            "\nmean:",round(mean.att,2),"sd:",round(sd.att,2))) +
    NULL
  plot.list[[i]] <- my.plot
}

all.plots <- patchwork::wrap_plots(plot.list,ncol = 3)

# -------------------------------------------------------------------------
ggsave("results/images/landscape_plots.pdf",
       plot = all.plots,
       device = cairo_pdf,
       width = 8, height = 6,dpi = 300)



