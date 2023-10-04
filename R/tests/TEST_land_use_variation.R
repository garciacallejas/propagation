
library(tidyverse)

# -------------------------------------------------------------------------

land.use.1990 <- read.csv2("results/land_use_tests/land_use_1990.csv")
land.use.2016 <- read.csv2("results/land_use_tests/land_use_2016.csv")

categories <- data.frame(category = sort(unique(land.use.2016$LUCID_2016)),
                         category_char = c("Natural forest","Pre-1990 planted forest",
                                           "Post-1989 planted forest", "Grassland with woody biomass",
                                           "Grassland - high producing", "Grassland - low producing",
                                           "Cropland - perennial","Cropland - annual", "Wetland - open water",
                                           "Wetland - vegetated non forest", "Settlements", "Other"
                                           ))

# -------------------------------------------------------------------------
ncells <- unique(land.use.1990$cell_id)
categories <- sort(unique(land.use.2016$LUCID_2016))

all.combinations <- expand.grid(cell_id = ncells,
                                category = categories,
                                variation = NA)

# -------------------------------------------------------------------------

for(i.cell in ncells){
  for(i.cat in 1:length(categories)){
    my.1990.obs <- land.use.1990$rel.area.1990[which(land.use.1990$cell_id == i.cell & 
                                                       land.use.1990$LUCID_1990 == categories[i.cat])]
    my.2016.obs <- land.use.2016$rel.area.2016[which(land.use.2016$cell_id == i.cell & 
                                                       land.use.2016$LUCID_2016 == categories[i.cat])]
    if(length(my.1990.obs) == 0){
      my.1990.obs <- 0
    }
    if(length(my.2016.obs) == 0){
      my.2016.obs <- 0
    }
    
    all.combinations$variation[which(all.combinations$cell_id == i.cell &
                                       all.combinations$category == categories[i.cat])] <- my.2016.obs - my.1990.obs
  }
}

all.combinations.summarised <- all.combinations %>%
  group_by(category) %>%
  summarise(mean.variation = mean(variation, na.rm = T),
            median.variation = median(variation, na.rm = T),
            sd.variation = sd(variation, na.rm = T)) %>%
  left_join(categories)

# -------------------------------------------------------------------------

categories.variation.densities <- ggplot(all.combinations, aes(x = variation)) + 
  geom_density(aes(fill = category)) + 
  facet_grid(category~., scales = "free_y") +
  xlim(-0.05,0.05) +
  NULL
# categories.variation.plot

categories.variation.mean <- ggplot(all.combinations.summarised, aes(x = category_char)) + 
  geom_point(aes(y = mean.variation)) +
  geom_errorbar(aes(ymin = mean.variation - sd.variation, ymax = mean.variation + sd.variation)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0)) +
  ylab("mean variation per cell") + xlab("") +
  # facet_grid(category~., scales = "free_y") +
  # xlim(-0.05,0.05) +
  NULL
# categories.variation.mean
ggsave(filename = "results/images/land_use_tests/land_use_variation_1990_2016.pdf",
       plot = categories.variation.mean,
       device = cairo_pdf,
       width = 7,height = 5,dpi = 300)

