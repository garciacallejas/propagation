


# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

disp.df <- read.csv2("results/dispersal_kernels.csv")

# -------------------------------------------------------------------------

disp.medians <- disp.df %>% group_by(exponential.rate,replicate) %>%
  summarise(median = median(dispersal.distance)) %>%
  group_by(exponential.rate) %>%
  summarise(avg.median = mean(median))

# -------------------------------------------------------------------------

disp.df <- disp.df %>% group_by(dispersal.category,exponential.rate, replicate) %>% 
  mutate(species.rank = rank(-dispersal.distance,ties.method = "first"))

dp <- ggplot(disp.df, aes(x = species.rank, y = dispersal.distance, 
                          group = replicate)) + 
  geom_line(alpha = .4) +
  facet_grid(.~exponential.rate) +
  labs(x = "species rank", y = "dispersal distance") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  NULL
# dp

ggsave(filename = paste("results/images/dispersal_distances.pdf",sep=""),
       plot = dp,
       device = cairo_pdf,
       width = 8,height = 3,dpi = 300)

