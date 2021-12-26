
library(tidyverse)

# -------------------------------------------------------------------------

com <- read.csv2("results/communicability_pairwise.csv")
eff <- read.csv2("results/interaction_effects_pairwise.csv")

df1 <- left_join(com,eff)

# plots

# -------------------------------------------------------------------------
# 1 - binary communicability - weighted communicability
p1 <- ggplot(df1,aes(x = scaled.binary.communicability,
                     y = weighted.communicability)) + 
  geom_point(aes(color = diag)) +
  geom_abline(linetype = "dashed", color = "grey", slope = 1) +
  theme_bw() +
  NULL
# p1

# first take: they should definitely be more related. Binary one seems to
# span more range, the weighted one only discriminates two broad groups
# (although it may be the case that that is the true underlying structure?)
# are these the diagonal points? yes they are

# -------------------------------------------------------------------------
# 2 - binary communicability - net effect

# cor(df1$scaled.binary.communicability,df1$net.effect)

p2 <- ggplot(df1,aes(x = scaled.binary.communicability,
                     y = net.effect)) + 
  geom_point(aes(color = diag)) +
  geom_abline(linetype = "dashed", color = "grey", slope = 1) +
  theme_bw() +
  # lims(y = c(-1,1)) +
  NULL
# p2



