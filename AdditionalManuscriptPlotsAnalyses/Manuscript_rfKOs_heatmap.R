# set up
setwd("Desktop/xanthahumol-master/")
library("dplyr")

# input files
thing1 <- readRDS("Saved_objects/dt_allCovar_randForests_imptKOs_sig_interactions.rds") %>%
  mutate("CovarKO" = paste0(Covar, "_", KO))
thing2 <- read.csv("../Tab4.csv") %>%
  mutate("CovarKO" = paste0(Cognitive.measure.associated.with, "_", KO))

# combine and clean
combothing <- thing1 %>%
  inner_join(y = thing2, by = "CovarKO")
# are there duplicates?
length(unique(combothing$CovarKO)) == nrow(combothing)
length(unique(thing1$CovarKO)) == nrow(combothing)
length(unique(thing2$CovarKO)) == nrow(combothing)
#yeah, from both files; remove them
combothing <- combothing %>% 
  group_by(CovarKO) %>%
  filter(!n() > 1)

# make heatmap
library(ggplot2)

ggplot(combothing, aes(Covar, KO.description)) +
  geom_tile(aes(fill = Intrxn.pval), colour = "white") +
  scale_fill_gradient(low = "steelblue", high = "grey") +
  theme_classic() +
  labs(x = "", y = "", fill = "Random Forest\np-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))





