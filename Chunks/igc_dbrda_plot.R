n <- "Bray-Curtis"
dbrda.obj <- dbrda.int.select.list[[n]]
permanova <- dbrda.int.anova.list[[n]]
permanova.dt <- tidy(permanova) %>% as.data.table()
sig.spat.vars <- sapply(permanova.dt[p.value <= 0.05]$term %>% str_split(":"), `[`, 2)

biplot.data <- get.biplot.data(ps.igc, dbrda.obj, plot.axes = c(1, 2))
sample.coords <- melt(
  biplot.data$sample.coords, 
  id.vars = c("Sample", "CAP1", "CAP2", "Treatment"), 
  measure.vars = sig.spat.vars,
  variable.name = "Spat.Var"
  )

vector.coords <- biplot.data$vector.coords[str_detect(Variable, paste(sig.spat.vars, collapse = "|"))]
vector.coords[, Spat.Var := sapply(str_split(Variable, ":"), `[`, 2)]
vector.coords[, Variable := sub("Treatment", "", Variable)]
centroid.coords0 <- biplot.data$centroid.coords
centroid.coords0[, Treatment := sub("Treatment", "", Variable)]
centroid.coords <- lapply(sig.spat.vars, function(var) {
  dt <- copy(centroid.coords0)
  dt[, Spat.Var := var]
}) %>% rbindlist()

caption <- paste(
  "dbRDA ordinations of", n, "distances for mouse microbiome samples. Circles represent individual samples, crosses represent group (diet) centroids. Both of these are colored by diet treatment. Black lines indicate vectors of greatest change for interaction terms from the PERMANOVA results (see previous section). Percent variance in", n, "distance explained is presented in parentheses in the axis labels."
)
dbrda.plot <- ggplot(sample.coords, aes(x = CAP1, y = CAP2)) + 
  geom_point(aes(color = Treatment), alpha = 0.6, size = 2) + 
  geom_point(data = centroid.coords, aes(color = Treatment), shape = 3, size = 4) +
  geom_segment(
    data = vector.coords,
    x = 0, y = 0,
    aes(xend = CAP1, yend = CAP2),
    arrow = arrow(length = unit(0.03, "npc"))
  ) + 
  geom_text_repel(
    data = vector.coords,
    aes(label = Variable),
    size = 3
  ) + 
  labs(
    title = paste(n, "- IGCs"),
    x = biplot.data$axes.labs[1], 
    y = biplot.data$axes.labs[2]
  ) + 
  facet_wrap(~ Spat.Var, ncol = 1) + 
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "top")
