# discarded_code_that_may_come_in_handy_later.R

cl <- makeCluster(maxCores, type = "FORK")
registerDoParallel(cl, maxCores)
dbrda.plot.list <- foreach(
  n = names(dbrda.select.list),
  .final = function(x) setNames(x, names(dbrda.select.list)),
  .verbose = FALSE
) %dopar% {
  dbrda.obj <- dbrda.select.list[[n]]
  permanova <- dbrda.anova.list[[n]]
  biplot.data <- get.biplot.data(ps.good, dbrda.obj, permanova)
  return(
    ggplot(biplot.data$sample.coords, aes(x = CAP1, y = CAP2)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_segment(
        data = biplot.data$vector.coords,
        x = 0, y = 0,
        aes(xend = CAP1, yend = CAP2),
        arrow = arrow(length = unit(0.03, "npc"))
      ) +
      geom_text_repel(
        data = biplot.data$vector.coords,
        aes(label = Variable),
        size = 3
      ) +
      labs(
        title = n,
        x = biplot.data$axes.labs[1],
        y = biplot.data$axes.labs[2]
      )
  )
}
stopCluster(cl)
plot_grid(plotlist = dbrda.plot.list, align = "hv", ncol = 1)

# ```{r igc-dbrda-plot-wIGCscores, echo=FALSE, fig.height=8}
n <- "Bray-Curtis"
igc.mat <- otu.matrix(ps.igc)
dbrda.obj <- dbrda.int.select.list[[n]]

igc.coords.file <- file.path(saveDir, "dt_igc_dbrda_coords.rds")
igc.coords <- redo.if("igc.analyses", save.file = igc.coords.file, {
  update(dbrda.obj, comm = igc.mat) %>%
    scores(choices = c(1, 2), display = "species") %>%
    divide_by(biplot.data$coord.scale) %>%
    as.data.table(keep.rownames = "IGC") %>%
    setkeyv("IGC")
})
qt <- 0.00001
quants <- igc.coords[
  ,
  .(sapply(.SD, function(x) quantile(x, probs = c(qt, 1-qt)))),
  .SDcols = c("CAP1", "CAP2")
]
igc.extrm.coords0 <- igc.coords[
  (CAP1 <= quants$CAP1[1] | CAP1 >= quants$CAP1[2])
  | (CAP2 <= quants$CAP2[1] | CAP2 >= quants$CAP2[2])
] %>% setkeyv("IGC")
igc.extrm.coords0[, ID := 1:nrow(igc.extrm.coords0)]
igc.extrm.coords <- igc.eggnog.map[igc.extrm.coords0]
setkey(igc.extrm.coords, ID)
caption <- paste0(
  "dbRDA ordinations of ", n, " distances for mouse microbiome samples. Circles represent individual samples, crosses represent group (diet) centroids. Both of these are colored by diet treatment. Black lines indicate vectors of greatest change for interaction terms from the PERMANOVA results (see previous section). Red numbers indicate 'species' scores for individual IGCs at the ", qt * 100, "th and ", 100 - (qt * 100), "th percentiles of either CAP score (see table below for what the IDs correspond to). Percent variance in Bray-Curtis distance explained is presented in parentheses in the axis labels."
)
print(
  dbrda.plot +
    geom_text(data = igc.extrm.coords, aes(label = ID), size = 2, color = "red") +
    gg_figure_caption(caption = caption)
)
# ```

# ```{r igc-dbrda-plot-wIGCscores-table, echo=FALSE, results='asis'}
gt(igc.extrm.coords[, .(ID, IGC, Preferred_name)], rowname_col = "ID") %>%
  tab_header(title = "IGCs Plot IDs") %>%
  as_raw_html()
# ```

# ```{r igc-dbrda-ceramide-plot-wIGCscores, eval=FALSE, echo=FALSE, fig.height=8}
n <- "Bray-Curtis"
igc.mat <- otu.matrix(ps.igc)
dbrda.obj <- dbrda.int.select.list[[n]]

igc.coords.file <- file.path(saveDir, "dt_igc_dbrda_coords.rds")
igc.coords <- redo.if("igc.analyses", save.file = igc.coords.file, {
  update(dbrda.obj, comm = igc.mat) %>%
    scores(choices = c(1, 2), display = "species") %>%
    divide_by(biplot.data$coord.scale) %>%
    as.data.table(keep.rownames = "IGC") %>%
    setkeyv("IGC")
})
qt <- 0.00001
quants <- igc.coords[
  ,
  .(sapply(.SD, function(x) quantile(x, probs = c(qt, 1-qt)))),
  .SDcols = c("CAP1", "CAP2")
]
igc.extrm.coords0 <- igc.coords[
  (CAP1 <= quants$CAP1[1] | CAP1 >= quants$CAP1[2])
  | (CAP2 <= quants$CAP2[1] | CAP2 >= quants$CAP2[2])
] %>% setkeyv("IGC")
igc.extrm.coords0[, ID := 1:nrow(igc.extrm.coords0)]
igc.extrm.coords <- igc.eggnog.map[igc.extrm.coords0]
setkey(igc.extrm.coords, ID)
caption <- paste0(
  "dbRDA ordinations of ", n, " distances for mouse microbiome samples. Circles represent individual samples, crosses represent group (diet) centroids. Both of these are colored by diet treatment. Black lines indicate vectors of greatest change for interaction terms from the PERMANOVA results (see previous section). Red numbers indicate 'species' scores for individual IGCs at the ", qt * 100, "th and ", 100 - (qt * 100), "th percentiles of either CAP score (see table below for what the IDs correspond to). Percent variance in Bray-Curtis distance explained is presented in parentheses in the axis labels."
)
print(
  dbrda.plot +
    geom_text(data = igc.extrm.coords, aes(label = ID), size = 2, color = "red") +
    gg_figure_caption(caption = caption)
)
# ```

# ```{r igc-dbrda-ceramide-plot-wIGCscores-table, eval=FALSE, echo=FALSE, results='asis'}
gt(igc.extrm.coords[, .(ID, IGC, Preferred_name)], rowname_col = "ID") %>%
  tab_header(title = "IGCs Plot IDs") %>%
  as_raw_html()
# ```
