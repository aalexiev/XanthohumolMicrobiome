dbrda.int.select.rds <- file.path(saveDir, "list_igc_selected_interactions_dbRDAs.rds")

dbrda.int.select.list <- redo.if("igc.analyses", save.file = dbrda.int.select.rds, {
  par.ordistep(
    full.dbrdas = dbrda.int0.list,
    selectDirection = "both",
    nCores = maxCores
  )
})

dbrda.int.anova.rds <- file.path(saveDir, "list_igc_interaction_dbRDA_permanovas.rds")

dbrda.int.anova.list <- redo.if("igc.analyses", save.file = dbrda.int.anova.rds, {
  par.anova.rda(
    dbrda.int.select.list, 
    by.what = "margin", 
    nCores = maxCores
    )
})

print.list(dbrda.int.anova.list)

