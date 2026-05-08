# map_kos.R
# take a phyloseq object with a KO table and map the KOs to modules and pathways with KEGGREST

require(KEGGREST)
require(phyloseq)
require(data.table)

map.kos <- function(ps) {
  kegg.kos <- keggList("ko")
  kegg.kos.dt <- data.table(
    "ko" = sub("ko:", "", names(kegg.kos)),
    "ko.name" = kegg.kos
  )
  setkey(kegg.kos.dt, ko)

  link.kos.mod <- keggLink("module", "ko")
  kegg.mods <- keggList("module")

  kegg.mods.dt <- data.table(
    "mod" = sub("md:", "", names(kegg.mods)),
    "mod.name" = kegg.mods
  )
  setkey(kegg.mods.dt, mod)

  kegg.link.dt <- data.table(
    "ko" = sub("ko:", "", names(link.kos.mod)),
    "mod" = sub("md:", "", link.kos.mod)
  )
  setkey(kegg.link.dt, mod)

  kegg.dt0 <- kegg.link.dt[kegg.mods.dt]
  setkey(kegg.dt0, ko)
  kegg.dt <- kegg.dt0[kegg.kos.dt]
  kegg.dt <- kegg.dt[, .(ko, ko.name, mod, mod.name)]
  setkey(kegg.dt, ko)
  # View(kegg.dt)

  my.kos.dt <- copy(kegg.dt)[taxa_names(ps)]
  saveRDS(my.kos.dt, file = file.path(saveDir, "my_ko_names_and_mod_name.rds"))
  return(my.kos.dt)
}
