saveDir <- "Saved_objects"
inDir <- "Input"
plotDir <- "Plots"
checkDirs <- lapply(c(saveDir, inDir, plotDir), function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
})

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, dpi = 150, cache = FALSE)
source("packages_sources.R")
source("Helper_scripts/map_kos.R")
source("Helper_scripts/random_forest_and_importance_functions.R")
theme_set(theme_cowplot())
my_theme <- theme_update(
  legend.position = "top",
  legend.box = "vertical",
  legend.box.just = "left",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  strip.text = element_text(size = 10),
  plot.caption = element_text(hjust = 0, size = 10),
  axis.text = element_text(size = 8)
)
maxCores <-ifelse(str_detect(getwd(), "stagamak"), 20, 3)

assign.redo(
  c(
    "make.igc.phyloseq",
    "make.kos.phyloseq",
    "igc.analyses",
    "kos.analyses",
    "igc.diet.random.forest",
    "kos.diet.random.forest",
    "igc.covar.random.forest",
    "kos.covar.random.forest"
  ),
  state = FALSE
)
toggle.redo(c("kos.analyses", "kos.diet.random.forest", "kos.covar.random.forest"))


sample.df <- read.csv(
  file.path(inDir, "SpatialLearning_metadata.csv"),
  header = T,
  row.names = 1
)
sample.df <- sample.df[complete.cases(sample.df), ]
sample.df$Diet <- factor(sample.df$Diet, levels = c("Control", "XN", "DXN", "TXN"))
sample.dt <- as.data.table(sample.df, keep.rownames = "Sample") %>% setkeyv("Sample")
spat.learn.vars <- names(sample.df)[-1]

igc.phyloseq.file <- file.path(saveDir, "phyloseq_with_igc_contig_counts.rds")
ps.igc <- redo.if("make.igc.phyloseq", igc.phyloseq.file, {
  igc.contig.mat0 <- read.table(
    file.path(inDir, "XN_count_table.txt"),
    header = T,
    sep = "\t",
    row.names = 1
  ) %>% as.matrix()
  gene.lens0 <- igc.contig.mat0[, 1]
  igc.contig.mat1 <- igc.contig.mat0[, -1]
  row.sums <- rowSums(igc.contig.mat1)
  col.sums <- colSums(igc.contig.mat1)
  keep.rows <- row.sums[row.sums > 0] %>% names()
  igc.contig.mat2 <- igc.contig.mat1[keep.rows, ]
  gene.lens1 <- gene.lens0[keep.rows]
  igc.contig.mat3 <- igc.contig.mat2 / gene.lens1
  col.names <- colnames(igc.contig.mat3)
  for (name in col.names) {
    igc.contig.mat3[, name] <- igc.contig.mat3[, name] / col.sums[name]
  }
  colnames(igc.contig.mat3) <- gsub("sample_", "s", colnames(igc.contig.mat3))
  igc.contig.mat4 <- t(igc.contig.mat3)
  igc.contig.mat <- igc.contig.mat4[row.names(sample.df), ]
  phyloseq(
    sample_data(sample.df),
    otu_table(igc.contig.mat, taxa_are_rows = FALSE)
  )
})

kos.phyloseq.file <- file.path(saveDir, "phyloseq_with_ko-collapsed_contig_counts.rds")
ps.kos <- redo.if("make.kos.phyloseq", kos.phyloseq.file, {
  igc.mat0 <- otu.matrix(ps.igc)
  map.tbl0 <- read.table(file.path(inDir, "XN_gene_to_ko_mapping.tsv"), header = F, sep = "\t") %>%
    as.data.table()
  names(map.tbl0) <- c("Contig", "GeneID", "KO", "Length")

  map.tbl1 <- map.tbl0[KO != ""]
  multi.kos <- str_subset(map.tbl1$KO, ",")

  cl <- makeCluster(maxCores, type = "FORK", outfile = "")
  registerDoParallel(cl, maxCores)
  new.rows.dt <- foreach(
    multi.ko = multi.kos,
    .final = rbindlist,
    .verbose = TRUE
  ) %dopar% {
    kos <- str_split(multi.ko, ",")[[1]]
    rows <- map.tbl1[KO == multi.ko]
    kos.dt <- lapply(kos, function(ko) {
      rows$KO <- ko
      return(rows)
    }) %>% rbindlist()
    return(kos.dt)
  }
  stopCluster(cl)

  map.tbl2 <- rbind(map.tbl1[!str_detect(KO, ",")], new.rows.dt)
  nrow(map.tbl1)
  nrow(map.tbl2)

  keep.contigs <- unique(map.tbl2$Contig)
  keep.contigs <- keep.contigs[keep.contigs %in% colnames(igc.mat0)]

  map.tbl <- map.tbl2[Contig %in% keep.contigs]
  igc.mat <- igc.mat0[, keep.contigs]

  uniq.kos <- unique(map.tbl$KO)
  kos.mat <- NULL
  for (ko in uniq.kos) {
    cat(ko, sep = "\n")
    ko.contigs <- map.tbl[KO == ko]$Contig
    ko.igcs <- igc.mat[, ko.contigs, drop = F]
    ko.col <- rowSums(ko.igcs, na.rm = TRUE)
    kos.mat <- cbind(kos.mat, ko.col)
    colnames(kos.mat)[ncol(kos.mat)] <- ko
  }

  phyloseq(
    sample_data(sample.df),
    otu_table(kos.mat, taxa_are_rows = FALSE)
  )
})

prev.table <- function(tbl, nrow = 6, ncol = 6) {
  ncol <- ifelse(ncol > ncol(tbl), ncol(tbl), ncol)
  print(tbl[1:nrow, 1:ncol])
}
par.dbrdas <- function(dist.mats, dbrda.frm, sample.data, nCores, verbose = TRUE) {
  cl <- makeCluster(nCores, type = "FORK", outfile = "")
  registerDoParallel(cl, nCores)
  dbrda.list <- foreach(
    n = names(dist.mats),
    .final = function(x) setNames(x, names(dist.mats)),
    .verbose = verbose
  ) %dopar% {
    dist <- dist.mats[[n]]
    dbrda.obj <- capscale(as.formula(dbrda.frm), data = sample.data)
    return(dbrda.obj)
  }
  stopCluster(cl)
  return(dbrda.list)
}

ordi.log <- file.path(saveDir, "ordistep.log")
if (!file.exists(ordi.log)) {
  file.create(ordi.log)
}
par.ordistep <- function(full.dbrdas, selectDirection, nCores, seed = 42, verbose = TRUE) {
  cl <- makeCluster(nCores, type = "FORK")
  registerDoParallel(cl, nCores)
  dbrda.select.list <- foreach(
    n = names(full.dbrdas),
    .final = function(x) setNames(x, names(full.dbrdas)),
    .verbose = verbose
  ) %dopar% {
    possibleDirections <- c("both", "forward", "reverse")
    dbrda0 <- full.dbrdas[[n]]
    set.seed(seed)
    dbrda.select <- try(ordistep(dbrda0, direction = selectDirection), silent = T)
    if ("try-error" %in% class(dbrda.select)) {
      cat(
        paste0(
          "# ", Sys.time(), "\n",
          "\tOrdistep on full model for ", n, " distance with `direction = ",
          selectDirection, "` failed."
        ),
        file = ordi.log,
        sep = "\n",
        append = TRUE
      )
      for (newDirection in possibleDirections[possibleDirections != selectDirection]) {
        cat(
          paste0("\tTrying `direction = ", newDirection, "`..."),
          file = ordi.log,
          sep = " ",
          append = TRUE
        )
        set.seed(seed)
        dbrda.select <- try(ordistep(dbrda0, direction = newDirection), silent = T)
        if ("try-error" %in% class(dbrda.select)) {
          cat(paste0("Failed."), file = ordi.log, sep = "\n", append = TRUE)
        } else {
          cat(paste0("Success"), file = ordi.log, sep = "\n", append = TRUE)
          return(dbrda.select)
        }
      }
    } else {
      return(dbrda.select)
    }
  }
  stopCluster(cl)
  return(dbrda.select.list)
}

par.anova.rda <- function(
  dbrdas,
  by.what = c("term", "margin", "axis"),
  perm.model = c("reduced", "direct", "full"),
  nCores,
  seed = 42,
  verbose = TRUE
) {
  cl <- makeCluster(nCores, type = "FORK", outfile = "")
  registerDoParallel(cl, nCores)
  dbrda.anova.list <- foreach(
    n = names(dbrdas),
    .final = function(x) setNames(x, names(dbrdas)),
    .verbose = verbose
  ) %dopar% {
    dbrda.obj <- dbrdas[[n]]
    set.seed(seed)
    permanova <- anova(dbrda.obj, by = by.what, model = perm.model)
    return(permanova)
  }
  stopCluster(cl)
  return(dbrda.anova.list)
}

print.list <- function(l) {
  to.print <- lapply(names(l), function(n) {
    i <- l[[n]]
    cat(paste("###", n, "###"), sep = "\n")
    print(i)
    cat("############\n", sep = "\n")
  })
}
obj.size <- function(x) {paste("Object size:", format(object.size(x), units = "auto"))}
f.size <- function(x) {paste("File size:", round(file.size(x) / 1024^2, 1), "Mb")}
