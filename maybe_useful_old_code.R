# setup
bad.smpls <- row.names(sample_data(ps1)[!complete.cases(sample_data(ps1)), ])
ps.good <- prune_samples(!(sample_names(ps1) %in% bad.smpls), ps1)
sample.df <- sample.data.frame(ps.good)
learn.data <- sample.df[, spat.learn.vars]
learn.cors <- cor(learn.data)
high.cor.indices <- which(abs(learn.cors) < 1 & abs(learn.cors) > 0.6, arr.in = TRUE)
correlated.vars <- data.table()
# i <- 1; j <- 2
for (i in 1:(nrow(high.cor.indices) - 1)) {
  row1 <- high.cor.indices[i, ]
  for (j in (i + 1):nrow(high.cor.indices)) {
    row2 <- high.cor.indices[j, ]
    if (row1["row"] == row2["col"] & row1["col"] == row2["row"]) {
      correlated.vars <- rbind(
        correlated.vars,
        data.table(
          Var1 = rownames(high.cor.indices)[i],
          Var2 = rownames(high.cor.indices)[j],
          r = learn.cors[rownames(high.cor.indices)[i], rownames(high.cor.indices)[j]]
        )
      )
    }
  }
}
to.remove <- correlated.vars$Var2
# set.seed(42)
# spat.learn.vars <- sample(names(sample.df)[!(names(sample.df) %in% to.remove)])
spat.learn.vars <- names(sample.df)[!(names(sample.df) %in% to.remove)]
sample.df <- sample.df[, spat.learn.vars]
spat.learn.vars <- spat.learn.vars[-which(spat.learn.vars == "treat")]
sample_data(ps.good) <- sample_data(sample.df)

if (redo$remap.kos) {
  source("Helper_scripts/map_kos.R")
}


# <!-- no significant interactions with diet treatment. Sørensen is behaving oddly with the new set of terms, and I'm not totally sure why. You'd think a reduction in the total input terms would make things more clear, not less so. -->

  ## Looking at diet alone

  ## Generate dbRDA ordinations with all Spatial Learning covariates

  # ```{r diet-dbrda-generation, echo=TRUE}
dbrda.diet.rds <- file.path(saveDir, "list_of_diet_dbRDAs.rds")
dbrda.diet.frm <- "dist ~ treat"

if (analyses) {
  dbrda.diet.list <- par.dbrdas(
    dist.mats = dist.list,
    sample.data = sample.df,
    nCores = maxCores,
    dbrda.frm = dbrda.diet.frm
  )
  saveRDS(dbrda.diet.list, file = dbrda.diet.rds)
} else {
  dbrda.diet.list <- readRDS(dbrda.diet.rds)
}
# ```

## Tests of significance

# ```{r diet-dbrda-anovas}
dbrda.diet.anova.rds <- file.path(saveDir, "list_of_diet_dbRDA_permanovas.rds")

if (analyses) {
  dbrda.diet.anova.list <- par.anova.rda(
    dbrda.diet.list,
    by.what = "term",
    nCores = maxCores
  )
  saveRDS(dbrda.diet.anova.list, file = dbrda.diet.anova.rds)
} else {
  dbrda.diet.anova.list <- readRDS(dbrda.diet.anova.rds)
}

print.list(dbrda.diet.anova.list)
# ```

## Test interaction between Diet and Hid3 (Bray-Curtis only)

# ```{r}
dist.bc <- dist.list[["Bray-Curtis"]]
dbrda <- capscale(dist.bc ~ treat * Hid3, data = sample.df)
print(anova(dbrda, by = "term"))
# ```



# ```{r set-specific-KOs}
focal.kos <- c(
  "K00838", # ARO8; aromatic amino acid aminotransferase I / 2-aminoadipate transaminase [EC:2.6.1.57 2.6.1.39 2.6.1.27 2.6.1.5]
  "K05821", # ARO9; aromatic amino acid aminotransferase II [EC:2.6.1.58 2.6.1.28]
  "K10797", # enr; 2-enoate reductase
  "K01859", # E5.5.1.6; chalcone isomerase [EC:5.5.1.6]
  "K01667", #tnaA; tryptophanase [EC:4.1.99.1]
  "K00466", # iaaM; tryptophan 2-monooxygenase [EC:1.13.12.3]
  "K00766", # trpD; anthranilate phosphoribosyltransferase [EC:2.4.2.18]
  "K01593", # DDC, TDC; aromatic-L-amino-acid/L-tryptophan decarboxylase [EC:4.1.1.28 4.1.1.105]
  # No ILDH, Indolelactate dehydrogenase
  # No ILD, Indolelactate dehydratase
  "K04103", # ipdC; indolepyruvate decarboxylase [EC:4.1.1.74]
  "K00249", # ACADM, acd; acyl-CoA dehydrogenase [EC:1.3.8.7]
  "K16173", # acd; glutaryl-CoA dehydrogenase (non-decarboxylating) [EC:1.3.99.32]
  # No IaaD, Indoleacetate decarboxylase
  # No IaaDH, Indoleacetaldehyde dehydrogenase
  # No IaaR, Indoleacetaldehyde reductase
  "K21801", # iaaH; indoleacetamide hydrolase [EC:3.5.1.-]
  "K07415", #CYP2E1; cytochrome P450 family 2 subfamily E polypeptide 1 [EC:1.14.13.-]
  "K01014", #   SULT1A; aryl sulfotransferase [EC:2.8.2.1]
  "K01015", #   SULT2B1; alcohol sulfotransferase [EC:2.8.2.2]
  "K01016", #   SULT1E1, STE; estrone sulfotransferase [EC:2.8.2.4]
  "K01025", #   SULT1; sulfotransferase [EC:2.8.2.-]
  "K11822", #   SULT2A1; bile-salt sulfotransferase [EC:2.8.2.14]
  "K11823", #   SULT4A1; sulfotransferase 4A1 [EC:2.8.2.-]  NA  NA
  "K16949", #   SULT3A1; amine sulfotransferase [EC:2.8.2.3]
  "K22523" # SULT6B1; sulfotransferase 6B1 [EC:2.8.2.-]
)
gt::gt(unique(kos.mod.dt[ko %in% focal.kos, .(ko.name), keyby = "ko"]))
saveRDS(unique(kos.mod.dt[ko %in% focal.kos, .(ko.name), keyby = "ko"]), file = "Saved_objects/dt_XN_kos_of_interest.rds")
# ```

# ```{r specific-KOs-by-diet, fig.height=10}
keep.cols <- c("Sample", focal.kos)
ko.tbl <- otu.data.table(ps.igc)
ko.tbl <- ko.tbl[, ..keep.cols]
setkey(ko.tbl, Sample)
focal.tbl0 <- ko.tbl[sample.dt]
focal.tbl <- melt(
  focal.tbl0,
  measure.vars = focal.kos,
  variable.name = "KO",
  value.name = "Rel_abund"
)
ggplot(focal.tbl, aes(x = treat, y = Rel_abund)) +
  geom_quasirandom() +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "red", width = 0.8) +
  stat_summary(fun = "mean", geom = "point", size = 2, color = "red") +
  facet_wrap(~ KO, scales = "free_y", ncol = 3) +
  labs(x = "Diet Treatment", y = "Relative Abund.") +
  theme(axis.text.y = element_text(size = 8))
# ```

# ```{r spec-KOs-heatmap-mean, fig.height=8}
heatmap.tbl0 <- copy(focal.tbl0)[, .(Sample, treat, sapply(.SD, scale)), .SDcols = focal.kos]
heatmap.tbl1 <- melt(
  heatmap.tbl0,
  measure.vars = focal.kos,
  variable.name = "KO",
  value.name = "Scaled_abund"
)
heatmap.tbl <- heatmap.tbl1[
  ,
  .(Mean_abund = mean(Scaled_abund), Med_abund = median(Scaled_abund)),
  by = c("treat", "KO")
]
ggplot(heatmap.tbl, aes(y = KO, x = treat)) +
  geom_bin2d(aes(fill = Mean_abund), stat = "identity") +
  scale_fill_gradientn(name = "Mean scaled\nabundance", colors = c("blue", "white", "red")) +
  labs(x = "Diet treatment") +
  theme(axis.text.y = element_text(size = 8))
# ```

# ```{r spec-KOs-heatmap-median, fig.height=8}
ggplot(heatmap.tbl, aes(y = KO, x = treat)) +
  geom_bin2d(aes(fill = Med_abund), stat = "identity") +
  scale_fill_gradientn(name = "Median scaled\nabundance", colors = c("blue", "white", "red")) +
  labs(x = "Diet treatment") +
  theme(axis.text.y = element_text(size = 8))
# ```

# for the above heatmaps, to account for large differences in total abundance, I scaled abundances by KO (centered on zero, divided by st. dev.), then plotted the mean/median of these scale values by treatment.

# ```{r specific-KOs-behavior, eval=FALSE}
ko.data <- ko.tbl[sample.dt]
cl <- makeCluster(maxCores, type = "FORK")
registerDoParallel(cl, maxCores)
ko.id <- focal.kos[1]
sig.pairs <- foreach(
  ko.id = focal.kos,
  .combine = "rbind",
  .verbose = TRUE
) %dopar% {
  focal.data <- copy(ko.data)
  ko.stat <- ko.id
  trans.yn <- "N"
  test <- "glm"
  if (shapiro.test(focal.data[[ko.id]])$p.value <= 0.05) {
    ko.trans <- rcompanion::transformTukey(ko.data[[ko.id]], plotit = F, quiet = T)
    if (shapiro.test(ko.trans)$p.value <= 0.05) {
      test <- "cpglm"
    } else {
      focal.data[, KO_trans := ko.trans]
      ko.stat = "KO_trans"
      trans.yn <- "Y"
    }
  }

  if (test == "glm") {
    frm <- paste(ko.stat, "~", paste(spat.learn.vars, collapse = " + "))
    full.model <- lm(as.formula(frm0), data = focal.data)
    step.model <- MASS::stepAIC(full.model, direction = "both", trace = F)
    model.dt <- as.data.table(tidy(step.model))

  } else if (test == "cpglm") {
    best.mods <- lapply(1:100, function(i) {
      min.AIC <- 1e10
      test.vars <- NULL
      keep.vars <- NULL
      for (var in sample(spat.learn.vars, length(spat.learn.vars))) {
        test.vars <- c(test.vars, var)
        cat(paste(test.vars, collapse = " "), sep = "\n")
        frm <- paste(ko.stat, "~", paste(test.vars, collapse = " + "))
        mod <- cplm::cpglm(as.formula(frm), data = focal.data)
        if (mod$aic < min.AIC) {
          keep.vars <- c(keep.vars, var)
          min.AIC <- mod$aic
          best.mod <- mod
        }
        test.vars <- keep.vars
      }
      return(best.mod)
    })
    aics <- sapply(best.mods, function(x) x$aic)
    best.mod <- best.mods[which(aics == min(aics))][[1]]
    rm(best.mods)
    model.dt <- as.data.table(cplm::summary(best.mod)$coefficients, keep.rownames = "term")
    names(model.dt)[5] <- "p.value"

  } else {
    stop("Something went wrong, the `test` variable should be 'glm' or 'cpglm'")
  }

  return(
    data.table(
      KO = ko.id,
      Abund_trans = trans.yn,
      Covar = model.dt[!grepl("Intercept", term) & p.value <= 0.05]$term
    )
  )
}
stopCluster(cl)
View(sig.pairs)
# ```


# ```{r specific-KOs-behavior1, eval=FALSE}
data1 <- ko.tbl[sample.dt]
ggplot(data1, aes(x = K10797)) +
  geom_histogram(binwidth = 0.00001)
frm1 <- paste(focal.kos[1], "~", paste(spat.learn.vars, collapse = " + "))
full.model <- glm(as.formula(frm1), data = data1)
step.model <- MASS::stepAIC(full.model, direction = "both", trace = F)
type1.anova <- summary(step.model)
type2.anova <- car::Anova(step.model, type = 2)
# print(type1.anova)
print(type2.anova)
step.model.dt <- as.data.table(tidy(step.model))
sig.terms <- step.model.dt[p.value <= 0.05, term]
sig.terms <- sig.terms[!grepl("Intercept", sig.terms)]
# data1.1 <- melt(
#   data1,
#   measure.vars = sig.terms,
#   variable.name = "Covar",
#   value.name = "Score"
# )
ggplot(data1, aes_string(x = sig.terms, y = "K10797")) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", size = 0.5)
# ```

# ```{r specific-KOs-behavior2, fig.height=8, eval=FALSE}
data2 <- ko.tbl[sample.dt]
ggplot(data2, aes(x = K01859)) +
  geom_histogram(binwidth = 0.0000001)
darwinCores <- 40
if (redo.permutations) {
  cl <- makeCluster(darwinCores, type = "FORK")
  registerDoParallel(cl, darwinCores)
  best.mods <- foreach(
    i = 1:1000,
    .verbose = TRUE
  ) %dopar% {
    min.AIC <- 1e10
    test.vars <- NULL
    keep.vars <- NULL
    for (var in sample(spat.learn.vars, length(spat.learn.vars))) {
      test.vars <- c(test.vars, var)
      cat(paste(test.vars, collapse = " "), sep = "\n")
      frm <- paste(focal.kos[2], "~", paste(test.vars, collapse = " + "))
      mod <- cplm::cpglm(as.formula(frm), data = data2)
      cat(mod$aic, sep = "\n")
      if (mod$aic < min.AIC) {
        keep.vars <- c(keep.vars, var)
        min.AIC <- mod$aic
        best.mod <- mod
      }
      test.vars <- keep.vars
    }
    return(best.mod)
  }
  stopCluster(cl)
  saveRDS(best.mods, file = file.path(saveDir, "list_K01859_cpglms.rds"))
} else {
  best.mods <- readRDS(file.path(saveDir, "list_K01859_cpglms.rds"))
}
aics <- sapply(best.mods, function(x) x$aic)
lowest.aic.mods <- best.mods[which(aics == min(aics))]
# sapply(lowest.aic.mods, function(x) x$coefficients)
best.mod <- lowest.aic.mods[[1]]
rm(best.mods)

print(best.mod)

# data2.1 <- melt(
#   data2,
#   measure.vars = names(best.mod$coefficients)[-1],
#   variable.name = "Covar",
#   value.name = "Score"
# )
ggplot(data2, aes(x = Vis1, y = log(K01859 + 1e-8))) +
  geom_point() +
  geom_abline(
    intercept = best.mod$coefficients[1],
    slope = best.mod$coefficients[2],
    color = "red"
  )
# ```

# ```{r, eval=FALSE}
library(tidymodels)

sample.data <- copy(sample.dt)
setkey(sample.data, Sample)
ko.data <- otu.data.table(ps.igc)
ko.mat <- otu.matrix(ps.igc)
setkey(ko.data, Sample)
kos.mod.link <- copy(kos.mod.dt)[, .(ko, mod)]
setkey(kos.mod.link, ko)

# focal vars: Hid3 & Hid456
focal.dt1 <- copy(sample.data)[, .(Sample, Hid3)][ko.data]
focal.dt1[, Sample := NULL]
focal.dt2 <- copy(sample.data)[, .(Sample, Hid456)][ko.data]
focal.dt2[, Sample := NULL]

focal.dt1.split <- initial_split(focal.dt1, prop = 0.6)
focal.dt1.split

# focal.dt1.split %>%
#   training() %>%
#   glimpse()

### RUN ON DARWIN
system('echo $HOME')
focal.dt1.recipe <- training(focal.dt1.split) %>%
  recipe(Hid3 ~ .) %>%
  step_corr(all_predictors(), threshold = 0.5) %>%
  prep()


focal.dt2.split <- initial_split(focal.dt2, prop = 0.6)
### RUN ON DARWIN
focal.dt2.recipe <- training(focal.dt2.split) %>%
  recipe(Hid456 ~ .) %>%
  step_corr(all_predictors(), threshold = 0.5) %>%
  prep()
# ```
