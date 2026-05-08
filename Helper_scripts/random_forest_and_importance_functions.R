# random forest and importance functions

gen.diet.random.forests <- function(threshholds, physeq) {
  rf.list <- NULL
  for (prev.thresh in threshholds) {
    cat(paste("\nPrevalance threshold:", prev.thresh), sep = "\n")
    feature.mat <- otu.matrix(physeq)
    feature.mat[feature.mat > 0] <- 1
    feature.pa.sums <- colSums(feature.mat)
    keep.features <- names(
      feature.pa.sums[feature.pa.sums >= (prev.thresh * nrow(feature.mat))]
    )
    cat(paste0("\tNumber features passing threshold: ", length(keep.features)), sep = "\n")

    feature.dt <- as.data.table(
      otu.matrix(physeq)[, keep.features],
      keep.rownames = "Sample"
    ) %>% setkeyv("Sample")
    diet.dt <- sample.data.table(physeq)[, .(Sample, Diet)] %>% setkeyv("Sample")
    rf.data <- diet.dt[feature.dt]
    rf.data[, Sample := NULL]

    diet.recipe <- recipe(rf.data) %>%
      update_role(Diet, new_role = "outcome") %>%
      update_role(-Diet) %>%
      step_nzv(all_predictors())

    set.seed(42)
    rf.list[[paste0("Diet-prev", prev.thresh)]] <- caret::train(
      diet.recipe,
      data = as.data.frame(rf.data),
      method = "ranger",
      importance = "permutation",
      metric = "AUC",
      tuneGrid = expand.grid(
        mtry = c(100, 250, 500, 1000),
        splitrule = c("gini", "extratrees"),
        min.node.size = 1 # 1 for classification, 5 for regression according to `ranger` manual
      ),
      trControl = trainControl(
        method = "cv",
        number = 5,
        summaryFunction = multiClassSummary, # to calculate ROC
        classProbs = TRUE, # IMPORTANT for twoClassSummary
        savePredictions = TRUE,
        verboseIter = FALSE
      )
    )

    cat(paste0("DONE"), sep = "\n")
  }

  aucs <- sapply(rf.list, function(mod) {
    return(max(mod$results$AUC))
  }) %>% sort(decreasing = TRUE)

  return(rf.list[names(aucs)[1]])
}

gen.covar.random.forests <- function(
  covars,
  threshholds,
  physeq,
  feature.type,
  save.memory = FALSE
) {
  if (save.memory) {
    resDir <- "tmp"
    if (!dir.exists(resDir)) { dir.create(resDir) }
  } else {
    rf.list <- NULL
  }
  for (covar in covars) {
    if (save.memory) {
      covarDir <- file.path(resDir, covar)
      if (!dir.exists(covarDir)) { dir.create(covarDir) }
    } else {
      covar.rf.list <- NULL
    }
    for (prev.thresh in threshholds) {
      cat(paste("\nBehavior covariate:", covar), sep = "\n")
      cat(paste("\tPrevalance threshold:", prev.thresh), sep = "\n")
      if (save.memory) {
        covar.mod.res.file <- file.path(covarDir, paste0("prev", prev.thresh, ".rds"))
        if (file.exists(covar.mod.res.file)) {
          cat(paste("\t File", covar.mod.res.file, "exists, moving to next"), sep = "\n")
          next
        }
      }
      feature.mat <- otu.matrix(physeq)
      feature.mat[feature.mat > 0] <- 1
      feature.pa.sums <- colSums(feature.mat)
      keep.features <- names(
        feature.pa.sums[feature.pa.sums >= (prev.thresh * nrow(feature.mat))]
      )
      cat(paste0("\tNumber features passing threshold: ", length(keep.features)), sep = "\n")

      feature.dt <- as.data.table(
        otu.matrix(physeq)[, keep.features],
        keep.rownames = "Sample"
      ) %>% setkeyv("Sample")
      covar.cols <- c("Sample", covar, "Diet")
      covar.dt <- sample.data.table(physeq)[, ..covar.cols] %>% setkeyv("Sample")
      covar.dt[, Diet := as.integer(Diet)]
      rf.data <- covar.dt[feature.dt]
      rf.data[, Sample := NULL]

      covar.recipe <- recipe(rf.data) %>%
        update_role(matches(covar), new_role = "outcome") %>%
        update_role(-matches(covar)) %>%
        step_nzv(all_predictors())

      set.seed(42)
      covar.mod.res <- caret::train(
        covar.recipe,
        data = as.data.frame(rf.data),
        method = "ranger",
        importance = "permutation",
        # preProcess = c("knnImpute"),
        tuneGrid = expand.grid(
          mtry = c(100, 250, 500, 1000),
          splitrule = c("variance", "extratrees"),
          min.node.size = 5 # 1 for classification, 5 for regression according to `ranger` manual
        ),
        trControl = trainControl(
          method = "cv",
          number = 5,
          # summaryFunction = multiClassSummary, # to calculate ROC
          # classProbs = TRUE, # IMPORTANT for twoClassSummary
          savePredictions = TRUE,
          verboseIter = TRUE
        )
      )
      if (save.memory) {
        saveRDS(covar.mod.res, file = covar.mod.res.file)
      } else {
        covar.rf.list[[paste0("prev", prev.thresh)]] <- covar.mod.res
      }
      cat(paste0("DONE"), sep = "\n")
    }

    if (save.memory) {
      file.check <- length(
        list.files(path = resDir, pattern = paste0(covar, "-prev"), full.names = T)
      ) > 0
      if (file.check) {
        cat("\tFinal model already selected, moving to next covar", sep = "\n")
      } else {
        rmses <- sapply(list.files(path = covarDir, full.names = T), function(file) {
          mod <- readRDS(file)
          return(sqrt(mean((mod$pred$pred - mod$pred$obs)^2)))
        }) %>% sort(decreasing = FALSE)
        file.copy(
          from = names(rmses)[1],
          to = file.path(
            resDir,
            paste0(covar, "-", tail(str_split(names(rmses)[1], "/")[[1]], 1))
          )
        )
      }
    } else {
      rmses <- sapply(covar.rf.list, function(mod) {
        return(sqrt(mean((mod$pred$pred - mod$pred$obs)^2)))
      }) %>% sort(decreasing = FALSE)
      rf.list[[paste0(covar, "-", names(rmses)[1])]] <- covar.rf.list[[names(rmses)[1]]]
    }
  }
  if (!save.memory) {
    saveRDS(
      names(rf.list),
      file = paste0("Saved_objects/vec_selected_", feature.type, "_randForest_model_ids.rds")
    )
    return(rf.list)
  }
}

id.sig.important.features <- function(randForest.list = NULL, feature.type) {
  resDir <- "tmp"
  if (is.null(randForest.list)) {
    if (!dir.exists(resDir)) {
      stop("no random forest list was supplied, but no `tmp` directory exists.")
    }
    rds.files <- list.files(path = resDir, pattern = ".rds$", full.names = T)
    names(rds.files) <- str_remove_all(rds.files, "tmp/|.rds")
    covar.mod.names <- names(rds.files)
  } else {
    covar.mod.names <- names(randForest.list)
  }
  all.sig.impt.dt <- NULL
  for (covar.mod.name in covar.mod.names) {
    save.file <- file.path(resDir, paste0(covar.mod.name, "_sigImportFeats.rds"))
    if (file.exists(save.file) & is.null(randForest.list)) {
      all.sig.impt.dt <- readRDS(save.file)
      cat(paste("Output file for", covar.mod.name, "already exists, skipping..."), sep = "\n")
      next
    }
    covar <- str_split(covar.mod.name, "-")[[1]][1]
    prev.thresh <-  str_remove(str_split(covar.mod.name, "-")[[1]][2], "prev") %>%
      as.numeric()
    cat(paste("\n###", covar.mod.name), sep = "\n")
    if (is.null(randForest.list)) {
      covar.rf <- readRDS(rds.files[covar.mod.name])
    } else {
      covar.rf <- randForest.list[[covar.mod.name]]
    }
    x <- covar.rf$finalModel
    dependent.var <- covar
    data <- covar.rf$trainingData
    num.permutations <- 100

    dat_x <- data[, -c(which(names(data) == dependent.var))]
    dat_x <- dat_x[, names(x$variable.importance)] %>% as.matrix()
    set.seed(42)
    vimp <- foreach(
      i = 1:num.permutations,
      .final = rlist::list.cbind,
      .verbose = TRUE
    ) %dopar% {
      dat_y <- data[sample(nrow(data)), dependent.var]
      if (class(dat_y) == "character") {dat_y <- factor(dat_y)}
      ranger::ranger(
        x = dat_x,
        y = dat_y,
        num.trees = x$num.trees,
        mtry = x$mtry,
        min.node.size = x$min.node.size,
        importance = x$importance.mode,
        replace = x$replace,
        num.threads = 1
      )$variable.importance
    }

    pval <- sapply(1:nrow(vimp), function(i) {
      (sum(vimp[i, ] >= x$variable.importance[i]) + 1)/(ncol(vimp) + 1)
    })

    covar.impt.mat <- cbind(x$variable.importance, pval)
    colnames(covar.impt.mat) <- c("importance", "pvalue")

    covar.sig.impt.dt <- covar.impt.mat[covar.impt.mat[, "pvalue"] <= 0.05, ] %>%
      as.data.table(keep.rownames = feature.type)
    covar.sig.impt.dt[, Covar := covar]
    all.sig.impt.dt <- rbind(
      all.sig.impt.dt,
      covar.sig.impt.dt[order(importance, decreasing = T)]
    )
    if (is.null(randForest.list)) { saveRDS(all.sig.impt.dt, file = save.file) }
  }
  return(all.sig.impt.dt)
}

trans.and.cat.covars <- function(dt, covars, reference = FALSE) {
  ref.dt <- NULL
  for (covar in covars) {
    covar.ref.dt <- data.table(Covar = covar, Transformed = FALSE, Categorized = FALSE)
    if (any(dt[[covar]] < 0)) {
      lambda.sets <- set_names(c("neg", "pos"), c("neg", "pos"))
      transforms <- lapply(lambda.sets, function(lambda.set) {
        vals <- transformTukey(
          dt[[covar]],
          start = ifelse(lambda.set == "neg", -20, 1),
          end = ifelse(lambda.set == "neg", -1, 20),
          int = 1,
          quiet = TRUE,
          plotit = FALSE
        )
        lambda <- transformTukey(
          dt[[covar]],
          start = ifelse(lambda.set == "neg", -20, 1),
          end = ifelse(lambda.set == "neg", -1, 20),
          int = 1,
          quiet = TRUE,
          plotit = FALSE,
          returnLambda = TRUE
        )
        return(list(Vals = vals, Lambda = lambda))
      })
      choice <- sapply(lambda.sets, function(lambda.set) {
        shapiro.test(transforms[[lambda.set]]$Vals)$statistic %>% return()
      }) %>%
        which.max()
      if (transforms[[choice]]$Lambda != 1) {
        covar.ref.dt$Transformed <- TRUE
        dt[, (covar) := trans]
      }
    } else {
      if (shapiro.test(dt[[covar]])$p.value <= 0.05) {
        trans <- rcompanion::transformTukey(dt[[covar]], plotit = F, quiet = T)
        covar.ref.dt$Transformed <- TRUE
        if (shapiro.test(trans)$p.value > 0.05) {
          dt[, (covar) := trans]
        } else {
          if (covar == "P_Cross") {
            dt[, (covar) := factor(dt[[covar]])]
            covar.ref.dt$Transformed <- FALSE
          } else {
            dt[
              ,
              (covar) := factor(
                ifelse(dt[[covar]] == max(dt[[covar]]), "High", "Low"),
                levels = c("Low", "High")
              )
            ]
            covar.ref.dt$Categorized <- TRUE
          }
        }
      }
    }
    ref.dt <- rbind(ref.dt, covar.ref.dt)
  }
  setkey(ref.dt, Covar)
  if (!reference) {
    return(dt)
  } else {
    return(list(DT = dt, REF.DT = ref.dt))
  }
}


diet.impt.feature.kw.tests <- function(sig.feature.dt, physeq, feature.type) {
  setkeyv(sig.feature.dt, feature.type)
  feature.mat <- otu.matrix(physeq)
  mod.dt <- feature.mat[, sig.feature.dt[[feature.type]]] %>%
    as.data.table(keep.rownames = "Sample")
  mod.dt[, Diet := sample.data.table(physeq)$Diet]
  kw.pvals <- sapply(
    setNames(sig.feature.dt[[feature.type]], sig.feature.dt[[feature.type]]),
    function(feature) {
      cat(feature, sep = "\n")
      return(
        kruskal.test(x = mod.dt[[feature]], g = mod.dt$Diet)$p.value
      )
    })

  kw.adj.pvals <- p.adjust(kw.pvals, method = "bonferroni")
  sig.features <- sig.feature.dt[names(kw.adj.pvals[kw.adj.pvals <= 0.05])]
  sig.features[, KW.adj.pval := kw.adj.pvals[kw.adj.pvals <= 0.05]]
  sig.features <- sig.features[order(-importance, KW.adj.pval)]
  return(sig.features)
}

covar.impt.feature.regressions <- function(
  covars,
  sig.feature.dt,
  covar.dt,
  physeq,
  feature.type
) {
  dt <- NULL
  for (covar in covars) {
    covar.sig.impt.dt <- sig.feature.dt[Covar == covar]
    ordered.features <- covar.sig.impt.dt[[feature.type]]
    feature.mat <- otu.matrix(physeq)
    feature.dt <- feature.mat[, ordered.features[ordered.features != "Diet"]] %>%
      as.data.table(keep.rownames = "Sample") %>%
      setkeyv("Sample")
    mod.dt.cols <- c("Sample", "Diet", covar)
    mod.dt <- covar.dt[, ..mod.dt.cols][feature.dt]

    glm.fam <- "gaussian"
    if (is.factor(mod.dt[[covar]])) {glm.fam <- "binomial"}
    model.type <- ifelse(covar == "P_Cross", "polr", "glm")

    intrxn.pval <- 1
    idx <- 0
    while (intrxn.pval > 0.05 & idx < length(ordered.features)) {
      idx <- idx + 1
      cat(paste0("covar <- '", covar, "'; idx <- ", idx), sep = "\n")
      frm <- paste(covar, "~ Diet *", ordered.features[idx])
      # set.seed(42)
      if (model.type == "glm") {
        glm0 <- try(
          glm(
            as.formula(frm),
            family = glm.fam,
            data = mod.dt
          ),
          silent = TRUE
        )
      } else {
        glm0 <- try(
          MASS::polr(
            as.formula(frm),
            data = mod.dt,
            Hess = TRUE
          ),
          silent = TRUE
        )
      }
      print(glm0)
      if ("try-error" %in% class(glm0)) {
        next
      } else {
        # set.seed(42)
        pval <- car::Anova(glm0, type = 2) %>% tidy() %>% use_series(p.value) %>% tail(1)
        intrxn.pval <- ifelse(is.na(pval), 1, pval)
      }
    }
    if (intrxn.pval <= 0.05 | idx < length(ordered.features)) {
      kept.features <- ordered.features[idx]
      curr.glm <- glm0
      for (i in (idx + 1):length(ordered.features)) {
        cat(paste0("covar <- '", covar, "'; i <- ", i), sep = "\n")
        test.features <- c(kept.features, ordered.features[i])
        new.glm <- try(
          update(curr.glm, as.formula(paste(". ~ . + Diet *", ordered.features[i]))),
          silent = TRUE
        )
        if ("try-error" %in% class(new.glm)) {
          cat("try-error on model update", sep = "\n")
          next
        } else {
          pval <- try(
            car::Anova(new.glm, type = 2) %>% tidy() %>% use_series(p.value) %>% tail(1),
            silent = TRUE
          )
          if (class(pval) == "try-error") {
            cat("try-error on car::Anova")
            next
          } else {
            print(car::Anova(new.glm, type = 2))
            intrxn.pval <- ifelse(is.na(pval), 1, pval)
            if (!is.na(intrxn.pval)) {
              if (intrxn.pval <= 0.05) {
                curr.glm <- new.glm
                kept.features <- test.features
              }
            }
          }
        }
      }
      final.model <- Anova(curr.glm, 2) %>% tidy() %>% as.data.table()
      results <- data.table(
        Covar = covar,
        V2 = kept.features,
        Intrxn.pval = final.model[term %in% str_subset(term, "Diet:")]$p.value
      )
      names(results)[2] <- feature.type
      dt <- rbind(dt, results)
    }
  }
  return(dt)
}
