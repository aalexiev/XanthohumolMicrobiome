altmann.importance.pvalues <- function (
  x,
  num.permutations = 100,
  dependent.var = NULL,
  data = NULL,
  nCores = 4,
  user.seed = 42,
  external.cluster = FALSE
) {
  require(doParallel)
  require(foreach)
  require(ranger)
  require(rlist)

  dat_x <- data[, -c(which(names(data) == dependent.var))]
  dat_x <- dat_x[, names(x$variable.importance)] %>% as.matrix()

  set.seed(user.seed)
  if (!external.cluster) {
    cl <- makeCluster(nCores, type = "FORK", outfile="")
    registerDoParallel(cl, nCores)
  }
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
  if (!external.cluster) {
    stopCluster(cl)
  }

  pval <- sapply(1:nrow(vimp), function(i) {
    (sum(vimp[i, ] >= x$variable.importance[i]) + 1)/(ncol(vimp) + 1)
  })

  res <- cbind(x$variable.importance, pval)
  colnames(res) <- c("importance", "pvalue")
  return(res)
}
