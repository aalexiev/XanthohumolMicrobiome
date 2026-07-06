# Iterative linear models: bile acid ~ Diet * behavior
# With FDR correction across all tests

library(tidyverse)

# Read data
setwd("~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/ba_diet_behavior/")
dat <- read.csv("BAdietbeh_meta.csv", header = TRUE)

# Identify bile acid columns (start with "var") and behavior columns (everything after Diet)
ba_cols <- grep("^var", names(dat), value = TRUE)
diet_col_index <- which(names(dat) == "Diet")
beh_cols <- names(dat)[(diet_col_index + 1):ncol(dat)]

# Run all combinations
results_list <- list()

for (ba in ba_cols) {
  for (beh in beh_cols) {

    formula <- as.formula(paste(ba, "~ Diet *", beh))
    mod <- lm(formula, data = dat)
    s <- summary(mod)

    coefs <- as.data.frame(s$coefficients)
    coefs$term <- rownames(coefs)
    coefs$bile_acid <- ba
    coefs$behavior <- beh
    rownames(coefs) <- NULL

    results_list[[length(results_list) + 1]] <- coefs
  }
}

# Combine all results
all_results <- bind_rows(results_list)
names(all_results)[1:4] <- c("Estimate", "Std.Error", "t.value", "p.value")

# FDR correction across all p-values
all_results$FDR <- p.adjust(all_results$p.value, method = "fdr")

# Filter significant results (FDR < 0.05)
sig_results <- all_results %>%
  filter(FDR < 0.05) %>%
  filter(term != "(Intercept)") %>%
  select(bile_acid, behavior, term, Estimate, Std.Error, FDR)

# Save significant results
# write.csv(sig_results, "significant_results.csv", row.names = FALSE)

# Print summary
cat("Total tests:", nrow(all_results), "\n")
cat("Significant (FDR < 0.05):", nrow(sig_results), "\n")


############## now make graphs ##############
library(ggpubr)

# Read legend and create var -> label lookup
legend <- read.csv("bileAcids_and_lipids_legend.csv")
ba_lookup <- setNames(legend$Label, legend$varID)

# Get unique bile acid + behavior combos from significant results
sig_combos <- sig_results %>%
  distinct(bile_acid, behavior)

# Build individual plots
plot_list <- list()

for (i in 1:nrow(sig_combos)) {
  ba  <- sig_combos$bile_acid[i]
  beh <- sig_combos$behavior[i]

  # Get the nice bile acid name (fall back to varID if not in legend)
  ba_label <- ifelse(ba %in% names(ba_lookup), ba_lookup[ba], ba)

  # Get stats for this combo (all terms)
  stats <- sig_results %>%
    filter(bile_acid == ba, behavior == beh)

  # Build annotation text
  ann_lines <- paste0(stats$term, ": Coeff=", round(stats$Estimate, 3),
                      ", p=", formatC(stats$FDR, format = "e", digits = 2))
  ann_text <- paste(ann_lines, collapse = "\n")

  # Make plot
  p <- ggplot(dat, aes(y = .data[[ba]], x = Diet)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_boxplot(aes(fill = Diet)) +
    # annotate("text", x = Inf, y = Inf, label = ann_text,
    #          hjust = 1.05, vjust = 1.2, size = 2.5, color = "grey30") +
    labs(y = ba_label, x = "Treatment") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")

  plot_list[[length(plot_list) + 1]] <- p
}

# Combine all plots into one figure
if (length(plot_list) > 0) {
  ncols <- 1
  combined <- ggarrange(plotlist = plot_list, ncol = ncols,
                        nrow = ceiling(length(plot_list) / ncols))

  ggsave("significant_BA_plots.tiff", combined,
         width = 5 * ncols,
         height = 4 * ceiling(length(plot_list) / ncols),
         dpi = 300)

  cat("Saved", length(plot_list), "plots to significant_BA_plots.tiff\n")
} else {
  cat("No significant results to plot.\n")
}



