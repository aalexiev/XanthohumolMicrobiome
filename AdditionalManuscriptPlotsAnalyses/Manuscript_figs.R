## This script is to make a heatmap of KO's that are predictive of behavior metrics
# Keaton previously found these KO's using random forest models
# then he put those KO abudances and behavior metrics into linear regressions to get significance

# read in his output results file
# setwd("/Users/alexieva/Documents/Projects/Analysis/xanthahumol-master")
library(dplyr)
# res <- readRDS("Saved_objects/dt_allCovar_randForests_imptKOs_sig_interactions.rds") %>%
#   inner_join(tab4[1:2], by = "KO") %>%
#   distinct(.keep_all = TRUE)

# write.table(res, "~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/rf_regress_results_summary.txt",
#             quote = F, sep = "\t")

# and then I lost tab4 wherever that thing is...
res <- read.delim("~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/rf_regress_results_summary_original.txt",
                  sep = "\t")


## testing some functions in Keaton's code using dummy data so I can add something to a function
# library(car)
# test <- data.frame(r = rnorm(50), d = rnorm(50), f = rnorm(50),
#                    g = as.factor(rep(c("A", "B", "C", "D", "E"), 10)), h = sample(0:1, 50, replace = T))
# fit <- glm(r~d*f, family = "gaussian", data = test)
# car::Anova(fit, type = "II") %>% tidy()
# sum_fit <- summary(fit)
# sum_fit_clean <- sum_fit$coefficients %>% as.data.frame() %>% rownames_to_column("term") %>%
#   dplyr::select(Estimate = 2, everything()) %>% as.data.table()
# sum_fit_clean[term %in% str_subset(term, "d:")]$Estimate
# fit <- MASS::polr(g~d*f, data = test, Hess = TRUE)
# summary(fit)
# fit <- glm(h~d*f, family = "binomial", data = test)

########################## make a heatmap ##########################
## first of p-values
res <- res %>%
  mutate(signif = case_when(Intrxn.pval <= 0.05 ~ "*",
                                      TRUE ~ "ns")) %>%
  dplyr::filter(signif != "ns")

heat_pvals <- ggplot(data = res, aes(x = Covar, y = KO.description)) +
  geom_tile(aes(fill = Intrxn.pval), colour = "black", linetype = 1) + # here I use IRR but more often people use lm coefficients or significance
  theme_classic() +
  scale_fill_gradient2(mid = "#d0c098", high = "#1ab8dc",
                       na.value = "grey50",
                       guide = "colourbar") +
  labs(y = "", x = "Model term, interacting with diet") +
  # geom_text(aes(label = signif)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_blank())

# ggsave("Plots/heatmap_plot_pvals.png",
#        width = 16,
#        height = 8,
#        dpi = 300)

## then of the chi squared terms
heat_pvals <- ggplot(data = res, aes(x = Covar, y = KO.description)) +
  geom_tile(aes(fill = Chisq_var), colour = "black", linetype = 1) + # here I use IRR but more often people use lm coefficients or significance
  theme_classic() +
  scale_fill_gradient2(mid = "#d0c098", high = "#1ab8dc",
                       na.value = "grey50",
                       guide = "colourbar") +
  labs(y = "", x = "Model term, interacting with diet") +
  geom_text(aes(label = signif)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_blank())

# ggsave("Plots/heatmap_plot_chisq.png",
#        width = 16,
#        height = 8,
#        dpi = 300)


########################## Heatmap Tom edits ##########################
# get data frame of KO's and XN's and behaviors


# run kendall's correlation (non-parametric, good for small samples and some outliers)
# between each KO and each XN metabolite and behavior combo (write a loop for this)


# make heatmap with legend in axis
## then of the chi squared terms
heat_pvals <- ggplot(data = res, aes(x = Covar, y = KO.description)) +
  geom_tile(aes(fill = Chisq_var), colour = "black", linetype = 1) + # here I use IRR but more often people use lm coefficients or significance
  theme_classic() +
  scale_fill_gradient2(mid = "#d0c098", high = "#1ab8dc",
                       na.value = "grey50",
                       guide = "colourbar") +
  labs(y = "", x = "Model term, interacting with diet") +
  geom_text(aes(label = signif)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_blank())

# ggsave("Plots/heatmap_plot_chisq.png",
#        width = 16,
#        height = 8,
#        dpi = 300)


######################## some file clean up and dunn's test ##########################

## now figure out running a dunn's test on the kruskall wallis Keaton ran for KO's associated with diet
# kw.sig.kos.dt <- readRDS("Saved_objects/dt_diet_randForests_imptKOs_sig_effects.rds")
# metadata <- read.csv("Input/SpatialLearning_metadata copy.csv")
# kos <- readRDS("Input/phyloseq.rds")
# topkos <- list("K15780", "K09961", "K06079", "K19353", "K03503", "K06297", "K01729",
#                "K18197", "K00841", "K08169", "K17239", "K08723", "K00179", "K00180",
#                "K01023", "K01667", "K06445", "K07516", "K10797", "K18244")
#
# library(phyloseq)
library(phyloseqCompanion)
meta_from_phylo <- sample.data.frame(phylo) %>%
  rownames_to_column("ID")
# kos_taxtab <- otu.data.table(kos) %>%
#   dplyr::select(as.character(append(topkos, "Sample")))
meta_from_phylo$ID <- gsub("sa", "s", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s1$", "s01", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s2$", "s02", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s3$", "s03", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s4$", "s04", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s5$", "s05", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s6$", "s06", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s7$", "s07", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s8$", "s08", meta_from_phylo$ID)
meta_from_phylo$ID <- gsub("s9$", "s09", meta_from_phylo$ID)
# kos_taxtab2 <- kos_taxtab %>%
#   inner_join(metadata, by = "Sample") %>%
#   column_to_rownames("Sample")
# write.csv(meta_from_phylo, "Metadata_from_phyloseq.csv")

# actually just use plot.dt object from README.Rmd chunk #60
plot.dt
plot.dt_long <- plot.dt %>%
  pivot_wider(id_cols = c(Sample, Diet), names_from = KO, values_from = Abund) %>%
  column_to_rownames("Sample")

# run dunn test
# krusk <- apply(kos_taxtab2, 2, function(x) kruskal.test(Diet ~ topkos$x, data = kos_taxtab2))
kruskal.test(Diet ~ K15780, data = plot.dt_long)
ggplot(data = plot.dt_long, aes(x = Diet, y = K15780)) +
  geom_point()

install.packages('dunn.test')
library(dunn.test)
dunn.test(plot.dt_long$K15780, plot.dt_long$Diet, method = "bonferroni")
dunn_res <- c()
for (i in 1:13) {
  name <- names(plot.dt_long)[1+i]
  dunn_res[[name]] <- dunn.test(plot.dt_long[,1+i], plot.dt_long$Diet, method = "bonferroni")
}

# now the ko's of interest from fig 1D
interestkos <- list("K00179", "K00180",
               "K01023", "K01667", "K06445",
               "K07516", "K10797", "K18244")
xn.ko.plot.dt_long <- xn.ko.plot.dt %>%
  rownames_to_column("Sample") %>%
  pivot_wider(id_cols = c(Sample, Diet), names_from = KO, values_from = Abund) %>%
  column_to_rownames("Sample")

dunn_res2 <- c()
for (i in 1:13) {
  name <- names(xn.ko.plot.dt_long)[1+i]
  dunn_res2[[name]] <- dunn.test(xn.ko.plot.dt_long[,1+i], xn.ko.plot.dt_long$Diet, method = "bonferroni")
}

# check KO names to add to manscript text
ko_names <- readRDS("Saved_objects/my_ko_names_and_mod_name.rds")


######################## get random forest info ##########################
library(RedoControl)
# diet only rfs - these are classification models
library(randomForest)

selected.diet.rf # this is the object from keaton's readme that holds the diet rfs, has info on model parameters
selected.diet.rf$`Diet-prev0.4`$finalModel # this gave me OOB error

# rfs of all covariates and model fit results
rf.list

######################## bile acids checks ##########################
bas <- read.csv("~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/bile_acids_micr.csv")
KO_tab <- read.table("Input/Merged_Relative_abundances.tab") %>%
  dplyr::select(bas$KO)
library("phyloseq")
library("tidyverse")
phylo <- readRDS("~/Documents/Projects/Analysis/XN_keaton/xanthahumol-summary/Input/phyloseq.rds")
# Extract abundance matrix from the phyloseq object
OTU1 <- as(otu_table(phylo), "matrix")
# transpose if necessary
if(taxa_are_rows(phylo)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf <- as.data.frame(OTU1) %>%
  dplyr::select(bas$KO) %>%
  rownames_to_column("ID")
metadata <- read.csv("~/Documents/Projects/Analysis/XN_keaton/xanthahumol-summary/Input/SpatialLearning_metadata.csv")
OTUdf$ID <- gsub("sa", "s", OTUdf$ID)
OTUdf$ID <- gsub("s1$", "s01", OTUdf$ID)
OTUdf$ID <- gsub("s2$", "s02", OTUdf$ID)
OTUdf$ID <- gsub("s3$", "s03", OTUdf$ID)
OTUdf$ID <- gsub("s4$", "s04", OTUdf$ID)
OTUdf$ID <- gsub("s5$", "s05", OTUdf$ID)
OTUdf$ID <- gsub("s6$", "s06", OTUdf$ID)
OTUdf$ID <- gsub("s7$", "s07", OTUdf$ID)
OTUdf$ID <- gsub("s8$", "s08", OTUdf$ID)
OTUdf$ID <- gsub("s9$", "s09", OTUdf$ID)
OTUdf2 <- OTUdf %>%
  inner_join(metadata, by = "ID") %>%
  column_to_rownames("ID")
OTUdf2 <- OTUdf2[-47,]

# graph each KO over diet and cognition
source("Helper_scripts/ko_diet_plots.R")
plots <- plot_kos_by_diet(OTUdf2, diet_col = "Diet")
plots[[1]]
plots[[2]]

# try to link the BA abundance from old data to gut microbiome metrics
library(tidyr)
library(vegan)
# mice <- read.csv("Input/individuals_treatment.csv")
BA_metadata <- read.csv("Input/bileAcids_and_lipids_legend.csv") %>%
  dplyr::filter(Data_type == "normalized-logTrans-pareto",
                Group == "Secondary",
                Category == "bile acids feces" | Category == "bile acids liver") %>%
  dplyr::select(-"Data_type")

BA_tab <- read.csv("Input/reformatted_bileAcids_and_lipids_data.csv") %>%
  dplyr::select(c("Sample", BA_metadata$varID))
BA_tab <- BA_tab[1:42,]

OTUdf3 <- OTUdf2 %>%
  rownames_to_column("Sample") %>%
  dplyr::filter(Sample %in% BA_tab$Sample) %>%
  column_to_rownames("Sample")

bc_div <- vegdist(OTUdf3[,1:8], method = "bray")

BA_tab_w_metadata <- BA_tab %>%
  rename(ID = Sample) %>%
  inner_join(metadata, by = "ID") %>%
  inner_join(meta_from_phylo[,c(1,46:47)], by = "ID")


## stats
# microbiome beta diversity ~ metabolite abund
# BAtabtest <- BA_tab_w_metadata %>% dplyr::filter(ID == "s01")
names <- names(BA_tab_w_metadata[,2:12])
res_perms <- c()
for (col in names) {
  formula <- as.formula(paste("bc_div ~", names, " + ", names, ":Diet"))
  perm <- adonis2(formula, data = BA_tab_w_metadata)
  if (perm$`Pr(>F)`[1] <= 0.05) {
    res_perms[[col]] <- perm
  } else {
    res_perms[[col]] <- "Not Significant"
  }
}
# none significant; neither was just BA in formula

# microbiome alpha diveristy ~ metabolite abund + abund:diet
names <- names(BA_tab_w_metadata[,2:12])
res_perms <- c()
for (col in names) {
  formula <- as.formula(paste("shannon ~", names, " + ", names, ":Diet"))
  perm <- summary(lm(formula, data = BA_tab_w_metadata))
  if (any(perm$coefficients[2:5,4] <= 0.05)) {
    res_perms[[col]] <- perm
  } else {
    res_perms[[col]] <- "Not Significant"
  }
}
# none significant for richness or shannon

## exploratory graphs
BA_tab_w_meta2 <- BA_tab_w_metadata %>%
  column_to_rownames("ID") %>%
  pivot_longer(cols = 2:12,
               names_to = "varID",
               values_to = "abundance")
# example graph
ggplot(data = BA_tab_w_metadata, aes(x = richness, y = var001, color = Diet)) +
  geom_point() +
  geom_smooth(method = 'lm')

## redo with BH correction 
# do these bile acids differ with treatment?
ba_names <- names(BA_tab_w_metadata[, 2:12])   # 11 secondary bile acids

# fit lm(BA ~ Diet) per bile acid; capture the omnibus treatment p-value (F-test)
ba_lm_res <- lapply(setNames(ba_names, ba_names), function(col) {
  frm <- as.formula(paste(col, "~ Diet"))       # was: paste(names, ...) -> bug
  fit <- lm(frm, data = BA_tab_w_metadata)
  sm  <- summary(fit)
  fstat <- sm$fstatistic
  omnibus_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  list(fit = fit, coefficients = sm$coefficients, omnibus_p = unname(omnibus_p))
})

# BH-correct the omnibus treatment p-values across all 11 bile acids
raw_p <- sapply(ba_lm_res, `[[`, "omnibus_p")
bh_p  <- p.adjust(raw_p, method = "BH")

ba_p_table <- data.frame(
  bile_acid = ba_names,
  raw_p     = raw_p,
  BH_p      = bh_p,
  row.names = NULL
)
ba_p_table <- ba_p_table[order(ba_p_table$BH_p), ]
print(ba_p_table)

# graphs
BA_legend <- read.csv("Input/bileAcids_and_lipids_legend.csv") %>%
  dplyr::filter(Group == "Secondary") %>%
  dplyr::select(c("Label", "varID")) %>%
  mutate(combo = paste0(varID, "_", Label))

# check order is same
BA_legend$varID == names(BA_tab_w_metadata[,2:12])

BA_tab_w_meta2 <- BA_tab_w_metadata %>%
  rename_with(~ BA_legend$combo, all_of(names(BA_tab_w_metadata[,2:12])))

# --- Characterize which contrast drives each significant bile acid ---

# 1. Which bile acids are significant after BH?
sig_bas <- ba_p_table$bile_acid[ba_p_table$BH_p < 0.05]

# 2. Pull the Diet contrast rows from each significant bile acid's coefficient table
contrast_summary <- do.call(rbind, lapply(sig_bas, function(col) {
  coefs <- ba_lm_res[[col]]$coefficients
  diet_rows <- grep("^Diet", rownames(coefs))          # XN/DXN/TXN vs control
  data.frame(
    bile_acid = col,
    label     = BA_legend$Label[match(col, BA_legend$varID)],
    contrast  = sub("^Diet", "", rownames(coefs)[diet_rows]),
    estimate  = coefs[diet_rows, 1],
    std_error = coefs[diet_rows, 2],
    p_value   = coefs[diet_rows, 4],
    row.names = NULL
  )
}))

# 3. Order for readability: by bile acid, then by contrast p-value
contrast_summary <- contrast_summary[
  order(contrast_summary$bile_acid, contrast_summary$p_value), ]

print(contrast_summary, digits = 3)

##################### BA KO rf regressions #####################
KOdf_all <- as.data.frame(OTU1) %>%
  rownames_to_column("Sample")
KOdf_all$Sample <- gsub("sa", "s", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s1$", "s01", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s2$", "s02", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s3$", "s03", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s4$", "s04", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s5$", "s05", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s6$", "s06", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s7$", "s07", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s8$", "s08", KOdf_all$Sample)
KOdf_all$Sample <- gsub("s9$", "s09", KOdf_all$Sample)
KOdf_all2 <- KOdf_all %>%
  inner_join(BA_tab, by = "Sample")
KOdf_all2 <- KOdf_all2[-47,]

# write.csv(KOdf_all2, "~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/BA_KO_dfforrf.csv")

## create random forest model to test which taxa are associated with cumulative SFN_NIT
# default uses 500 trees
# regression model

# source this function that will run rf's and create/save importance graphs of
# top 20 KO's,
source("rf_ba_ko/rf_bile_acid_function.R")
run_rf_bile_acids(train_df = TrainSet,
                  test_df = ValidSet,
                  seed = 3)
# Saved: /Users/alexieva/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/model_fit_summary.csv
# Saved: /Users/alexieva/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/significant_genes_per_bile_acid.csv
# Saved: /Users/alexieva/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/full_importance_all_genes.csv


# run lm's looping through combos of bile acids and important KO's from the rf
library(broom)

# these are the significant BAs and KOs
rfsignifKO_BA <- read.csv("rf_ba_ko/significant_genes_per_bile_acid.csv")

# Get top 20 KO-bile acid pairs per bile acid from RF results
top20_pairs <- rfsignifKO_BA %>%
  group_by(BileAcid) %>%
  slice_max(order_by = Importance_pctIncMSE, n = 20) %>%
  ungroup() %>%
  select(BileAcid, Gene) %>%
  rename(bile_acids = BileAcid, KO = Gene)

# filter the original KO BA data frame for these
rfKO_BA_tab <- KOdf_all2 %>%
  dplyr::select(c(unique(rfsignifKO_BA$BileAcid), unique(rfsignifKO_BA$Gene)))

# do a quick check of all the bile acid distributions to make sure they're gaussian
ggplot(data = reshape2::melt(rfKO_BA_tab[,1:11]),
       aes(x = value)) +
  geom_density(fill = "skyblue") +
  stat_function(
    fun  = dnorm,
    args = list(mean = mean(reshape2::melt(rfKO_BA_tab[,1:11])$value, na.rm = TRUE),
                sd   = sd(reshape2::melt(rfKO_BA_tab[,1:11])$value,   na.rm = TRUE)),
    color    = "red",
    linewidth = 1
  ) +
  facet_wrap("variable")

# lets statistically check for normality since var001 and var003 are not so normal looking
reshape2::melt(rfKO_BA_tab[,1:11]) %>%
  group_by(variable) %>%
  summarise(W = shapiro.test(value)$statistic,
            p = shapiro.test(value)$p.value)
# p non-significant means it is normal
# var001, var003, var016, and var022 are non-normal distributions

# need to either transform these into normality or find a better fitting distribution
reshape2::melt(rfKO_BA_tab[,1:11]) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(variable) %>%
  summarise(W_orig = shapiro.test(value)$statistic,
            p_orig = shapiro.test(value)$p.value,
            W_log  = shapiro.test(log1p(value))$statistic,
            p_log  = shapiro.test(log1p(value))$p.value) %>%
  mutate(normal = ifelse(p_orig >= 0.05, "original", ifelse(p_log >= 0.05, "log1p", "non-normal")))


reshape2::melt(rfKO_BA_tab[,1:11]) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(variable) %>%
  summarise(p_orig = shapiro.test(value)$p.value,
            p_log  = shapiro.test(log1p(value))$p.value,
            p_sqrt = shapiro.test(sqrt(abs(value)))$p.value,
            p_sq   = shapiro.test(value^2)$p.value,
            p_ref  = shapiro.test(log1p(max(value) - value))$p.value) %>%
  mutate(across(starts_with("p_"), ~ round(.x, 4)))

# transform the four "non-normal" distributed BA's with a square root transformation
rfKO_only_tab <- rfKO_BA_tab[,12:221] %>%
  rownames_to_column("ID") %>%
  reshape2::melt()

trans_rfKO_BA_tab <- rfKO_BA_tab[,1:11] %>%
  mutate(var001 = sqrt(abs(var001)),
         var003 = sqrt(abs(var003)),
         var016 = sqrt(abs(var016)),
         var022 = sqrt(abs(var022))) %>%
  rownames_to_column("ID") %>%
  reshape2::melt() %>%
  inner_join(rfKO_only_tab, by = "ID") %>%
  inner_join(rownames_to_column(as.data.frame(metadata[,2]), "ID"), by = "ID") %>%
  rename(bile_acids = variable.x,
         ba_levels = value.x,
         KO = variable.y,
         abundance = value.y,
         Diet = 6) %>%
  semi_join(top20_pairs, by = c("KO", "bile_acids"))

## run lm model with below code with bile acids and KOs that were important

results_tidy <- trans_rfKO_BA_tab %>%
  group_by(KO, bile_acids) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(abundance ~ ba_levels + Diet + ba_levels:Diet, data = .)),
    tidied = map(model, tidy),
    glanced = map(model, glance)
  ) %>%
  unnest(tidied, glanced, .drop = TRUE)

print(results_tidy)

# clean up results to only focus on significance at bile_acids or bile acids X diet level
ba_focused <- results_tidy %>%
  filter((term == "ba_levels" | str_starts(term, "ba_levels:Diet")) & p.value < 0.05) %>%
  select(KO, bile_acids, term, estimate, std.error, p.value) %>%
  mutate(
    term_type = ifelse(term == "ba_levels", "main_effect", "interaction"),
    interaction_diet = ifelse(term == "ba_levels", NA,
                              str_extract(term, "Diet[^\\s,]+"))
  ) %>%
  arrange(KO, bile_acids, p.value)

print("\nBile acid-focused results (main effects + interactions, excluding Diet-only terms):")
print(ba_focused)
# write.csv(ba_focused, "rf_ba_ko/bileacid_focused_significant_results.csv", row.names = FALSE)

ba_focused %>% count(bile_acids)
nrow(unique(ba_focused[, c("KO", "bile_acids")]))
# 36 unique combos of significant KOs to bile acids
KO_bas <- data.frame("Feature" = unique((as.data.frame(ba_focused %>% count(bile_acids)))[,1]))
ko_names_bafocused <- KO_bas %>%
  inner_join(ko_names[,1:2], by = "Feature")


## now make graphs of all these
# ============================================================
# KO Abundance ~ Bile Acid Level Plots
# Plot 1: Main effects only (term == "ba_levels")
# Plot 2: Diet interaction plots (term == "ba_levels:DietXXX")
# ============================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)

# --- INPUTS ---
# ba_focused   : the significant results data frame (loaded below)
# trans_rfKO_BA_tab : long-format data frame with columns:
#                     Sample, bile_acids, ba_levels, KO, abundance, Diet

# ============================================================
# SHARED AESTHETICS
# ============================================================
shapes_list    <- c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0, 6, 7, 9, 10, 11, 12, 13, 14)
linetypes_list <- c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash")
# diet_colors    <- c("CTR"  = "#999999",
#                     "XN"   = "#E69F00",
#                     "DXN"  = "#56B4E9",
#                     "TXN"  = "#009E73")

# ============================================================
# PLOT 1: MAIN EFFECTS (term == "ba_levels", no diet interaction)
# One facet per bile acid; different KOs shown with shape + linetype
# Estimate printed in corner
# ============================================================

main_pairs <- ba_focused %>%
  filter(term_type == "main_effect") %>%
  select(KO, bile_acids, estimate, p.value) %>%
  distinct()

# Join raw data to significant main-effect KO-bile acid pairs
main_data <- trans_rfKO_BA_tab %>%
  inner_join(main_pairs, by = c("KO", "bile_acids"))

# Build estimate label per bile acid (all KOs stacked)
main_labels <- main_pairs %>%
  group_by(bile_acids) %>%
  summarise(
    label = paste(paste0(KO, ": β=", formatC(estimate, format = "e", digits = 2)),
                  collapse = "\n"),
    .groups = "drop"
  )

# Assign shape and linetype per KO (consistent across facets)
ko_levels_main <- unique(main_data$KO)
ko_aes_main <- data.frame(
  KO       = ko_levels_main,
  ko_shape = shapes_list[seq_along(ko_levels_main) %% length(shapes_list) + 1],
  ko_lty   = linetypes_list[seq_along(ko_levels_main) %% length(linetypes_list) + 1]
)
main_data <- main_data %>% left_join(ko_aes_main, by = "KO")

# Build per-panel estimate labels (one per KO-bile acid combo)
main_labels_grid <- main_pairs %>%
  mutate(label = paste0("coeff=", formatC(estimate, format = "e", digits = 2)))

plot1 <- ggplot(main_data, aes(x = ba_levels, y = abundance)) +
  geom_point(alpha = 0.5, size = 1.2, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  geom_text(data    = main_labels_grid,
            mapping = aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.4,
            size = 5, inherit.aes = FALSE) +
  facet_wrap(KO ~ bile_acids, scales = "free") +
  labs(
    x     = "Bile Acid Level",
    y     = "KO Abundance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text.x     = element_text(face = "bold", size = 8),
    strip.text.y     = element_text(face = "bold", size = 7, angle = 0),
    panel.spacing    = unit(0.3, "lines"),
    plot.title       = element_text(face = "bold"),
    legend.position  = "none",
    text      = element_text(size = 16),
  )

# ggsave("rf_ba_ko/BA_KO_main_effects.tiff", plot1,
#        width  = 14,
#        height = 10,
#        dpi = 300)

# ============================================================
# PLOT 2: DIET INTERACTIONS
# One facet per bile acid; colored lines per Diet
# Asterisk appended to diet label if that diet has a significant interaction
# Estimate printed in corner
# ============================================================

# Get KO-bile acid pairs that have at least one significant interaction term
interaction_pairs <- ba_focused %>%
  filter(term_type == "interaction") %>%
  select(KO, bile_acids, estimate, p.value, interaction_diet) %>%
  distinct()

# The significant diet interactions per KO-bile acid (for asterisk labeling)
sig_interactions <- interaction_pairs %>%
  mutate(diet = str_remove(interaction_diet, "^Diet")) %>%  # e.g. "DietXN" -> "XN"
  select(KO, bile_acids, diet, estimate, p.value)

# Join raw data to interaction KO-bile acid pairs
int_ko_ba <- interaction_pairs %>%
  select(KO, bile_acids) %>% distinct()

int_data <- trans_rfKO_BA_tab %>%
  inner_join(int_ko_ba, by = c("KO", "bile_acids"))

# Build estimate label per bile acid (all significant interaction KOs)
int_labels <- sig_interactions %>%
  group_by(bile_acids) %>%
  summarise(
    label = paste(paste0(KO, " (", diet, "): β=", formatC(estimate, format = "e", digits = 2)),
                  collapse = "\n"),
    .groups = "drop"
  )

# Assign shape + linetype per KO
ko_levels_int <- unique(int_data$KO)
ko_aes_int <- data.frame(
  KO       = ko_levels_int,
  ko_shape = shapes_list[seq_along(ko_levels_int) %% length(shapes_list) + 1],
  ko_lty   = linetypes_list[seq_along(ko_levels_int) %% length(linetypes_list) + 1]
)
int_data <- int_data %>% left_join(ko_aes_int, by = "KO")

# Build asterisk annotation positions:
# Place "*" at the max ba_levels value of the significant diet's data,
# at the predicted y from lm(abundance ~ ba_levels) for that KO-bile acid-diet combo
asterisk_data <- sig_interactions %>%
  rename(bile_acids = bile_acids, Diet = diet) %>%
  inner_join(int_data, by = c("KO", "bile_acids", "Diet")) %>%
  group_by(KO, bile_acids, Diet) %>%
  summarise(
    x_pos = max(ba_levels, na.rm = TRUE),
    y_pos = {
      fit <- lm(abundance ~ ba_levels, data = pick(everything()))
      predict(fit, newdata = data.frame(ba_levels = max(ba_levels, na.rm = TRUE)))
    },
    .groups = "drop"
  )

# Build per-panel estimate labels (one per KO-bile acid-diet combo)
# int_labels_grid <- sig_interactions %>%
#   group_by(KO, bile_acids) %>%
#   summarise(
#     label = paste(paste0(diet, ": coeff=", formatC(estimate, format = "e", digits = 2)),
#                   collapse = "\n"),
#     .groups = "drop"
#   )

plot2 <- ggplot(int_data, aes(x = ba_levels, y = abundance, color = Diet)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Diet), linewidth = 0.7) +
  geom_text(data    = asterisk_data,
            mapping = aes(x = x_pos, y = y_pos, color = Diet, label = "*"),
            size = 10, fontface = "bold", hjust = -0.2,
            inherit.aes = FALSE) +
  # geom_text(data    = int_labels_grid,
  #           mapping = aes(label = label),
  #           x = -Inf, y = Inf,
  #           hjust = -0.05, vjust = 1.4,
  #           size = 3, color = "black", inherit.aes = FALSE) +
  facet_wrap(KO ~ bile_acids, scales = "free") +
  coord_cartesian(clip = "off") +
  labs(
    x        = "Bile Acid Level",
    y        = "KO Abundance",
    color    = "Treatment"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text.x    = element_text(face = "bold", size = 14),
    strip.text.y    = element_text(face = "bold", size = 14, angle = 0),
    panel.spacing   = unit(0.3, "lines"),
    text = element_text(size = 16),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )


# ggsave("rf_ba_ko/BA_KO_diet_interactions.tiff", plot2,
#        width  = 18,
#        height = 10,
#        dpi = 300)

## make a list of KO's


################### make rf models with pathways predicting bile acids #############
# used KEGG api to convert the KO sample data frame into a pathway sample data frame of abundances
# now use this to run random forest modeling to see which pathways predict bile acid levels



################### ceramides graphs ########################
# get and clean up this file keaton has of ceramide metadata so I can add to the object below
cerm_metadata <- read.csv("Input/bileAcids_and_lipids_legend.csv") %>%
  dplyr::filter(Label == "C22 ceramide") %>%
  dplyr::select(c("Label", "varID", "Category")) %>%
  rename(Covar = varID)

# get ceramide levels
cermlvl_tab <- read.csv("Input/reformatted_bileAcids_and_lipids_data.csv") %>%
  dplyr::select(c("Sample", cerm_metadata$Covar))

# ceramide random forest with added metadata, filtered for C22, which was found in liver and brain
cerm_rfs <- readRDS("Saved_objects/dt_allCeramide_randForests_imptFeature_sig_interactions.rds") %>%
  dplyr::filter(Set == "KOs") %>%
  inner_join(cerm_metadata, by = "Covar") %>%
  dplyr::filter(Label == "C22 ceramide")

# now filter KO abundance table by these KOs and add ceramide values
cerm_ko_tab <- read.table("Input/Merged_Relative_abundances.tab") %>%
  dplyr::select(cerm_rfs$Feature) %>%
  rownames_to_column("Sample") %>%
  mutate(across(Sample, ~ paste0("s", .x))) %>%
  pivot_longer(cols = cerm_rfs$Feature, names_to = "KO", values_to = "abundance")
cerm_ko_tab2 <- cerm_ko_tab %>%
  inner_join(cermlvl_tab, by = "Sample") %>%
  rename(Brain = var033, Liver = var180) %>%
  pivot_longer(cols = c("Brain", "Liver"), names_to = "ceramide", values_to = "ceramide_levels") %>%
  rename(ID = Sample) %>%
  inner_join(metadata[,1:2], by = "ID")

# stats
# run a quick lm for Diet and ceramide levels and interactions
library(broom)

results_tidy <- cerm_ko_tab2 %>%
  group_by(KO, ceramide) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(abundance ~ ceramide_levels + Diet + ceramide_levels:Diet, data = .)),
    tidied = map(model, tidy),
    glanced = map(model, glance)
  ) %>%
  unnest(tidied, glanced, .drop = TRUE)

print(results_tidy)

# clean up results to only focus on significance at ceramide or ceramideXdiet level
ceramide_focused <- results_tidy %>%
  filter((term == "ceramide_levels" | str_starts(term, "ceramide_levels:Diet")) & p.value < 0.05) %>%
  select(KO, ceramide, term, estimate, std.error, p.value) %>%
  mutate(
    term_type = ifelse(term == "ceramide_levels", "main_effect", "interaction"),
    interaction_diet = ifelse(term == "ceramide_levels", NA,
                              str_extract(term, "Diet[^\\s,]+"))
  ) %>%
  arrange(KO, ceramide, p.value)

print("\nCeramide-focused results (main effects + interactions, excluding Diet-only terms):")
print(ceramide_focused)
# write.csv(ceramide_focused, "Saved_objects/ceramide_focused_significant_results.csv", row.names = FALSE)

# graph abundances of each
ggplot(data = cerm_ko_tab2, aes(x = ceramide_levels, y = abundance, color = ceramide, group = Diet)) +
  geom_point(aes(shape = Diet)) +
  geom_smooth(aes(color = ceramide, linetype = Diet), method = lm) +
  theme_classic() +
  facet_wrap(~ KO + ceramide, scales = "free") +
  labs(x = "C22 ceramide levels", y = "KO abundance", color = "Tissue")

# ggsave("~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/rfKOs_ceramidediet.tiff",
#        width = 12, height = 7,
#        dpi = 300)

# ceramide focused only
cerm_ko_signif <- cerm_ko_tab2 %>%
  dplyr::filter(KO == c("K03338", "K20268"))

ggplot(data = cerm_ko_signif, aes(x = ceramide_levels, y = abundance, color = Diet)) +
  geom_point() +
  geom_smooth(method = lm, se = F) +
  theme_classic() +
  facet_wrap("KO", scales = "free") +
  labs(x = "C22 ceramide levels", y = "KO abundance", color = "Treatment")

# ggsave("~/Documents/Projects/Writing/Research manuscripts/Xanthohumol/XN fig prototypes/rfKOs_ceramide22diet_signif.tiff",
#        width = 7, height = 5,
#        dpi = 300)




