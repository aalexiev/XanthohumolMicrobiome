#initialization things
# setwd("/Users/susie/OneDrive/Documents/R/sharpton-lab/XN-feeding/")
library("dplyr")
library("phyloseq")
library("ggplot2")
library("vegan")

#import AUC metadata
meta1 <- read.table("AUC_metadata.tsv",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)
#import Feed Efficiency metadata
meta2 <- read.table("FeedEfficiency_metadata.tsv",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)
#import Spatial Learning metadata
meta3 <- read.table("SpatialLearning_metadata.tsv",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)
#import Treatment metadata
meta4 <- read.table("Treatment_metadata.tsv",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)

#merging metadata tables into a meta-metatable
metadata_table <- merge(meta4, meta1,
                        by = "ID",
                        all = TRUE,
                        sort = TRUE)
metadata_table <- merge(metadata_table, meta2,
                        by = "ID",
                        all = TRUE,
                        sort = TRUE)
metadata_table <- merge(metadata_table, meta3,
                        by = "ID",
                        all = TRUE,
                        sort = TRUE)
#rename ID column to "Sample_ID" to be consistent
colnames(metadata_table)[1] <- "Sample_ID"

#import protein/ko data
ko_metadata <- read.table("Merged_Metadata.tab",
                               sep = "\t",
                               header = TRUE,
                               stringsAsFactors = FALSE)
ko_relative_abundance <- read.table("Merged_Relative_abundances.tab",
                                    sep = " ",
                                    header = TRUE,
                                    row.names = NULL,
                                    stringsAsFactors = FALSE)
#again, make sample ID column read "Sample_ID"
colnames(ko_metadata)[1] <- "Sample_ID"
colnames(ko_relative_abundance)[1] <- "Sample_ID"

#large combined metadata table
metadata_table <- merge(metadata_table, ko_metadata,
                        by = "Sample_ID",
                        all = TRUE,
                        sort = TRUE)

#test histogram
glucose_histo <- hist(metadata_table$gluc,
                      main = "Histogram of Glucose Feed",
                      xlab = "Glucose",
                      border = "black",
                      col = "pink",
                      prob = TRUE)
lines(density(metadata_table$gluc))

LWBW_histo <- hist(metadata_table$LWBW,
                   main = "Histogram of LWBW",
                   xlab = "LWBW",
                   border = "black",
                   col = "light blue",
                   prob = TRUE)

#phyloseq work

ko_abundance_matrix <- data.matrix(ko_relative_abundance)
#convert tables to phyloseq "sub-objects"
ko_table <- otu_table(ko_abundance_matrix, taxa_are_rows = FALSE)
##drop Sample_ID colunmn
drops <- c( "Sample_ID" )
ko_table_clean <- as.data.frame(ko_table)
ko_table_clean <- ko_table_clean[ , !(colnames(ko_table_clean) %in% drops) ]
ko_table <- otu_table(ko_table_clean, taxa_are_rows = FALSE)
#metadata
ko_phylo_meta <- sample_data(metadata_table)
#merge these into phyloseq object
txn_phyloseq <- phyloseq(ko_table, ko_phylo_meta)
saveRDS(txn_phyloseq, file = "Input/phyloseq.rds")
#fractional abundance table - ultimately not used
fraction_txn_phyloseq <- transform_sample_counts(txn_phyloseq, function(OTU) OTU/sum(OTU))

# convert treatment into a factor
## Get our dataframe
df <- data.frame(sample_data(txn_phyloseq))
df$treat <- factor(df$treat)
sample_data(txn_phyloseq) <- df

#ordination and plotting in phyloseq - looking at treatment first
##NMDS and Bray, with treatment
txn.ord <- ordinate(txn_phyloseq, method = "NMDS", distance = "bray")
p1 <- plot_ordination(txn_phyloseq, txn.ord,
                      type = "samples",
                      color = "treat",
                      title = "NMDS+Bray, with treatment")
print(p1)

##PCoA and Bray, with treatment
txn2.ord <- ordinate(txn_phyloseq, method = "PCoA", distance = "bray")
p2 <- plot_ordination(txn_phyloseq, txn2.ord,
                      type = "samples",
                      color = "treat",
                      title = "PCoA+Bray, with treatment")
print(p2)

#NMDS and jaccard
txn3.ord <- ordinate(txn_phyloseq, method = "NMDS", distance = "jaccard")
p3 <- plot_ordination(txn_phyloseq, txn3.ord,
                      type = "samples",
                      color = "treat",
                      title = "NMDS+jaccard, with treatment")
print(p3)
#PCoA and jaccard
txn4.ord <- ordinate(txn_phyloseq, method = "PCoA", distance = "jaccard")
p4 <- plot_ordination(txn_phyloseq, txn4.ord,
                      type = "samples",
                      color = "treat",
                      title = "PCoA+jaccard, with treatment")
print(p4)

#statistical tests - PERMANOVA and envfit
##changing this back to dataframe
df_ko_metadata <- as(sample_data(txn_phyloseq), "data.frame")
##PERMANOVA for treatment
adonis(distance(txn_phyloseq, method = "bray") ~ treat,
       data = df_ko_metadata)
treat_permanova <- adonis(distance(txn_phyloseq, method = "bray") ~ treat,
       data = df_ko_metadata)
View(treat_permanova$aov.tab)
View(treat_permanova$coef.sites)
##envfit - finds LWBW "*" and ptime.x and ptime.y "."
txn_envfit <- envfit(txn.ord, df_ko_metadata, na.rm = TRUE)
View(txn_envfit$vectors$pvals)
##PERMANOVA for LWBW
LWBW_permanova <- adonis(distance(txn_phyloseq, method = "bray") ~ LWBW,
       data = df_ko_metadata)
View(LWBW_permanova$aov.tab)
View(LWBW_permanova$coef.sites)

##envfit for NMDS/jaccard
txn_envfit2 <- envfit(txn3.ord, df_ko_metadata, na.rm = TRUE)
LWBW_permanova2 <- adonis(distance(txn_phyloseq, method = "jaccard") ~ LWBW,
                          data = df_ko_metadata)
View(LWBW_permanova2$aov.tab)

##ordination again, NMDS+Bray, colored by LWBW
p_LWBW <- plot_ordination(txn_phyloseq, txn.ord,
                             type = "samples",
                             color = "LWBW",
                             title = "NMDS+Bray, with LWBW")
print(p_LWBW)
