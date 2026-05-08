Linking metagenome composition and behavioral responses
================
Keaton Stagaman
2020-05-01

## Background

Xanthohumol (XN), a flavonoid produced by hops, has been shown to
mitigate the effects of metabolic syndrome due to high-fat diets (HFD)
in various animal models. However, XN can spontaneously form a stable
isomer, isoxanthohumol (IX). Both animals and gut microbiome
constituents produce enzymes that can transform IX into
8-prenylnaringenin (8-PN), the most potent phytoestrogen known to date.
Alternatives to XN, i.e., hydrogenated derivatives of XN: α,β-dihydro-XN
(DXN) and tetrahydro-XN (TXN), show negligible affinity for estrogen
receptors, and cannot be metabolically converted into 8-PN. These
compounds have been shown to have similar effects on metabolic syndrome,
in including improving spatial learning outcomes in mice. Here, we
endeavor to determine how XN, DXN, and TXN influence the functional
potential (metagenome) of mice being fed high-fat diets, and whether
these such effects associate with specific spatial learning outcomes. We
assigned shotgun metagenomic sequences using two methods: the first was
to align them to reference sequences (KOs) in the Kyoto Encyclopedia of
Gene and Genomes (KEGG) database; the second was to align the sequences
to a integrated gene catalog (IGC) of the mouse metagenomes. The
following analyses utilize both of these data sets.

## Effects of Diet alone

### Differences in overall composition (beta-diversity)

![](Plots/diet_ordinations.png)

![](Plots/diet_permanova_results.png)

The presence of XN or its derivatives affects the composition of the
potential functional capacity of the gut microbiome in mice (Figure 1 &
Table 1). The distinctions appear to be clearer when using the IGC data
set, rather than the KO set. In either case, however, XN treatment
appears to result in metagenomic composition more similar to the
controls than DXN or TXN treatments do. It *appears* that the control
(HFD) metagenomes have much less beta-dispersion (average distance from
the centroid) than the diet treatments.

### Taxonomic (Kraken2) composition

![](Plots/diet_taxonomic_ordination.png)

![](Plots/kraken_permanova_results.png)

We used Kraken 2 to assign taxonomy to the mouse metagenomic sequences
to assess whether administration of XN and its derivatives affects not
just the potential functional capacity of the mouse gut microbiome, but
also which microbial taxa are present. While there are statistically
significant differences for all three distance metics measured (Table
2), the differences in the ordinations indicate some interesting
results. Bray-Curtis, which is weighted by abundance, implies that there
many high-abundance taxa share across the treatments, while Sørensen
(presence-absence) indicates a strong effect of diet treatment on the
rarer taxa. That is to say, administration of XN, DXN, and TXN to mice
being fed HFDs appears to primarily affect the presence or absence of
rare taxa, and have a smaller (but still stastically signficant) effect
on the more abundant taxa. As with the functional potential, taxonomic
composition of mice supplemented with XN apppears to be more similar to
the HFD controls than supplementation with DXN or TXN.

### Associations with individual functional annotations

![](Plots/diet_rf_sig_impt_features.png)

![](Plots/diet_rf_sig_impt_kos_tbl.png)

After to assessing the overall composition of the mouse metagenomes, we
wanted to determine if there association with diet treatment and
specific metagenomic functions. To do so, we created two random forest
models, the first using KO abundances and the second using IGC
abundances, to predict diet treatment. We assessed feature (KO or IGC
abudance) importance through a non-parametric (permutational) method and
ploted the abundance-by-diet relationships of the top 12 significantly
important features for each data set (Figure 3). Higher level (module)
assignments, if determined, are presented in Table 3.

<!-- fix code and fill out text (KOs previously shown to associate with XN) -->

![](Plots/diet_xn_kos_of_interest.png)

![](Plots/diet_xn_kos_of_interest_tbl.png)

Previous work had identified particular KOs that associated with XN diet
supplementation. We search for KOs in this work that matched any of a
set of keywords from that prior work and plotted their abundances by
diet treatment.

## Diet and Spatial Learning covariate interactions

### Differences in overall composition (beta-diversity)

#### Methods

1.  Generate dbRDA ordinations with all Spatial Learning covariates
2.  Use `ordistep` to select Spatial Learning covariates that explain
    the most variance in beta-diversity
3.  Significance assessment on `ordistep`-selected dbRDAs with permanova
4.  Add in diet term (main effect)
5.  More tests of significance
6.  Add in diet term interactions to full models
7.  Use `ordistep` to select Spatial Learning covariates by Diet
    interactions that explain the most variance in beta-diversity
8.  Significance assessment on `ordistep`-selected dbRDAs with permanova

![](Plots/diet_by_covars_ordinations.png)

![](Plots/diet_by_covars_permanova_results.png)

Distance-based redundancy analysis (dbRDA) revealed a number of
significant interactions between diet treatment and spatial learning
covariates in predicting differences in metagenomic composition (Figure
4 & Table 4). Significant interactions, in this case, implies that the
association between a given spatial learning covariate score and
metagenomic composition is dependent on which diet treatment (control,
XN, DXN, or TXN) the mouse received. Of note, selection on models
predicting Sørensen (presence-absence) scores return no significant
interactions between diet and spatial learning covariates, while
Bray-Curtis and Canberra (both abundance-weighted) did. This implies
that the spatial learning covariates associate with changes in abundance
of ceratain (probably overall abundant) microbial taxa, rather than the
presence or absence of particular microbial taxa.

### Associations with individual functional annotations

As with diet treatment, we created a set of random forest models to
predict spatial learning covariate scores from feature (KO or IGC)
abundances. We then used the models to assess which features were
significanly important. Using only these important features, we built
regression models to assess whether feature abundance by diet
interactions significantly predicted spatial learning covariate scores.
As the table and figure below indicate, there were a number of IGCs and
KOs that significantly associated with spatial learning covariates in a
diet-dependent manner. For example, K21600 (csoR, ricR; CsoR family
transcriptional regulator, copper-sensing transcriptional repressor) had
a positive association with L_Vis1 scores in control and XN-supplemented
mice, but a negative association with L_Vis1 scores in DXN- and
TXN-supplmented mice, and K00566 (mnmA, trmU; tRNA-uridine
2-sulfurtransferase) had a negative association with L_Hid6 scores in
control samples, but positive associations with L_Hi6 scores for all
three supplemented samples.

This [CSV table](Plots/diet_by_covars_rf_sig_impt_features_table.csv)
reports all significant associations between metagenomic functional
features and behavioral covariate by diet interactions ([PNG
version](Plots/diet_by_covars_rf_sig_impt_features_table.png))

All significant associations between spatial learning covariates and
diet by KO/IGC abundance interactions can be found in [this
directory](Plots/Diet_by_covar_rf_sig_impt_features_plots).

## Diet and ceramide interactions

One possible contribution to the effects of HFDs on metabolic syndrome
is the increased concentrations of sphingolipids (including ceramides)
in the brain and liver. We measured ceramide concentrations in both of
these tissues of the mouse samples, and determined whether there were
significant associtation between these concentrations and metagenomic
composition in a diet-dependent manner.

### Differences in overall composition (beta-diversity)

*Same methodology as above*

![](Plots/diet_by_ceramide_ordinations.png)

![](Plots/diet_by_ceramide_permanova_results.png)

The KOs dataset yielded the most significant results between metagenomic
composition and ceramide concentration by diet treatment interactions
(Figure 6 & Table 6). Specifically brain ceramides C20, C22, C22:1, and
C24; and liver ceramides C16 and C22 exhibited significant associations
with differences in KO composition between samples in a diet-dependent
manner. Brain ceramides C20 and C22, and liver ceramide C16 showed
significant associations with differences in IGC composition of the gut
microbiome.

## Summary

Overall, we have discovered that supplementing a high fat diet with
xanthohumol, or either of its derivatives that we examined here, affects
the composition of both the potential functional capacity and taxonomy
of the gut microbiome in mice. Diet supplementation has the greatest
predictive power of any of the covariates we examined in this study
(Supplemental Tables). Furthermore, we were able to identify
associations between both spatial learning covariates and brain and
liver ceramide concentrations and the abundances of specific functions
(KOs and IGCs) in the gut metagenome. The associations we highlighted
here may be of particular interest as they were dependent on the
particular dietary supplement supplied to the mice, indicating that
related molecules like XN, DXN, and TXN can have significanly different
effects on the functional potential of gut microbiome as well as
concentrations of biologically important compounds in host tissues.
