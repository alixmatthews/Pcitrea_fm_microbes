This directory includes R scripts and associated .csv files which are needed to calculate electivity indices per sample. 

#### Data files in this directory and used for this step:
- `dietR_16S_20230502.R` calculates the different indices using `*curated*.csv` files (see below for their specifics)
- `16S_ITS_COMBO_VanderploegScavia_20230525.R` plots V-S indices for both 16S and ITS data in a single script for ease of downstream editing

#### Specific .csv files listed:
**NOTE:** on all of these .csv files, the first column is the bird ID (listed as a number instead of the "name", e.g., "P01" is "Ja" from original manifest file), each of the following columns is the microbial taxonomy at different levels (phylum, family, genus), and the matrix is composed of the relative abundance of each taxonomic unit for each sample. 

PHYLUM LEVEL:

- `16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at PHYLUM level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)

FAMILY LEVEL:

- `16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at FAMILY level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)

GENUS LEVEL: 
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at GENUS level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)
- `genus_vs_long_plotting_20230522.csv`: output file from R script used for final plotting. See description of attributes below.


---

In `dietR/20230501/electivity_20230522` directory, includes R script (`electivity_16S_20230522.R`) for "pickiness score" analysis, averaging over the top 10 most available resources and associated electivity index per sample by resource (input found here: `16S_electivity_top10avail_20230522.csv`)


---
## Attribute descriptions
`genus_vs_long_plotting_20230522.csv`: this table provides a reference for understanding the statistical results related to ASV electivity, including measures of mean, variability, confidence intervals, and categorization.

| Variable Name        | Description                                                        | Units            | Notes on Categorical Values             |
|----------------------|--------------------------------------------------------------------|------------------|-----------------------------------------|
| `ASV`                | Amplicon sequence variant (ASV) identifier                         | N/A              | N/A                                     |
| `mean.vs`            | Mean value of V-S electivity for the ASV                               | Mean score       | N/A                                     |
| `sd.vs`              | Standard deviation of V-S electivity for the ASV                       | Score (unitless) | N/A                                     |
| `n.vs`               | Number of samples used to calculate the V-S electivity for the ASV     | Count            | N/A                                     |
| `se.vs`              | Standard error of the mean V-S electivity for the ASV                  | Score (unitless) | N/A                                     |
| `lower.ci.vs`        | Lower bound of the 95% confidence interval for V-S electivity          | Score (unitless) | N/A                                     |
| `upper.ci.vs`        | Upper bound of the 95% confidence interval for V-S electivity          | Score (unitless) | N/A                                     |
| `electivity_category`| Categorization of the electivity based on CI bounds  | N/A              | "against" = mites selected against the ASV, "for" = mites selected for the ASV, "neutral = mites were neutral towards the ASV (CIs crossed 0)           |

---
