This directory includes R scripts and associated .csv files which are needed to calculate electivity indices per sample. 

#### Data files in this directory and used for this step:
- `dietR_16S_20230502.R` calculates the different indices using `*curated*.csv` files (see below for their specifics)
- `16S_ITS_COMBO_VanderploegScavia_20230525.R` plots V-S indices for both 16S and ITS data in a single script for ease of downstream editing

#### Specific .csv files listed:

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
- `genus_vs_long_plotting_20230522.csv`: output file from R script used for final plotting


---

In `dietR/20230501/electivity_20230522` directory, includes R script (`electivity_16S_20230522.R`) for "pickiness score" analysis, averaging over the top 10 most available resources and associated electivity index per sample by resource (input found here: `16S_electivity_top10avail_20230522.csv`)