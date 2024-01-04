# R analyses for diversity indices and dietary selection

In this main directory is the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script for ease of editing downstream.

Specific files in this main directory:
- `Master_results_20230413.csv`: master dataframe of 16S results of 500+ read results and those using a cutoff of 1451 reads
- `Plotting_16S_ITS_combos_500_data.csv`: master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_500_data_20230903.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_rarefied_data.csv`: master dataframe of 16S and ITS results using only 1451+ read results
- `Plotting_16S_ITS_combos_rarefied_data_20230627.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 1451+ read results

---

In `dietR/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample. 

- `dietR_16S_20230502.R` calculates the different indices using `*curated*.csv` files (SEE BELOW)
- `16S_ITS_COMBO_VanderploegScavia_20230525.R` plots V-S indices for both 16S and ITS data in a single script

Specific .csv files listed:

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




