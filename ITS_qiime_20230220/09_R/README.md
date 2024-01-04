# R analyses for diversity indices and dietary selection

In main directory is the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script. Also includes `ITS_forward_master_results_20230502.csv`, which is a master dataframe of ITS results only (500plus and 'rarefied') that includes a column with Pielou, whereas `ITS_forward_master_results_20230413.csv` is identical but does not include Pielou column (but is used in `10_Maaslin` so need to include here)

Another R script is `Plotting_16S_ITS_combos_500_data_20230903.R`, which is 16S and ITS plotting combos using only data with 500+ seqs, requiring `Plotting_16S_ITS_combos_500_data.csv` as input.

A final R script herein is `Plotting_16S_ITS_combos_rarefied_data_20230627.R` which is 16S and ITS combo pllotting using 'rarefied' data, requiring `Plotting_16S_ITS_combos_rarefied_data.csv` as input.


---


In `dietR/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample. 

- `dietR_ITS_20230514.R` calculates the different indices using the `*curated*.csv` files (SEE BELOW)
- `16S_ITS_COMBO_VanderploegScavia_20230525.R` plots V-S indices for both 16S and ITS data in a single script



Specific .csv files listed:

PHYLUM LEVEL:

- `1ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at PHYLUM level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)

FAMILY LEVEL:

- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at FAMILY level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)

GENUS LEVEL: 
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2.csv`: Relative abundance of 500+ seq samples at GENUS level, dataset has been curated for MaAsLin2 (see `07_TaxaBarplots/README.md`) for more information on how this file was created.
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`: Same as above, but only includes FEATHER (available) samples that were paired (by SampleID number)
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv`: Same as above, but only includes MITES (consumed) samples that were paired (by SampleID number)
- `genus_vs_long_plotting_20230522.csv` is the output from previous R script is needed for input in following R script (ITS data only)



---


In `dietR/20230501/electivity_20230522` directory, includes R script (`electivity_ITS_20230522.R`) for "pickiness score" analysis, averaging over the top 10 most available resources and associated electivity index per sample by resource (input found here: `ITS_electivity_top10avail_20230522.csv`)







