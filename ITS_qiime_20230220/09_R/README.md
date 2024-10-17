# R analyses for diversity indices and dietary selection

This directory includes the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script, as well as ITS-only csv files (for a few downstream analyses)

#### Specific files in this directory:
- `ITS_forward_master_results_20230413.csv`: master dataframe of ITS results only (500plus and 'rarefied') that excludes a column with Pielou and is used in `10_Maaslin` directory, so need to include here
- `ITS_forward_master_results_20230502.csv`: master dataframe of ITS results only (500plus and 'rarefied') that includes a column with Pielou
- `Plotting_16S_ITS_combos_500_data.csv`: master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_500_data_20230903.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_rarefied_data.csv`: master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison)
- `Plotting_16S_ITS_combos_rarefied_data_20230627.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison)


---

In `dietr/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample (see more in the README file there)
