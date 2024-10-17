# R analyses for diversity indices and dietary selection

This directory includes the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script for ease of editing downstream.

#### Specific files in this directory:
- `Master_results_20230413.csv`: master dataframe of 16S results of 500+ read results and those using a cutoff of 1451 reads
- `Plotting_16S_ITS_combos_500_data.csv`: master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_500_data_20230903.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 500+ read results
- `Plotting_16S_ITS_combos_rarefied_data.csv`: master dataframe of 16S and ITS results using only 1451+ read results (not used in final analyses, just a comparison)
- `Plotting_16S_ITS_combos_rarefied_data_20230627.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 1451+ read results (not used in final analyses, just a comparison)

---

In `dietR/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample (see more in the README file there)
