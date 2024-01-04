# Maaslin analyses, results, and figures

All conducted in R. Everything within the `./20230425` subdirectory.

Analyses: `MaAsLin2_16S_20230425.R`
Inputs for analyses: 
- `16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv`: phyla level
- `16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv`: family level
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv`: genus level


Results: `16S_ITS_COMBO_RESULTS_phyla_family_genus.xlsx`
Includes 16S and ITS results (all results at the phyla, family, and genus levels)

Figures: `16S_ITS_COMBO_MaAsLin2_figures.R`
Includes 16S and ITS figure-making

Calculations of prevalence and abundance at the genus level can be found in the files below: 
Dataframe has been x-y transformed from `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv` and then separated by mites/feathers to calculate prevalence and abundance, n, etc. of each genus for each group: `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2_PREVALENCEABUNDANCE.xlsx`... differentially abundant ones based on Maaslin results are summarized here: `differentially_abundant_16S_ITS_PREVALENCEABUNDANCE.xlsx`
