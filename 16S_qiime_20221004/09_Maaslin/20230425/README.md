Maaslin analyses, results, and figures

#### Data files in this directory and used for this step (and attribute descriptions):

- `MaAsLin2_16S_20230425.R`: maaslin analyses in R for 16S
- `16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv`: phyla level input file for analyses
- `16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv`: family level input file for analyses
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv`: genus level input file for analyses
**NOTE:** For these three files above, the first column is the sample ID, followed by columns of the taxonomic units at different levels of taxonomy (e.g., level-2 file is the phyla level). Values within the matrix are the relative abundance of each of the taxonomic units for each sample.

- `16S_ITS_COMBO_RESULTS_phyla_family_genus.xlsx`: results file, includes 16S and ITS results (all results at the phyla, family, and genus levels)
**NOTE:** This file has 6 different tabs (for each of the taxonomic levels)

| Variable Name        | Description                                                      | Units            | Notes on attributes           |
|----------------------|------------------------------------------------------------------|------------------|-----------------------------------------|
| `feature`            | Name or identifier of the feature under analysis                 | N/A              | taxonomic level differs for each tab                                     |
| `metadata`           | Metadata associated with the feature                             | N/A              | e.g., category or sample information    |
| `value`              | Value of the measured feature or statistic                       | N/A              | always "mites"                                     |
| `coef`               | Coefficient from statistical model | Numeric (unitless)| N/A                                     |
| `stderr`             | Standard error of the coefficient                                 | Numeric          | N/A                                     |
| `N`                  | Number of samples used in the analysis                           | Count            | N/A                                     |
| `N.not.0`            | Number of samples with non-zero values for the feature          | Count            | N/A                                     |
| `pval`               | P-value for statistical significance                            | p-value (unitless) | N/A                                    |
| `qval`               | Adjusted p-value (e.g., FDR correction)    | p-value (unitless) | N/A                                    |

  
- `16S_ITS_COMBO_MaAsLin2_figures.R`: making the figures file, includes 16S and ITS figure-making
  
- `16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2_PREVALENCEABUNDANCE.xlsx`: calculations of prevalence and abundance at the genus level (which was x-y transformed from the genus level input files for analysis) then separated by mites/feathers to calculate prev/abund, n, mean, min, max, etc. of each genus in each group

- `differentially_abundant_16S_ITS_PREVALENCEABUNDANCE.xlsx`: differentially abundant summary based on maaslin results
  
