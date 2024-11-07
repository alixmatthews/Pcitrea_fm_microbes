# R analyses for diversity indices and dietary selection

This directory includes the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script, as well as ITS-only csv files (for a few downstream analyses)

#### Specific files in this directory and attribute descriptions:
- `ITS_forward_master_results_20230413.csv`: master dataframe of ITS results only (500plus and 'rarefied') and is used in `10_Maaslin` directory

| Variable Name                        | Description                                                           | Units            | Notes on Categorical Values             |
|--------------------------------------|-----------------------------------------------------------------------|------------------|-----------------------------------------|
| `sample-id`                          | Unique identifier for each sample                                     | N/A              | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`                     | Indicates if sample is biological or a control                        | N/A              | `biological` for biological, `control` for control samples |
| `category_broad`                     | Broad category for sample classification                              | N/A              | e.g., `mites`, `feather`, `control` |
| `category_specific`                  | Specific category for sample classification                           | N/A              | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`         |
| `bird_id`                            | Unique identifier for each bird                                       | N/A              | N/A                                     |
| `bird_sex`                           | Sex of the bird                                                       | N/A              | `Male` for male, `Female` for female            |
| `bird_age`                           | Age of the bird                                                       | Years            | `ASY` is after-second year, `SY` is second year                    |
| `capture_date`                       | Date when the sample was collected                                    | Date (MM/DD/YYYY) | N/A                                     |
| `num_mites_feather`                  | Number of mites observed on feathers                                  | Count            | N/A                                     |
| `num_mites_total`                    | Total number of mites observed across the bird                        | Count            | N/A                                     |
| `mite_sp_fm`                         | Genus of mites found on the feathers                               | N/A              |       |
| `shannon_500`                        | Shannon diversity index for ITS sequences at depth 500                | Index (unitless) | Measures species diversity              |
| `observed_features_500`              | Number of observed features for ITS sequences at depth 500            | Count            | N/A                                     |
| `chao_500`                           | Chao1 estimator of species richness for ITS sequences at depth 500    | Count            | Estimator including rare species        |
| `pielou_500`                         | Pielou’s evenness index for ITS sequences at depth 500                | Index (unitless) | Measures evenness of species distribution|
| `braycurtis_pc1_500`                 | Bray-Curtis principal component 1 for ITS sequences at depth 500      | Arbitrary units  | N/A                                     |
| `braycurtis_pc2_500`                 | Bray-Curtis principal component 2 for ITS sequences at depth 500      | Arbitrary units  | N/A                                     |
| `braycurtis_pc3_500`                 | Bray-Curtis principal component 3 for ITS sequences at depth 500      | Arbitrary units  | N/A                                     |
| `jaccard_pc1_500`                    | Jaccard principal component 1 for ITS sequences at depth 500         | Arbitrary units  | N/A                                     |
| `jaccard_pc2_500`                    | Jaccard principal component 2 for ITS sequences at depth 500         | Arbitrary units  | N/A                                     |
| `jaccard_pc3_500`                    | Jaccard principal component 3 for ITS sequences at depth 500         | Arbitrary units  | N/A                                     |
| `rpca_pc1_500`                       | RPCA (PC1) for ITS sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `rpca_pc2_500`                       | RPCA (PC2) for ITS sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `rpca_pc3_500`                       | RPCA (PC3) for ITS sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `shannon_995`                        | Shannon diversity index for ITS sequences at depth 995                | Index (unitless) | Measures species diversity              |
| `faith_pd_995`                       | Faith's phylogenetic diversity for ITS sequences at depth 995         | PD units         | Measures phylogenetic diversity         |
| `observed_features_995`              | Number of observed features for ITS sequences at depth 995            | Count            | N/A                                     |
| `braycurtis_pc1_995`                 | Bray-Curtis principal component 1 for ITS sequences at depth 995      | Arbitrary units  | N/A                                     |
| `braycurtis_pc2_995`                 | Bray-Curtis principal component 2 for ITS sequences at depth 995      | Arbitrary units  | N/A                                     |
| `braycurtis_pc3_995`                 | Bray-Curtis principal component 3 for ITS sequences at depth 995      | Arbitrary units  | N/A                                     |
| `jaccard_pc1_995`                    | Jaccard principal component 1 for ITS sequences at depth 995         | Arbitrary units  | N/A                                     |
| `jaccard_pc2_995`                    | Jaccard principal component 2 for ITS sequences at depth 995         | Arbitrary units  | N/A                                     |
| `jaccard_pc3_995`                    | Jaccard principal component 3 for ITS sequences at depth 995         | Arbitrary units  | N/A                                     |


- `ITS_forward_master_results_20230502.csv`: same as `ITS_forward_master_results_20230413.csv` in other words, master dataframe of ITS results only (500plus and 'rarefied'), but this dataset also includes a column with Pielou

- `Plotting_16S_ITS_combos_500_data.csv`: master dataframe of 16S and ITS results using only 500+ read results

| Variable Name                        | Description                                                           | Units            | Notes on Categorical Values             |
|--------------------------------------|-----------------------------------------------------------------------|------------------|-----------------------------------------|
| `sample-id`                          | Unique identifier for each sample                                     | N/A              | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`                     | Indicates if sample is biological or a control                        | N/A              | `biological` for biological, `control` for control samples |
| `category_broad`                     | Broad category for sample classification                              | N/A              | e.g., `mites`, `feather`, `control` |
| `category_specific`                  | Specific category for sample classification                           | N/A              | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`         |
| `bird_id`                            | Unique identifier for each bird                                       | N/A              | N/A                                     |
| `shannon_500_16s`                    | Shannon diversity index for 16S sequences at depth 500                | Index (unitless) | Measures species diversity              |
| `chao_500_16s`                       | Chao1 estimator of species richness for 16S sequences at depth 500    | Count            | Estimator including rare species        |
| `RPCA_pc1_500_16S`                   | RPCA (PC1) for 16S sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `RPCA_pc2_500_16S`                   | RPCA (PC2) for 16S sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `shannon_500_its`                    | Shannon diversity index for ITS sequences at depth 500                | Index (unitless) | Measures species diversity              |
| `chao_500_its`                       | Chao1 estimator of species richness for ITS sequences at depth 500    | Count            | Estimator including rare species        |
| `RPCA_pc1_500_its`                   | RPCA (PC1) for ITS sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `RPCA_pc2_500_its`                   | RPCA (PC2) for ITS sequences at depth 500                             | Arbitrary units  | N/A                                     |



- `Plotting_16S_ITS_combos_500_data_20230903.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 500+ read results


- `Plotting_16S_ITS_combos_rarefied_data.csv`: master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison)

| Variable Name                        | Description                                                           | Units            | Notes on Categorical Values             |
|--------------------------------------|-----------------------------------------------------------------------|------------------|-----------------------------------------|
| `sample-id`                          | Unique identifier for each sample                                     | N/A              | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`                     | Indicates if sample is biological or a control                        | N/A              | `biological` for biological, `control` for control samples |
| `category_broad`                     | Broad category for sample classification                              | N/A              | e.g., `mites`, `feather`, `control` |
| `category_specific`                  | Specific category for sample classification                           | N/A              | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`         |
| `bird_id`                            | Unique identifier for each bird                                       | N/A              | N/A                                     |
| `shannon_1451`                       | Shannon diversity index for 16S sequences at depth 1451               | Index (unitless) | Measures species diversity              |
| `observed_features_1451`             | Number of observed features for 16S sequences at depth 1451           | Count            | N/A                                     |
| `RPCA_pc1_0_16S`                     | RPCA (PC1) for 16S sequences at depth 0                                | Arbitrary units  | N/A                                     |
| `RPCA_pc2_0_16S`                     | RPCA (PC2) for 16S sequences at depth 0                                | Arbitrary units  | N/A                                     |
| `RPCA_pc1_1451`                      | RPCA (PC1) for 16S sequences at depth 1451                             | Arbitrary units  | N/A                                     |
| `RPCA_pc2_1451`                      | RPCA (PC2) for 16S sequences at depth 1451                             | Arbitrary units  | N/A                                     |
| `shannon_995`                        | Shannon diversity index for ITS sequences at depth 995                | Index (unitless) | Measures species diversity              |
| `observed_features_995`              | Number of observed features for ITS sequences at depth 995            | Count            | N/A                                     |
| `RPCA_pc1_0_ITS`                     | RPCA (PC1) for ITS sequences at depth 0                                | Arbitrary units  | N/A                                     |
| `RPCA_pc2_0_ITS`                     | RPCA (PC2) for ITS sequences at depth 0                                | Arbitrary units  | N/A                                     |
| `RPCA_pc1_995`                       | RPCA (PC1) for ITS sequences at depth 995                              | Arbitrary units  | N/A                                     |
| `RPCA_pc2_995`                       | RPCA (PC2) for ITS sequences at depth 995                              | Arbitrary units  | N/A                                     |



- `Plotting_16S_ITS_combos_rarefied_data_20230627.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison)


---

In `dietr/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample (see more in the README file there)
