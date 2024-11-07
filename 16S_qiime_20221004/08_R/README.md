# R analyses for diversity indices and dietary selection

This directory includes the plotting of alpha and beta diversity indices, for both 'rarefied' data and 500 seqs+ only data. Includes both ITS and 16S results in a single R script for ease of editing downstream.

#### Specific files in this directory and attribute descriptions:
- `Master_results_20230413.csv`: master dataframe of 16S results of 500+ read results and those using a cutoff of 1451 reads


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
| `shannon_500`                        | Shannon diversity index for 16S sequences at depth 500                | Index (unitless) | Measures species diversity              |
| `faith_pd_500`                       | Faith's phylogenetic diversity for 16S sequences at depth 500         | Index (unitless)         | Measures phylogenetic diversity         |
| `observed_features_500`              | Number of observed features for 16S sequences at depth 500            | Count            | N/A                                     |
| `chao_500`                           | Chao1 estimator of species richness for 16S sequences at depth 500    | Count            | Estimator including rare species        |
| `braycurtis_pc1_500`                 | Bray-Curtis principal component 1 for 16S sequences at depth 500      | Arbitrary units  | N/A                                     |
| `braycurtis_pc2_500`                 | Bray-Curtis principal component 2 for 16S sequences at depth 500      | Arbitrary units  | N/A                                     |
| `braycurtis_pc3_500`                 | Bray-Curtis principal component 3 for 16S sequences at depth 500      | Arbitrary units  | N/A                                     |
| `jaccard_pc1_500`                    | Jaccard principal component 1 for 16S sequences at depth 500         | Arbitrary units  | N/A                                     |
| `jaccard_pc2_500`                    | Jaccard principal component 2 for 16S sequences at depth 500         | Arbitrary units  | N/A                                     |
| `jaccard_pc3_500`                    | Jaccard principal component 3 for 16S sequences at depth 500         | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc1_500`          | Unweighted UniFrac PC1 for 16S sequences at depth 500                 | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc2_500`          | Unweighted UniFrac PC2 for 16S sequences at depth 500                 | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc3_500`          | Unweighted UniFrac PC3 for 16S sequences at depth 500                 | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc1_500`            | Weighted UniFrac PC1 for 16S sequences at depth 500                   | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc2_500`            | Weighted UniFrac PC2 for 16S sequences at depth 500                   | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc3_500`            | Weighted UniFrac PC3 for 16S sequences at depth 500                   | Arbitrary units  | N/A                                     |
| `rpca_pc1_500`                       | RPCA (PC1) for 16S sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `rpca_pc2_500`                       | RPCA (PC2) for 16S sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `rpca_pc3_500`                       | RPCA (PC3) for 16S sequences at depth 500                             | Arbitrary units  | N/A                                     |
| `shannon_1451`                       | Shannon diversity index for 16S sequences at depth 1451               | Index (unitless) | Measures species diversity              |
| `faith_pd_1451`                      | Faith's phylogenetic diversity for 16S sequences at depth 1451        | Index (unitless)         | Measures phylogenetic diversity         |
| `observed_features_1451`             | Number of observed features for 16S sequences at depth 1451           | Count            | N/A                                     |
| `braycurtis_pc1_1451`                | Bray-Curtis principal component 1 for 16S sequences at depth 1451     | Arbitrary units  | N/A                                     |
| `braycurtis_pc2_1451`                | Bray-Curtis principal component 2 for 16S sequences at depth 1451     | Arbitrary units  | N/A                                     |
| `braycurtis_pc3_1451`                | Bray-Curtis principal component 3 for 16S sequences at depth 1451     | Arbitrary units  | N/A                                     |
| `jaccard_pc1_1451`                   | Jaccard principal component 1 for 16S sequences at depth 1451        | Arbitrary units  | N/A                                     |
| `jaccard_pc2_1451`                   | Jaccard principal component 2 for 16S sequences at depth 1451        | Arbitrary units  | N/A                                     |
| `jaccard_pc3_1451`                   | Jaccard principal component 3 for 16S sequences at depth 1451        | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc1_1451`         | Unweighted UniFrac PC1 for 16S sequences at depth 1451                | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc2_1451`         | Unweighted UniFrac PC2 for 16S sequences at depth 1451                | Arbitrary units  | N/A                                     |
| `unweightedunifrac_pc3_1451`         | Unweighted UniFrac PC3 for 16S sequences at depth 1451                | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc1_1451`           | Weighted UniFrac PC1 for 16S sequences at depth 1451                  | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc2_1451`           | Weighted UniFrac PC2 for 16S sequences at depth 1451                  | Arbitrary units  | N/A                                     |
| `weightedunifrac_pc3_1451`           | Weighted UniFrac PC3 for 16S sequences at depth 1451                  | Arbitrary units  | N/A                                     |
| `rpca_pc1_1451`                      | RPCA (PC1) for 16S sequences at depth 1451                            | Arbitrary units  | N/A                                     |
| `rpca_pc2_1451`                      | RPCA (PC2) for 16S sequences at depth 1451                            | Arbitrary units  | N/A                                     |
| `rpca_pc3_1451`                      | RPCA (PC3) for 16S sequences at depth 1451                            | Arbitrary units  | N/A                                     |



- `Plotting_16S_ITS_combos_500_data.csv`: master dataframe of 16S and ITS results using only 500+ read results.

| Variable Name            | Description                                                        | Units           | Notes on Categorical Values             |
|--------------------------|--------------------------------------------------------------------|-----------------|-----------------------------------------|
| `sample-id`                | Unique identifier for each sample                                 | N/A             | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`           | Indicates if sample is biological or a control                    | N/A             | `biological` for biological, `control` for control samples |
| `category_broad`           | Broad category for sample classification                          | N/A             | e.g., `mites`, `feather`, `control`  |
| `category_specific`        | Specific category for sample classification                       | N/A             | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`          |
| `bird_id`                  | Unique identifier for each bird                                   | N/A             | "names" for each bird in the study (easier and more fun way to keep track)                                     |
| `shannon_500_16s`        | Shannon diversity index for 16S sequences at read depth of 500         | Index (unitless)| Measures species diversity              |
| `chao_500_16s`           | Chao1 estimator of species richness for 16S sequences at read depth of 500| Count           | Estimates total richness including rare species |
| `RPCA_pc1_500_16S`       | RPCA component 1 (PC1) score using a read cutoff of 500 (bacteria, 16S) | Arbitrary units | N/A                                     |
| `RPCA_pc2_500_16S`       | RPCA component 2 (PC2) score using a read cutoff of 500 (bacteria, 16S)                   | Arbitrary units | N/A                                     |
| `shannon_500_its`        | Shannon diversity index for ITS sequences at read depth of 500         | Index (unitless)| Measures species diversity              |
| `chao_500_its`           | Chao1 estimator of species richness for ITS sequences at read depth 500| Count           | Estimates total richness including rare species |
| `RPCA_pc1_500_its`       | RPCA component 1 (PC1) score using a read cutoff of 500 (fungi, ITS)                   | Arbitrary units | N/A                                     |
| `RPCA_pc2_500_its`       | RPCA component 2 (PC2) score using a read cutoff of 500 (fungi, ITS)                    | Arbitrary units | N/A                                     |



- `Plotting_16S_ITS_combos_500_data_20230903.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 500+ read results

  
- `Plotting_16S_ITS_combos_rarefied_data.csv`: master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison).


| Variable Name              | Description                                                        | Units           | Notes on Categorical Values             |
|----------------------------|--------------------------------------------------------------------|-----------------|-----------------------------------------|
| `sample-id`                | Unique identifier for each sample                                 | N/A             | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`           | Indicates if sample is biological or a control                    | N/A             | `biological` for biological, `control` for control samples |
| `category_broad`           | Broad category for sample classification                          | N/A             | e.g., `mites`, `feather`, `control`  |
| `category_specific`        | Specific category for sample classification                       | N/A             | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`          |
| `bird_id`                  | Unique identifier for each bird                                   | N/A             | "names" for each bird in the study (easier and more fun way to keep track)                                     |
| `shannon_1451`             | Shannon diversity index for bacteria (16S) using a read cutoff value of 1451                              | Index (unitless)| Measures species diversity              |
| `observed_features_1451`   | Number of observed features (ASVs) for bacteria (16S) using a read cutoff value of 1451                   | Count           | N/A                                     |
| `RPCA_pc1_0_16S`           | RPCA component 1 (PC1) score using a read cutoff of 0 (bacteria, 16S)   | Arbitrary units | N/A                                     |
| `RPCA_pc2_0_16S`           | RPCA component 2 (PC2) score using a read cutoff of 0 (bacteria, 16S)                                  | Arbitrary units | N/A                                     |
| `RPCA_pc1_1451`            | RPCA component 1 (PC1) score using a read cutoff of 1451 (bacteria, 16S)                               | Arbitrary units | N/A                                     |
| `RPCA_pc2_1451`            | RPCA component 2 (PC2) score using a read cutoff of 1451 (bacteria, 16S)                              | Arbitrary units | N/A                                     |
| `shannon_995`              | Shannon diversity index for fungi (ITS) using a read cutoff value of 995                              | Index (unitless)| Measures species diversity              |
| `observed_features_995`    | Number of observed features (ASVs) for fungi (ITS) using a read cutoff value of 995                    | Count           | N/A                                     |
| `RPCA_pc1_0_ITS`           | RPCA component 1 (PC1) score using a read cutoff of 0 (fungi, ITS)                                  | Arbitrary units | N/A                                     |
| `RPCA_pc2_0_ITS`           | RPCA component 2 (PC2) score using a read cutoff of 0 (fungi, ITS)                                  | Arbitrary units | N/A                                     |
| `RPCA_pc1_995`             | RPCA component 1 (PC1) score using a read cutoff of 995 (fungi, ITS)                                | Arbitrary units | N/A                                     |
| `RPCA_pc2_995`             | RPCA component 2 (PC2) score using a read cutoff of 995 (fungi, ITS)                                | Arbitrary units | N/A                                     |


  
- `Plotting_16S_ITS_combos_rarefied_data_20230627.R`: R script which plots the data from the master dataframe of 16S and ITS results using only 1451+ and 995+ read results (not used in final analyses, just a comparison)

---

In `dietR/20230501` directory, includes R scripts and associated .csv files which are needed to calculate electivity indices per sample (see more in the README file there)

