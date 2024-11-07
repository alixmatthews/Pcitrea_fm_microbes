This directory includes R script and associated .csv file which are used for "pickiness score" analysis

#### Data files in this directory and used for this step:

- `16S_electivity_top10avail_20230522.csv`: top 10 most available resources and associated electivity index per sample by resource. See below for more details on this table.
- `electivity_16S_20230522.R`: R script for "pickiness score" analysis


| Variable Name                                                        | Description                                                            | Units       | Notes on Categorical Values               |
|----------------------------------------|---------------------------------------------------------|-------------|-------------------------------------------|
| `sample-id`                                                          | Unique identifier for each sample                                       | N/A         | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                       |
| `bird_sex`                                                           | Sex of the bird                                                         | N/A         | `Male` for male, `Female` for female              |
| `bird_age`                                                           | Age of the bird                                                         | Years       | e.g., `ASY` for after-second year; `SY` for second year                       |
| `num_mites_feather`                                                  | Number of mites observed on the individual feather that was collected (that mites were living on)                                    | Count       | N/A                                       |
| `num_mites_total`                                                    | Total number of mites observed across the bird                          | Count       | N/A                                       |
| `paired_for_electivity`                                              | Indicates if the sample was paired for electivity testing               | N/A         | `yes` or `no`, though these should all be `yes` since this is the final test                             |
| `VS_electivity_*`                                                     | V-S electivity score for each taxonomic unit (top 10 most available microbes) | Score | represent taxonomic classifications at the genus level. These columns contain the V-S score of specific taxa detected in the samples during the VS electivity analyses. |
| `VS_electivity_AVERAGE`                                                     | Average of the V-S electivity scores of the top 10 most available microbes | Score | N/A |



