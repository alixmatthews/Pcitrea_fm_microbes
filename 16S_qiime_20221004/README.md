This directory contains files and subdirectories for the 16S data analyses (cleaned up from private repos). Below is a description of the file to begin the analyses (`16S_sample_metadata.tsv`)

#### Files in this directory
- `16S_sample_metadata.tsv` is the original metadata file. Modified in some downstream steps (e.g., for plotting)

#### Variable descriptions
`16S_sample_metadata.tsv`

| Variable Name       | Description                                      | Units         | Notes             |
|---------------------|--------------------------------------------------|---------------|-----------------------------------------|
| `sample-id`         | Unique identifier for each sample                | N/A           | Samples that are paired share the same number value (e.g., PM01 and PF01 were mites and the associated feather, respectively, collected from the same bird)                                     |
| `bio_or_control`    | Indicates if sample is biological or a control   | N/A           | `biological` for biological, `control` for control samples |
| `category_broad`    | Broad category for sample classification         | N/A           | e.g., `mites`, `feather`, `control` |
| `category_specific` | Specific category for sample classification      | N/A           | e.g., `mites`, `feather`, `wash_control`, `kit_control`, `area_control`, `water_control`, `pcr_control`         |
| `bird_id`           | Unique identifier for each bird                  | N/A           | "names" for each bird in the study (easier and more fun way to keep track)                                     |
| `bird_sex`          | Sex of the bird                                  | N/A           | `Male` for male, `Female` for female            |
| `bird_age`          | Age of the bird                                  | N/A         | e.g., `ASY` for after-second year; `SY` for second year                    |
| `capture_date`      | Date the bird was captured                       | MM/DD/YY    | N/A                                     |
| `num_mites_fm`      | Number of mites observed on the feather      | Count         | feather (PF*) and mite (PM*) samples from the same bird have the same value here                                     |
| `num_mites_m`       | Number of mites observed on the feather       | Count         | only mite (PM*) samples have a value here (feather [PF*] samples have "na")                                  |
| `mite_sp_fm`        | Species (genus level) of mites found on the feather        | N/A           | feather (PF*) and mite (PM*) samples from the same bird have the same value here, is either `Amerodectes`, `Proctophyllodes`, or `Amerodectes and Proctophyllodes`                          |
| `mite_sp_m`         | Species (genus level) of mites found on the feather          | N/A           | only mite (PM*) samples have a value here (feather [PF*] samples have "na"), is either `Amerodectes`, `Proctophyllodes`, or `Amerodectes and Proctophyllodes`                      |

---
