# Decontam

Here we are running a decontam pipeline in R.

#### Data files here:
- `decontam_20221005.R`: decontam R pipeline
- `16S_filtered-all_samples-feature_table-with-taxonomy.biom`: biom file with all samples, features, and taxonomy
- `DecontamResults_prev05.csv`: results after running pipeline
- `ContamASVsPrev05ToKeep.txt`: filtered results from  `DecontamResults_prev05.csv` that are ASVs to keep
- `ContamASVsPrev05ToRemove.txt`: filtered results from  `DecontamResults_prev05.csv` that are ASVs to remove (they are true contaminants)
