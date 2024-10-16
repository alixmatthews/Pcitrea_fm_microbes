Here we are running a decontam pipeline in R.

#### Data file in this directory and used for this step:
- `decontam_20221005.R`: decontam R pipeline
- `16S_filtered-all_samples-feature_table-with-taxonomy.biom`: biom file with all samples, features, and taxonomy
- `DecontamResults_prev05.csv`: results after running pipeline
- `ContamASVsPrev05ToKeep.txt`: filtered results from  `DecontamResults_prev05.csv` that are ASVs to keep
- `ContamASVsPrev05ToRemove.txt`: filtered results from  `DecontamResults_prev05.csv` that are ASVs to remove (they are true contaminants)


## Specific steps to get to that point

- First need to prepare input files for decontam
- I need a feature table converted to a .biom file and the metadata file for import into R for running decontam
- This is a multi-step process

---


### 1. First export feature table as a .biom file

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output/only_bacteria

qiime tools export \
  --input-path 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered.qza \
  --output-path exported-feature-table
  
cd exported-feature-table

mv feature-table.biom 16S_filtered-all_samples-feature_table.biom
```

### 2. Then export taxonomy info

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output

qiime tools export \
  --input-path taxonomy/taxonomy_16S_grouped-silva.qza \
  --output-path only_bacteria/exported-feature-table

cd only_bacteria/exported-feature-table

mv taxonomy.tsv 16S_filtered-all_samples-feature_taxonomy.tsv
```


### 3. Change the headers of this .tsv file

```
nano 16S_filtered-all_samples-feature_taxonomy.tsv
```

- Edit the column names and change
`Feature ID` to `#OTUID`
`Taxon` to `taxonomy`
`Confidence` to `confidence`

- ctrl + X to save, yes, save as same name


### 4. Combine taxonomy and metadata and remake a new biom file to import into R

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output/only_bacteria


biom add-metadata \
  --input-fp exported-feature-table/16S_filtered-all_samples-feature_table.biom \
  --output-fp exported-feature-table/16S_filtered-all_samples-feature_table-with-taxonomy.biom \
  --observation-metadata-fp exported-feature-table/16S_filtered-all_samples-feature_taxonomy.tsv \
  --sc-separated taxonomy
```


### 5. Now you're ready to run through the `decontam_20221005.R` file with `16S_filtered-all_samples-feature_table-with-taxonomy.biom`!

- Tried prevalence thresholds of 0.1 and 0.5, but I think 0.5 is what we should stick with.


### 6. Check out the results file: `DecontamResults_prev05.csv` 
- Now you can make a contaminants file from this and then remove these in qiime2 
- Need to filter the .csv by contaminant/TRUE
- Copy these TRUE contaminants (first column) over to a text file
- Insert a header named "#OTUID"
- Save this file as `ContamASVsPrev05ToRemove.txt`
- Ended up also doing the inverse as a sanity check (contaminant/FALSE) = `ContamASVsPrev05ToKeep.txt`


### 7. Move back into qiime and filter these contaminant features out!

- First, it'd be a good idea to make a directory with decontam results, so copy these over to the new directory (and whatever other files you think might be nice to have over there)

#### To remove true contaminants
```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004

qiime feature-table filter-features \
   --i-table 02_Run1Run2_combo/dada2_output/only_bacteria/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered.qza \
   --m-metadata-file 03_decontam/ContamASVsPrev05ToRemove.txt \
   --p-exclude-ids \
   --o-filtered-table 03_decontam/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam.qza

# Summarize the decontam'ed data
qiime feature-table summarize \
  --i-table 03_decontam/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam.qza \
  --o-visualization 03_decontam/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam.qzv \
  --m-sample-metadata-file 16S_sample_metadata.tsv
```



### Now move on to separating the bio from control samples for further analysis! (`04_SeparateSamples`)!

