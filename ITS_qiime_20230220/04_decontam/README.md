Here we are running a decontam pipeline in R.

#### Data file in this directory and used for this step:
- `decontam_20230220.R`: decontam R pipeline
- `ITS_fungi-all_samples-feature_table-with-taxonomy.biom`: biom file with all samples, features, and taxonomy
- `DecontamResults_prev05.csv`: results after running pipeline (see description of table at end of this README file)
- `ContamASVsPrev05ToRemove.txt`: filtered results from  `DecontamResults_prev05.csv` that are ASVs to remove (they are true contaminants), list


## Specific steps to get to that point

- First need to prepare input files for decontam
- I need a feature table converted to a .biom file and the metadata file for import into R for running decontam
- This is a multi-step process

---

### 1. First export feature table as a .biom file

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/03_taxonomy/only_fungi

qiime tools export \
  --input-path ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qza \
  --output-path exported-feature-table
  
cd exported-feature-table

mv feature-table.biom ITS_fungi-all_samples-feature_table.biom
```

---

### 2. Then export taxonomy info

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/03_taxonomy

qiime tools export \
  --input-path taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --output-path only_fungi/exported-feature-table

cd only_fungi/exported-feature-table

mv taxonomy.tsv ITS_fungi-all_samples-feature_taxonomy.tsv
```

---


### 3. Change the headers of this .tsv file

```
nano ITS_fungi-all_samples-feature_taxonomy.tsv
```

- Edit the column names and change
`Feature ID` to `#OTUID`
`Taxon` to `taxonomy`
`Confidence` to `confidence`

- ctrl + X to save, yes, save as same name



---

### 4. Combine taxonomy and metadata and remake a new biom file to import into R

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/03_taxonomy/only_fungi

biom add-metadata \
  --input-fp exported-feature-table/ITS_fungi-all_samples-feature_table.biom \
  --output-fp exported-feature-table/ITS_fungi-all_samples-feature_table-with-taxonomy.biom \
  --observation-metadata-fp exported-feature-table/ITS_fungi-all_samples-feature_taxonomy.tsv \
  --sc-separated taxonomy
```

---

### 5. Now you're ready to run through the `decontam_20230220.R` file with `ITS_fungi-all_samples-feature_table-with-taxonomy.biom`!

- Tried prevalence thresholds of 0.1 and 0.5, but I think 0.5 is what we should stick with.




---

### 6. Check out the results file: `DecontamResults_prev05.csv` 
- Now you can make a contaminants file from this and then remove these in qiime2 
- Need to filter the .csv by contaminant/TRUE
- Copy these TRUE contaminants (first column) over to a text file
- Insert a header named "#OTUID"
- Save this file as `ContamASVsPrev05ToRemove.txt`


---

### 7. Move back into qiime and filter these contaminant features out!

- First, it'd be a good idea to add decontam results to decontam directory, so make a new dir and copy these results over to the decontam directory (and whatever other files you think might be nice to have over there)

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
mkdir 04_decontam
```

#### To remove true contaminants
```
# pwd: /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

qiime feature-table filter-features \
   --i-table 03_taxonomy/only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qza \
   --m-metadata-file 04_decontam/ContamASVsPrev05ToRemove.txt \
   --p-exclude-ids \
   --o-filtered-table 04_decontam/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam.qza

# Summarize the decontam'ed data
qiime feature-table summarize \
  --i-table 04_decontam/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam.qza \
  --o-visualization 04_decontam/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam.qzv \
  --m-sample-metadata-file ITS_forward_metadata.tsv
```






---

### Now move on to separating the bio from control samples for further analysis! (`05_SeparateSamples`)


---

### Variable descriptors
`DecontamResults_prev05.csv` - first columm is the ASV/OTU ID number


| Variable Name | Description                                                          | Units | Notes                |
|---------------|----------------------------------------------------------------------|-------|-------------------------------------------|
| `pa.pos`      | Proportion of abundance of each ASV in positive 'biological' samples                       | Count|         |
| `pa.neg`      | Proportion of abundance of each ASV in negative 'control' samples               | Count|         |
| `contaminant` | Classification of whether the sequence is a likely contaminant       | Boolean | `TRUE` for contaminant, `FALSE` otherwise |


