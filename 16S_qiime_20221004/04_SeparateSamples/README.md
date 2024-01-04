# Separate samples

- Separating biological from control samples
- Removing biological samples with <500 feature count
- This is a pretty straightforward step, but want it to be separate from the alpha and beta diversity steps


#### Files here include:
- `16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.fasta`: fasta file containing features (aka, ASVs; MD5 identifier and representative sequence) for filtered, decontaminated, biological samples that have 500+ seqs associated with them
- `16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv`: a table showing filtered, decontaminated, biological samples that have 500+ seqs associated with them, their read counts (values), and associated ASV
- `16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus_transposed.tsv`: the same as `16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv`, but tranposed x-y... a table showing filtered, decontaminated, biological samples that have 500+ seqs associated with them, their read counts (values), and associated ASV


---

## 1. Separate biological from control samples

---
### 1.1. Pull out biological samples in the table
- In other words, exclude the control samples

- Let's get a new directory made first
```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/03_decontam

mkdir ../04_SeparateSamples
```


- Now separate those samples!
```
qiime feature-table filter-samples \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --p-where "bio_or_control='control'" \
  --p-exclude-ids true \
  --o-filtered-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza
```
   
- Take a look at this table and see if it did what we expected.

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

## Summarize the filtered biological samples data table
qiime feature-table summarize \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qzv \
  --m-sample-metadata-file ../16S_sample_metadata.tsv
```

- It sure did!






### 1.2. Pull out control samples in the table
- In other words, exclude the biological samples

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/03_decontam

qiime feature-table filter-samples \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --p-where "bio_or_control='biological'" \
  --p-exclude-ids true \
  --o-filtered-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-control_samples.qza
```
   
- Take a look at this table and see if it did what we expected.

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

## Summarize the filtered biological samples data table
qiime feature-table summarize \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-control_samples.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-control_samples.qzv \
  --m-sample-metadata-file ../16S_sample_metadata.tsv
```

- Sure 'nuff!


### 1.3. Pull out biological samples in the rep-seqs

```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

# pull 'em out
qiime feature-table filter-seqs \
  --i-data ../02_Run1Run2_combo/dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --o-filtered-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples.qza
  
# take a look
qiime feature-table tabulate-seqs \
  --i-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples.qzv
```



### 1.4. Pull out control samples in the rep-seqs

```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

# pull 'em out
qiime feature-table filter-seqs \
  --i-data ../02_Run1Run2_combo/dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-control_samples.qza \
  --o-filtered-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-control_samples.qza
  
# take a look
qiime feature-table tabulate-seqs \
  --i-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-control_samples.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-control_samples.qzv
```




### 1.5. Prep ASV table (ASV counts by sample) for downstream analyses

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

qiime tools export --input-path 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza --output-path table

# that results a biom table. Then convert it into tsv

biom convert --to-tsv -i table/feature-table.biom -o table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.tsv

nano table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.tsv

# remove `# Constructed from biom file`
# change ``# OTU ID` to `index`

awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS $i: $i) } END{ for (i in a) print a[i] }' table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.tsv > table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_transposed.tsv

# for some strange reason, this puts the header in a weird place, but the rest is fine.
# manually replace the header where it's supposed to go, copy and paste over to BBEdit, and reimport back into the server
```


---
## 2. Separate biological samples that have <500 features
---

### 2.1. Remove bio samples that have <500 feature count

- Got the sample IDs from viewing the previous table

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples
```

```
qiime feature-table filter-samples \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --p-where "\"sample-id\" IN ('PF05', 'PM24', 'PM29', 'PF14', 'PF15', 'PF04', 'PM19', 'PF19', 'PM06', 'PF28', 'PM08', 'PF25', 'PM05', 'PF20', 'PF08')" \
  --p-exclude-ids \
  --o-filtered-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza
```
   
- Take a look at this table and see if it did what we expected.

```
## Summarize the filtered biological samples with >500 features data table
qiime feature-table summarize \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qzv \
  --m-sample-metadata-file ../16S_sample_metadata.tsv
```

- It sure did!


### 2.2. Pull out biological samples >500 feature count in the rep-seqs

```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

# pull 'em out
qiime feature-table filter-seqs \
  --i-data ../02_Run1Run2_combo/dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --o-filtered-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.qza
  
# take a look
qiime feature-table tabulate-seqs \
  --i-data 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.qza \
  --o-visualization 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.qzv
```


### 2.3. Prep ASV table (ASV counts by sample) for downstream analyses

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

qiime tools export --input-path 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza --output-path table

# that results a biom table. Then convert it into tsv

biom convert --to-tsv -i table/feature-table.biom -o table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv

nano table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv

# remove `# Constructed from biom file`
# change ``# OTU ID` to `index`

awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS $i: $i) } END{ for (i in a) print a[i] }' table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv > table/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus_transposed.tsv

# for some strange reason, this puts the header in a weird place, but the rest is fine.
# manually replace the header where it's supposed to go, copy and paste over to BBEdit, and reimport back into the server
```



---


### Now move on to `05_AlphaBetaDiversity` (rarefied data and `500plus` samples) and `06_RPCA` (unrarefied data and `500plus` samples)!
