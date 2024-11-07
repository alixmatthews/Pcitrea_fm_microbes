# Separate samples

- Separating biological from control samples
- Removing biological samples with <500 feature count
- This is a pretty straightforward step, but want it to be separate from the alpha and beta diversity steps


#### Files here include:
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.tsv`: a table showing filtered, decontaminated, biological samples (columns) that have 500+ seqs associated with them, their read counts (values in matrix, count data), and associated ASV (rows)
- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus_transposed.tsv`: the same as `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.tsv`, but tranposed x-y... a table showing filtered, decontaminated, biological samples (rows) that have 500+ seqs associated with them, their read counts (values in matrix, count data), and associated ASV (columns)


---

## 1. Separate biological from control samples

---

### 1.1. Pull out biological samples in the table
- In other words, exclude the control samples

- Let's get a new directory made first
```
# pwd: /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

mkdir 05_SeparateSamples
```


- Now separate those samples!
```
qiime feature-table filter-samples \
  --i-table 04_decontam/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam.qza \
  --m-metadata-file ITS_forward_metadata.tsv \
  --p-where "bio_or_control='control'" \
  --p-exclude-ids true \
  --o-filtered-table 05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza
```


   
- Take a look at this table and see if it did what we expected.

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples

## Summarize the filtered biological samples data table
qiime feature-table summarize \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --o-visualization ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qzv \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv
```

- It sure did!










### 1.2. Pull out control samples in the table
- In other words, exclude the biological samples

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

qiime feature-table filter-samples \
  --i-table 04_decontam/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam.qza \
  --m-metadata-file ITS_forward_metadata.tsv \
  --p-where "bio_or_control='biological'" \
  --p-exclude-ids true \
  --o-filtered-table 05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-control_samples.qza
```
   
- Take a look at this table and see if it did what we expected.

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples

## Summarize the filtered biological samples data table
qiime feature-table summarize \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-control_samples.qza \
  --o-visualization ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-control_samples.qzv \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv
```

- Sure 'nuff!


---


### 1.3. Prep ASV table (ASV counts by sample) for downstream analyses

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples

qiime tools export --input-path ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza --output-path table

# that results a biom table. Then convert it into tsv

biom convert --to-tsv -i table/feature-table.biom -o table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.tsv

nano table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.tsv

# remove `# Constructed from biom file`
# change ``# OTU ID` to `index`

awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS $i: $i) } END{ for (i in a) print a[i] }' table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.tsv > table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_transposed.tsv

# for some strange reason, this puts the header in a weird place, but the rest is fine.
# manually replace the header where it's supposed to go, copy and paste over to BBEdit, and reimport back into the server
```











---
## 2. Separate biological samples that have <500 features
---

### 2.1. Remove bio samples that have <500 feature count

- Got the sample IDs from viewing the previous table

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples
```

```
qiime feature-table filter-samples \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --p-where "\"sample-id\" IN ('PM13', 'PM03', 'PM06', 'PF03', 'PM09', 'PM08', 'PM04', 'PF25', 'PF10', 'PF13', 'PF30')" \
  --p-exclude-ids \
  --o-filtered-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza
```
   
- Take a look at this table and see if it did what we expected.

```
## Summarize the filtered biological samples with >500 features data table
qiime feature-table summarize \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --o-visualization ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qzv \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv
```

- It sure did!


### 2.2. Pull out biological samples >500 feature count in the rep-seqs

```
# pwd: /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples


# pull 'em out
qiime feature-table filter-seqs \
  --i-data ../02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qza \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --o-filtered-data ITS_forward_afterQtrim-pmin1-dada2_rep_seq-fungi-decontam-bio_samples_500plus.qza
  
# take a look
qiime feature-table tabulate-seqs \
  --i-data ITS_forward_afterQtrim-pmin1-dada2_rep_seq-fungi-decontam-bio_samples_500plus.qza \
  --o-visualization ITS_forward_afterQtrim-pmin1-dada2_rep_seq-fungi-decontam-bio_samples_500plus.qzv
```


### 2.3. Prep ASV table (ASV counts by sample) for downstream analyses

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/05_SeparateSamples

qiime tools export --input-path ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza --output-path table

# that results a biom table. Then convert it into tsv

biom convert --to-tsv -i table/feature-table.biom -o table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.tsv

nano table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.tsv

# remove `# Constructed from biom file`
# change ``# OTU ID` to `index`

awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS $i: $i) } END{ for (i in a) print a[i] }' table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.tsv > table/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus_transposed.tsv

# for some strange reason, this puts the header in a weird place, but the rest is fine.
# manually replace the header where it's supposed to go, copy and paste over to BBEdit, and reimport back into the server
```


---


### Now move on to Alpha and beta diversity (rarefied data and `500plus` samples) and RPCA (unrarefied data and `500plus` samples)!

