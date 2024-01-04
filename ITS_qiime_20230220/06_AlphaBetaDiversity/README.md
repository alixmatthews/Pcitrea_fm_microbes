# Alpha and beta diversity statistics
#### Biological samples only

## 1. Preparation
- Get some directories set up
- Phylogeny to use will be the 'global phylogeny' that was estimated (after dada2, before UNITE/mitochondrial, archael, etc. removal)
  - anything that is not included in the input table in which the phylogeny is run against will be dropped, so do not need to estimate another phylogeny using the dropped samples. Just use the original one!
  - **Important to note:** ITS phylogeny is NOT phylogenetically informative due to high genetic diversity among unrelated taxa. So cannot use/rely on phylogenetic beta metrics.

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
mkdir 06_AlphaBetaDiversity
```


### Before moving forward...
- Here is where we'll want to play around with the --p-sampling-depth to see where we can retain the most samples without losing resolution of any patterns we uncover at a higher depth.
- Need to inspect the `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qzv` sample detail to do so.
- Upon inspection, we don't have a lot of data to work with... 
  - If we do 995, we have 13 mite samples and 24 feather samples
  - We will also consider all samples that have >500 seq counts



## 2. Use 995 as --p-sampling-depth


### 2a. Calculate core diversity stats

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 02_dada2/phylogeny/ITS_forward_afterQtrim-pmin1-dada2_repseqs-rooted_tree.qza \
  --i-table 05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --p-sampling-depth 995 \
  --m-metadata-file ITS_forward_metadata.tsv \
  --output-dir 06_AlphaBetaDiversity/forward-core-metrics-results-995
```




- Export these PCOA results

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/06_AlphaBetaDiversity

qiime tools export \
  --input-path forward-core-metrics-results-995/bray_curtis_pcoa_results.qza \
  --output-path forward-core-metrics-results-995/exported
  
cd forward-core-metrics-results-995/exported
mv ordination.txt bray_curtis_pcoa_results.txt
cd ../..
  
qiime tools export \
  --input-path forward-core-metrics-results-995/jaccard_pcoa_results.qza \
  --output-path forward-core-metrics-results-995/exported
  
cd forward-core-metrics-results-995/exported
mv ordination.txt jaccard_pcoa_results.txt
cd ../..

# skipping UniFrac because they are phylogeneticlally-informed and not meaningful for ITS with how I have done this...
```



### 2b. Alpha diversity
- Alpha diversity means community variation within a community (sample)

#### Shannon Index
- Compare community richness and diversity (quantitative measure that accounts for both abundance and evenness) 

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/06_AlphaBetaDiversity

qiime diversity alpha-group-significance \
  --i-alpha-diversity forward-core-metrics-results-995/shannon_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization forward-core-metrics-results-995/shannon-group-significance.qzv
```



#### Observed Features
- Compare community richness (a qualitative measure of community richness) 

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity forward-core-metrics-results-995/observed_features_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization forward-core-metrics-results-995/observed_features-group-significance.qzv
```



####  Faith's phylogenetic diversity
- Compare community richness (qualitiative biodiversity measure accounting for phylogenetic relatedness between features)

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity forward-core-metrics-results-995/faith_pd_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization forward-core-metrics-results-995/faith-pd-group-significance.qzv
```








### 2c. Beta diversity
- Beta diversity quantifies (dis-)similarites between communities (samples)

#### Jaccard distance
- A qualitative measure of community dissimilarity

- PERMANOVA first (default)
```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/06_AlphaBetaDiversity

qiime diversity beta-group-significance \
  --i-distance-matrix forward-core-metrics-results-995/jaccard_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization forward-core-metrics-results-995/jaccard-significance-category_broad-permanova.qzv \
  --p-pairwise
```








- Permdisp
  - Checking to see if the differences within groups are smaller than between groups (i.e., dispersion within group)

```
qiime diversity beta-group-significance \
  --i-distance-matrix forward-core-metrics-results-995/jaccard_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization forward-core-metrics-results-995/jaccard-significance-category_broad-permdisp.qzv \
  --p-pairwise
```






#### Bray-Curtis distance
- A quantitative measure of community dissimilarity

- PERMANOVA first (default)
```
qiime diversity beta-group-significance \
  --i-distance-matrix forward-core-metrics-results-995/bray_curtis_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization forward-core-metrics-results-995/bray_curtis-significance-category_broad-permanova.qzv \
  --p-pairwise
```




- Permdisp
```
qiime diversity beta-group-significance \
  --i-distance-matrix forward-core-metrics-results-995/bray_curtis_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization forward-core-metrics-results-995/bray_curtis-significance-category_broad-permdisp.qzv \
  --p-pairwise
```





---
---
---










## 3. Biological samples only that have at least 500 feature count in the table

- Cannot do Faith's PD or Weighted UniFrac or Unweighted UniFrac because these are phylogenetically informed indices and these ITS data are not phylogetically informative.

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220/06_AlphaBetaDiversity
mkdir metrics_500plus
```

### 3.1. Shannon

```
qiime diversity alpha \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric shannon \
  --o-alpha-diversity metrics_500plus/shannon_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/shannon_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization metrics_500plus/shannon_group_significance.qzv
```

### 3.2. Observed features

```
qiime diversity alpha \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric observed_features \
  --o-alpha-diversity metrics_500plus/observed_features_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/observed_features_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization metrics_500plus/observed_features_group_significance.qzv
```

### 3.3 Chao1
```
qiime diversity alpha \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric chao1 \
  --o-alpha-diversity metrics_500plus/chao1_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/chao1_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization metrics_500plus/chao1_group_significance.qzv
```


### 3.4 Pielou's
```
qiime diversity alpha \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric pielou_e \
  --o-alpha-diversity metrics_500plus/pielou_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/pielou_vector.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization metrics_500plus/pielou_group_significance.qzv
```

### 3.5 Jaccard

```
qiime diversity beta \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric jaccard \
  --o-distance-matrix metrics_500plus/jaccard_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/jaccard_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/jaccard-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/jaccard_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization metrics_500plus/jaccard-significance-category_broad-permdisp.qzv \
  --p-pairwise
  
qiime diversity pcoa \
  --i-distance-matrix metrics_500plus/jaccard_distance_matrix.qza \
  --o-pcoa metrics_500plus/jaccard_pcoa_results.qza

qiime tools export \
  --input-path metrics_500plus/jaccard_pcoa_results.qza \
  --output-path metrics_500plus/exported
  
cd metrics_500plus/exported
mv ordination.txt jaccard_pcoa_results.txt
cd ../..
```


### 3.6. Bray-Curtis

```
qiime diversity beta \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-metric braycurtis \
  --o-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/braycurtis-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization metrics_500plus/braycurtis-significance-category_broad-permdisp.qzv \
  --p-pairwise
  
qiime diversity pcoa \
  --i-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza \
  --o-pcoa metrics_500plus/braycurtis_pcoa_results.qza

qiime tools export \
  --input-path metrics_500plus/braycurtis_pcoa_results.qza \
  --output-path metrics_500plus/exported
  
cd metrics_500plus/exported
mv ordination.txt braycurtis_pcoa_results.txt
cd ../..
```






