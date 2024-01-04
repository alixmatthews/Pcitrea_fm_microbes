# Alpha and beta diversity statistics
#### Biological samples only

## 1. Preparation
- Get some directories set up
- Phylogeny to use will be the 'global phylogeny' that was estimated once the runs were combined (after dada2 and before mitochondrial, archael, etc. removal)
  - anything that is not included in the input table in which the phylogeny is run against will be dropped, so do not need to estimate another phylogeny using the dropped samples. Just use the original one!

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004
mkdir 05_AlphaBetaDiversity
```

### Before moving forward...
- Here is where we'll want to play around with the --p-sampling-depth to see where we can retain the most samples without losing resolution of any patterns we uncover at a higher depth.
- Need to inspect the `16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qzv` sample detail to find the value before the biggest drop off in feature counts.
- Upon inspection, a few potential drop offs. Will keep 1451 as a potential value to 'rarefy' by, but we are going to just keep samples with >500 feature counts.



## Use 1451 as --p-sampling-depth


### Calculate core diversity stats

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/05_AlphaBetaDiversity

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ../02_Run1Run2_combo/dada2_output/phylogeny/16S_FR_Run1Run2_afterQtrim-rooted_tree.qza \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --p-sampling-depth 1451 \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --output-dir merged-core-metrics-results-1451
```





- Export these PCOA results

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/05_AlphaBetaDiversity

qiime tools export \
  --input-path merged-core-metrics-results-1451/bray_curtis_pcoa_results.qza \
  --output-path merged-core-metrics-results-1451/exported
  
cd merged-core-metrics-results-1451/exported
mv ordination.txt bray_curtis_pcoa_results.txt
cd ../..
  
qiime tools export \
  --input-path merged-core-metrics-results-1451/jaccard_pcoa_results.qza \
  --output-path merged-core-metrics-results-1451/exported
  
cd merged-core-metrics-results-1451/exported
mv ordination.txt jaccard_pcoa_results.txt
cd ../..

qiime tools export \
  --input-path merged-core-metrics-results-1451/unweighted_unifrac_pcoa_results.qza \
  --output-path merged-core-metrics-results-1451/exported
  
cd merged-core-metrics-results-1451/exported
mv ordination.txt unweighted_unifrac_pcoa_results.txt
cd ../..
  
qiime tools export \
  --input-path merged-core-metrics-results-1451/weighted_unifrac_pcoa_results.qza \
  --output-path merged-core-metrics-results-1451/exported
  
cd merged-core-metrics-results-1451/exported
mv ordination.txt weighted_unifrac_pcoa_results.txt
```



### Alpha diversity
- Alpha diversity means community variation within a community (sample)

#### Shannon Index
- Compare community richness and diversity (quantitative measure that accounts for both abundance and evenness) 

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity merged-core-metrics-results-1451/shannon_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization merged-core-metrics-results-1451/shannon-group-significance.qzv
```




#### Observed Features
- Compare community richness (a qualitative measure of community richness) 

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity merged-core-metrics-results-1451/observed_features_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization merged-core-metrics-results-1451/observed_features-group-significance.qzv
```




####  Faith's phylogenetic diversity
- Compare community richness (qualitiative biodiversity measure accounting for phylogenetic relatedness between features)

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity merged-core-metrics-results-1451/faith_pd_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization merged-core-metrics-results-1451/faith-pd-group-significance.qzv
```



### Beta diversity
- Beta diversity quantifies (dis-)similarites between communities (samples)

#### Jaccard distance
- A qualitative measure of community dissimilarity

- PERMANOVA first (default)
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/jaccard_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization merged-core-metrics-results-1451/jaccard-significance-category_broad-permanova.qzv \
  --p-pairwise
```

- Permdisp
  - Checking to see if the differences within groups are smaller than between groups (i.e., dispersion within group)

```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/jaccard_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization merged-core-metrics-results-1451/jaccard-significance-category_broad-permdisp.qzv \
  --p-pairwise
```

#### Bray-Curtis distance
- A quantitative measure of community dissimilarity

- PERMANOVA first (default)
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/bray_curtis_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization merged-core-metrics-results-1451/bray_curtis-significance-category_broad-permanova.qzv \
  --p-pairwise
```

- Permdisp
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/bray_curtis_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization merged-core-metrics-results-1451/bray_curtis-significance-category_broad-permdisp.qzv \
  --p-pairwise
```


#### Unweighted unifrac
- A qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features

- PERMANOVA first (default)
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization merged-core-metrics-results-1451/unweighted-unifrac-significance-category_broad-permanova.qzv \
  --p-pairwise
```




- Perdisp
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization merged-core-metrics-results-1451/unweighted-unifrac-significance-category_broad-permdisp.qzv \
  --p-pairwise
```




#### Weighted unifrac
- A quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features

- PERMANOVA first (default)
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization merged-core-metrics-results-1451/weighted-unifrac-significance-category_broad-permanova.qzv \
  --p-pairwise
```



- Perdisp
```
qiime diversity beta-group-significance \
  --i-distance-matrix merged-core-metrics-results-1451/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization merged-core-metrics-results-1451/weighted-unifrac-significance-category_broad-permdisp.qzv \
  --p-pairwise
```





## Biological samples only that have at least 500 feature count in the table

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/05_AlphaBetaDiversity
mkdir metrics_500plus
```

### Faith's PD

```
qiime diversity alpha-phylogenetic \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --i-phylogeny ../02_Run1Run2_combo/dada2_output/phylogeny/16S_FR_Run1Run2_afterQtrim-rooted_tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity metrics_500plus/faith_pd_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/faith_pd_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization metrics_500plus/faith_pd_group_significance.qzv
```

### Shannon

```
qiime diversity alpha \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric shannon \
  --o-alpha-diversity metrics_500plus/shannon_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/shannon_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization metrics_500plus/shannon_group_significance.qzv
```

### Observed features

```
qiime diversity alpha \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric observed_features \
  --o-alpha-diversity metrics_500plus/observed_features_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/observed_features_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization metrics_500plus/observed_features_group_significance.qzv
```

### Chao1
```
qiime diversity alpha \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric chao1 \
  --o-alpha-diversity metrics_500plus/chao1_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/chao1_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization metrics_500plus/chao1_group_significance.qzv
```


### Pielou's
```
qiime diversity alpha \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric pielou_e \
  --o-alpha-diversity metrics_500plus/pielou_vector.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_500plus/pielou_vector.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization metrics_500plus/pielou_group_significance.qzv
```


###  Jaccard

```
qiime diversity beta \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric jaccard \
  --o-distance-matrix metrics_500plus/jaccard_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/jaccard_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/jaccard-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/jaccard_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
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


###  Bray-Curtis

```
qiime diversity beta \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-metric braycurtis \
  --o-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/braycurtis-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/braycurtis_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
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


### Weighted UniFrac

```
qiime diversity beta-phylogenetic \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --i-phylogeny ../02_Run1Run2_combo/dada2_output/phylogeny/16S_FR_Run1Run2_afterQtrim-rooted_tree.qza \
  --p-metric weighted_unifrac \
  --o-distance-matrix metrics_500plus/weighted_unifrac_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/weighted_unifrac-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization metrics_500plus/weighted_unifrac-significance-category_broad-permdisp.qzv \
  --p-pairwise
  
qiime diversity pcoa \
  --i-distance-matrix metrics_500plus/weighted_unifrac_distance_matrix.qza \
  --o-pcoa metrics_500plus/weighted_unifrac_pcoa_results.qza

qiime tools export \
  --input-path metrics_500plus/weighted_unifrac_pcoa_results.qza \
  --output-path metrics_500plus/exported
  
cd metrics_500plus/exported
mv ordination.txt weighted_unifrac_pcoa_results.txt
cd ../..
```


###  Unweighted UniFrac

```
qiime diversity beta-phylogenetic \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --i-phylogeny ../02_Run1Run2_combo/dada2_output/phylogeny/16S_FR_Run1Run2_afterQtrim-rooted_tree.qza \
  --p-metric unweighted_unifrac \
  --o-distance-matrix metrics_500plus/unweighted_unifrac_distance_matrix.qza

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --o-visualization metrics_500plus/unweighted_unifrac-significance-category_broad-permanova.qzv \
  --p-pairwise
  
qiime diversity beta-group-significance \
  --i-distance-matrix metrics_500plus/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization metrics_500plus/unweighted_unifrac-significance-category_broad-permdisp.qzv \
  --p-pairwise
  
qiime diversity pcoa \
  --i-distance-matrix metrics_500plus/unweighted_unifrac_distance_matrix.qza \
  --o-pcoa metrics_500plus/unweighted_unifrac_pcoa_results.qza

qiime tools export \
  --input-path metrics_500plus/unweighted_unifrac_pcoa_results.qza \
  --output-path metrics_500plus/exported
  
cd metrics_500plus/exported
mv ordination.txt unweighted_unifrac_pcoa_results.txt
cd ../..
```



