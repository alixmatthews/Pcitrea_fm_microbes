# RPCA

```
# on my own machine for the time being because cannot install gemelli or qurro on AHPCC
# AHPCC did install these, but they are for two different qiime versions...
conda activate qiime2-2022.2.v2
qiime info

# shows versions of plugins I'll be using in this section
# gemelli: 0.0.8
# qurro: 0.6.0

cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/05_SeparateSamples
```



## 1. Gemelli



### 1b. Use 994 as the p-min-sample-count
- Try with a few tweaks (recommendations from tutorial and manual)
- This is because using a depth of 995 results in 1 less sample the other diversity stats from qiime2, so the function of the cutoff for p-min-sample-count is "remove that value and BELOW" instead of "keep that value and ABOVE" (confirmed with 16S analyses this is truly the case)

```
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/05_SeparateSamples

qiime gemelli auto-rpca \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --p-min-sample-count 994 \
  --o-biplot ../08_RPCA/auto-rpca_biosamples_ordination_994.qza \
  --o-distance-matrix ../08_RPCA/auto-rpca_biosamples_distance_994.qza
  
qiime emperor biplot \
  --i-biplot ../08_RPCA/auto-rpca_biosamples_ordination_994.qza \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv \
  --m-feature-metadata-file ../03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --o-visualization ../08_RPCA/auto-rpca_biosamples_biplot_994.qzv
```




- Export results for plotting

```
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA

qiime tools export \
  --input-path auto-rpca_biosamples_ordination_994.qza \
  --output-path exported
  
cd exported
mv ordination.txt auto-rpca_biosamples_ordination_994.txt
cd ..
```












## 2. Qurro loading plots


```
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA
```


### 2b. Threshold of 994

```
qiime qurro loading-plot \
  --i-ranks auto-rpca_biosamples_ordination_994.qza \
  --i-table ../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv \
  --m-feature-metadata-file ../03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --o-visualization qurro_plot_994.qzv
```









## 3. Beta group significance

- Switch over to different qiime environment now; may be best to just do this in a new window


```
conda activate qiime2-2022.2.original

cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA
```




### 3b. Threshold of 994

```
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_994.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permanova \
  --o-visualization category_broad_beta_significance_permanova_994.qzv
  
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_994.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization category_broad_beta_significance_permdisp_994.qzv
```






---






# 4. RPCA analysis on bio samples with >500 feature count

```
# on my own machine
conda activate qiime2-2022.2.v2
```


```
# make a new dir for results
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA
mkdir bio_samples_500plus
```


```
# change to correct dir
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/05_SeparateSamples
```

## 4.1. Run gemmeli RPCA on 500plus bio samples
- I'm setting the p-min-sample-count to 500, but it is not really necessary since this table only includes samples with p-min-sample-counts of 500+ (putting as 499 because the logical vector in gemelli is "remove that value and BELOW" instead of "keep that value and ABOVE")

```
qiime gemelli auto-rpca \
  --i-table ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --p-min-sample-count 499 \
  --o-biplot ../08_RPCA/bio_samples_500plus/auto-rpca_biosamples_ordination_500.qza \
  --o-distance-matrix ../08_RPCA/bio_samples_500plus/auto-rpca_biosamples_distance_500.qza
  
  
qiime emperor biplot \
  --i-biplot ../08_RPCA/bio_samples_500plus/auto-rpca_biosamples_ordination_500.qza \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv \
  --m-feature-metadata-file ../03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --o-visualization ../08_RPCA/bio_samples_500plus/auto-rpca_biosamples_biplot_500.qzv
```



```
# export results for plotting

cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA/bio_samples_500plus

qiime tools export \
  --input-path auto-rpca_biosamples_ordination_500.qza \
  --output-path exported
  
cd exported
mv ordination.txt auto-rpca_biosamples_ordination_500.txt
cd ..
```



## 4.2. Qurro loading plots

```
cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA/bio_samples_500plus

qiime qurro loading-plot \
  --i-ranks auto-rpca_biosamples_ordination_500.qza \
  --i-table ../../05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --m-sample-metadata-file ../../ITS_forward_metadata.tsv \
  --m-feature-metadata-file ../../03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --o-visualization qurro_plot_500.qzv
```


## 4.3. Beta group significance

- It's best to start a new terminal and load the other qiime2 environment

```
conda activate qiime2-2022.2.original

cd /home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/08_RPCA/bio_samples_500plus

qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_500.qza \
  --m-metadata-file ../../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permanova \
  --o-visualization category_broad_beta_significance_permanova_500.qzv
  
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_500.qza \
  --m-metadata-file ../../ITS_forward_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization category_broad_beta_significance_permdisp_500.qzv
```







