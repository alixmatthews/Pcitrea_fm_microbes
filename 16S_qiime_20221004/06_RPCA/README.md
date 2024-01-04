# RPCA

- Needed to create a new qiime2 env because something broke on mine with this gemelli stuff....
- Only do this first part once (if needed)

```
pwd: /home/lse305/miniconda3/envs
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2.v2 --file qiime2-2022.2-py38-linux-conda.yml
rm qiime2-2022.2-py38-linux-conda.yml
conda activate qiime2-2022.2.v2
qiime --version # to check 2022.2

pip install gemelli
pip install qurro # this messes up the diversity stats because it uninstalls pandas 1.2.5 and installs pandas-0.25.3

# so this means qiime2-2022.2.v2 will be the one I can use with qurro only

# so also install an original version that you can activate for the diversity stats on own machine
pwd: /home/lse305/miniconda3/envs
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2.original --file qiime2-2022.2-py38-linux-conda.yml
rm qiime2-2022.2-py38-linux-conda.yml
conda activate qiime2-2022.2.original
qiime --version # to check 2022.2
```

- Now can move onto the analyses

```
# on my own machine for the time being because cannot install gemelli or qurro on AHPCC
# AHPCC did install these, but they are for two different qiime versions...
conda activate qiime2-2022.2.v2
qiime info

# shows versions of plugins I'll be using in this section
# gemelli: 0.0.8
# qurro: 0.6.0

cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/04_SeparateSamples
```

## 1. Gemelli



### Use 1450 as the p-min-sample-count
- Try with a few tweaks (recommendations from tutorial and manual)
- Doing 1450 instead of 1451 to try to capture the same # of samples as the other diversity stats (in qiime2)

```
qiime gemelli auto-rpca \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --p-min-sample-count 1450 \
  --o-biplot ../06_RPCA/auto-rpca_biosamples_ordination_1450.qza \
  --o-distance-matrix ../06_RPCA/auto-rpca_biosamples_distance_1450.qza
  
qiime emperor biplot \
  --i-biplot ../06_RPCA/auto-rpca_biosamples_ordination_1450.qza \
  --m-sample-metadata-file ../16S_sample_metadata.tsv \
  --m-feature-metadata-file ../02_Run1Run2_combo/dada2_output/taxonomy/taxonomy_16S_grouped-silva.qza \
  --o-visualization ../06_RPCA/auto-rpca_biosamples_biplot_1450.qzv
```



- Export results for plotting

```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/06_RPCA

qiime tools export \
  --input-path auto-rpca_biosamples_ordination_1450.qza \
  --output-path exported
  
cd exported
mv ordination.txt auto-rpca_biosamples_ordination_1450.txt
cd ..
```



## 2. Qurro loading plots


```
cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA
```



### Threshold of 1450

```
qiime qurro loading-plot \
  --i-ranks auto-rpca_biosamples_ordination_1450.qza \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --m-sample-metadata-file ../16S_sample_metadata.tsv \
  --m-feature-metadata-file ../02_Run1Run2_combo/dada2_output/taxonomy/taxonomy_16S_grouped-silva.qza \
  --o-visualization qurro_plot_1450.qzv
```




## 3. Beta group significance

- Switch over to different qiime environment now; may be best to just do this in a new window


```
conda activate qiime2-2022.2.original

cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA
```



### 3cii. Threshold of 1450

```
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_1450.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permanova \
  --o-visualization category_broad_beta_significance_permanova_1450.qzv
  
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_1450.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization category_broad_beta_significance_permdisp_1450.qzv
```



---




# RPCA analysis on bio samples with >500 feature count

```
# on my own machine
conda activate qiime2-2022.2.v2
```


```
# make a new dir for results
cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA
mkdir bio_samples_500plus
```


```
# change to correct dir
cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/04_SeparateSamples
```

## Run gemmeli RPCA on 500plus bio samples
- I'm setting the p-min-sample-count to 500, but it is not really necessary since this table only includes samples with p-min-sample-counts of 500+ (putting as 499 because the logical vector in gemelli is "remove that value and BELOW" instead of "keep that value and ABOVE")

```
qiime gemelli auto-rpca \
  --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --p-min-sample-count 499 \
  --o-biplot ../06_RPCA/bio_samples_500plus/auto-rpca_biosamples_ordination_500.qza \
  --o-distance-matrix ../06_RPCA/bio_samples_500plus/auto-rpca_biosamples_distance_500.qza
  
  
qiime emperor biplot \
  --i-biplot ../06_RPCA/bio_samples_500plus/auto-rpca_biosamples_ordination_500.qza \
  --m-sample-metadata-file ../16S_sample_metadata.tsv \
  --m-feature-metadata-file ../02_Run1Run2_combo/dada2_output/taxonomy/taxonomy_16S_grouped-silva.qza \
  --o-visualization ../06_RPCA/bio_samples_500plus/auto-rpca_biosamples_biplot_500.qzv
```


```
# export results for plotting

cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA/bio_samples_500plus

qiime tools export \
  --input-path auto-rpca_biosamples_ordination_500.qza \
  --output-path exported
  
cd exported
mv ordination.txt auto-rpca_biosamples_ordination_500.txt
cd ..
```



## 4.2. Qurro loading plots

```
cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA/bio_samples_500plus

qiime qurro loading-plot \
  --i-ranks auto-rpca_biosamples_ordination_500.qza \
  --i-table ../../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
  --m-sample-metadata-file ../../16S_sample_metadata.tsv \
  --m-feature-metadata-file ../../02_Run1Run2_combo/dada2_output/taxonomy/taxonomy_16S_grouped-silva.qza \
  --o-visualization qurro_plot_500.qzv
```


## 4.3. Beta group significance

- It's best to start a new terminal and load the other qiime2 environment

```
conda activate qiime2-2022.2.original

cd /home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/06_RPCA/bio_samples_500plus

qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_500.qza \
  --m-metadata-file ../../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permanova \
  --o-visualization category_broad_beta_significance_permanova_500.qzv
  
qiime diversity beta-group-significance \
  --i-distance-matrix auto-rpca_biosamples_distance_500.qza \
  --m-metadata-file ../../16S_sample_metadata.tsv \
  --m-metadata-column category_broad \
  --p-method permdisp \
  --o-visualization category_broad_beta_significance_permdisp_500.qzv
```






