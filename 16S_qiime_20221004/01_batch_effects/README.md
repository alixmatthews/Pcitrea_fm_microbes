# 16S data analysis
# First checking for batch effects


#### Data files used for this step:
`16S_FR_Run1_manifest.tsv`: manifest file containing metadata and filepaths for Run1 only

`16S_FR_Run2_manifest.tsv`: manifest file containing metadata and filepaths for Run2 only

#### Start an interactive srun on the AHPCC

```
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --partition pcon06 --time 6:00:00 --pty /bin/bash
```

#### Load modules and check things

```
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate qiime2-2022.2

# look at the versions of all the plugins
qiime info 

# cd to working directory with all files
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/01_batch_effects
```

## Import and summarize each MiniSeq run individually
```
## import Run1
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path 16S_FR_Run1_manifest.tsv \
  --output-path 16S_FR_Run1_beforeQtrim-demux.qza

## summarize Run 1
qiime demux summarize \
--i-data 16S_FR_Run1_beforeQtrim-demux.qza \
--o-visualization 16S_FR_Run1_beforeQtrim-demux.qzv
```

```
## import Run2
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path 16S_FR_Run2_manifest.tsv \
  --output-path 16S_FR_Run2_beforeQtrim-demux.qza

## summarize Run 2
qiime demux summarize \
--i-data 16S_FR_Run2_beforeQtrim-demux.qza \
--o-visualization 16S_FR_Run2_beforeQtrim-demux.qzv
```





## Trim primers and summarize trimmed data for each MiniSeq run individually
Trimming primers will also consequently trim Illumina adapters

These primer seqs were obtained from the sequencing facility. They are the 515F/806R primer set for 16S (V4 region).

This takes less than 5 minutes

```
## trim Run 1
qiime cutadapt trim-paired \
  --p-cores 16 \
  --p-front-f ^GTGCCAGCMGCCGCGGTAA \
  --p-front-r ^GGACTACHVGGGTWTCTAAT \
  --p-error-rate 0.10 \
  --p-minimum-length 100 \
  --p-discard-untrimmed \
  --i-demultiplexed-sequences 16S_FR_Run1_beforeQtrim-demux.qza \
  --o-trimmed-sequences 16S_FR_Run1_afterQtrim-demux.qza

## summarize trimmed Run 1
qiime demux summarize \
--i-data 16S_FR_Run1_afterQtrim-demux.qza \
--o-visualization 16S_FR_Run1_afterQtrim-demux.qzv
```


```
## trim Run 2
qiime cutadapt trim-paired \
  --p-cores 16 \
  --p-front-f ^GTGCCAGCMGCCGCGGTAA \
  --p-front-r ^GGACTACHVGGGTWTCTAAT \
  --p-error-rate 0.10 \
  --p-minimum-length 100 \
  --p-discard-untrimmed \
  --i-demultiplexed-sequences 16S_FR_Run2_beforeQtrim-demux.qza \
  --o-trimmed-sequences 16S_FR_Run2_afterQtrim-demux.qza

## summarize trimmed Run 2
qiime demux summarize \
--i-data 16S_FR_Run2_afterQtrim-demux.qza \
--o-visualization 16S_FR_Run2_afterQtrim-demux.qzv
```




## DADA2 on the trimmed (still paired) data for each MiniSeq run individually
This takes about 7 minutes

This step will merge reads and remove poor quality bases based on parameters I choose.

Parameters in this step were either adjusted incrementally/experimentally or chosen from the `afterQtrim-demux.qzv` files. Final parameters were chosen based on examining the results in the stats tables files that are generated in this step. I chose parameters that allowed the most reads to be retained as possible and fewest proportion of chimeras. All parameters in this step should be chosen carefully. If you are checking for batch effects, the parameters for this step should be IDENTICAL when conducting this step on each individual run.

- `--p-chimera-method`: experiment with this 
- `--p-trunc-len-f`: check `*afterQtrim-demux.qzv` to determine
- `--p-trunc-len-r`: check `*afterQtrim-demux.qzv` to determine
- `--p-min-overlap`: experiment with this
- `--p-min-fold-parent-over-abundance`: experiment with this
- `--p-trunc-q`: experiment with this 

```
## make new dir for dada2 output
mkdir dada2_output
```

```
## DADA2 Run 1 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S_FR_Run1_afterQtrim-demux.qza \
  --p-n-threads 16 \
  --p-chimera-method consensus \
  --p-hashed-feature-ids \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 130 \
  --p-min-overlap 4 \
  --p-min-fold-parent-over-abundance 4.0 \
  --p-trunc-q 15 \
  --o-table dada2_output/16S_FR_Run1_afterQtrim-dada2_table.qza \
  --o-representative-sequences dada2_output/16S_FR_Run1_afterQtrim-dada2_rep_set.qza \
  --o-denoising-stats dada2_output/16S_FR_Run1_afterQtrim-dada2_stats.qza
  
## Generate stats tables for these new runs
qiime metadata tabulate \
  --m-input-file dada2_output/16S_FR_Run1_afterQtrim-dada2_stats.qza \
  --o-visualization dada2_output/16S_FR_Run1_afterQtrim-dada2_stats.qzv
  
```




```
## DADA2 Run 2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S_FR_Run2_afterQtrim-demux.qza \
  --p-n-threads 16 \
  --p-chimera-method consensus \
  --p-hashed-feature-ids \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 130 \
  --p-min-overlap 4 \
  --p-min-fold-parent-over-abundance 4.0 \
  --p-trunc-q 15 \
  --o-table dada2_output/16S_FR_Run2_afterQtrim-dada2_table.qza \
  --o-representative-sequences dada2_output/16S_FR_Run2_afterQtrim-dada2_rep_set.qza \
  --o-denoising-stats dada2_output/16S_FR_Run2_afterQtrim-dada2_stats.qza

qiime metadata tabulate \
  --m-input-file dada2_output/16S_FR_Run2_afterQtrim-dada2_stats.qza \
  --o-visualization dada2_output/16S_FR_Run2_afterQtrim-dada2_stats.qzv
```




## Merge and summarize denoised data 
This takes less than 1 minute

Because I have multiple sequencing runs, I need to proceed with merging the table and sequences from separate DADA2 runs.

```
## Let's set up new directories for merged data
mkdir Runs_merged_output
cd Runs_merged_output
mkdir dada2_output
cd ..

# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/01_batch_effects
```


```
## Merge tables
qiime feature-table merge \
  --i-tables dada2_output/16S_FR_Run1_afterQtrim-dada2_table.qza \
  --i-tables dada2_output/16S_FR_Run2_afterQtrim-dada2_table.qza \
  --o-merged-table Runs_merged_output/dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_table.qza

## Merge rep_seq sets
qiime feature-table merge-seqs \
  --i-data dada2_output/16S_FR_Run1_afterQtrim-dada2_rep_set.qza \
  --i-data dada2_output/16S_FR_Run2_afterQtrim-dada2_rep_set.qza \
  --o-merged-data Runs_merged_output/dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_rep_set.qza

## Summarize the merged data
qiime feature-table summarize \
  --i-table Runs_merged_output/dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_table.qza \
  --o-visualization Runs_merged_output/dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_table.qzv \
  --m-sample-metadata-file ../16S_FR_Run1Run2_manifest.tsv
```



Take a look at merged table to see where the biggest drop off in sequences are in order to pick the --p-sampling-depth in `qiime diversity core-metrics-phylogenetic` step down below ("Alpha and beta diversity analysis on merged data")




## Estimate a phylogeny for the diversity stats step on merged data
This takes 2 minutes


```
## Get some new directories set
cd Runs_merged_output
mkdir phylogeny

# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/01_batch_effects/Runs_merged_output
```


```
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 16 \
  --i-sequences dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_rep_set.qza \
  --o-alignment phylogeny/16S_FR_Run1_Run2_afterQtrim-aligned_merged_dada2_rep_set.qza \
  --o-masked-alignment phylogeny/16S_FR_Run1_Run2_afterQtrim-masked_merged_dada2_rep_set.qza \
  --o-tree phylogeny/16S_FR_Run1_Run2_afterQtrim-unrooted_tree.qza \
  --o-rooted-tree phylogeny/16S_FR_Run1_Run2_afterQtrim-rooted_tree.qza
```



## Alpha and beta diversity analysis on merged data
This will give us PCoA plots to check for batch effects

```
# First using --p-sampling-depth of 10348 because that is lowest value before the biggest drop off

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny/16S_FR_Run1_Run2_afterQtrim-rooted_tree.qza \
  --i-table dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_table.qza \
  --p-sampling-depth 10348 \
  --m-metadata-file ../../16S_FR_Run1Run2_manifest.tsv \
  --output-dir merged-core-metrics-results-10348
```



```
# Let's try a lower depth so we can see more samples from both runs

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny phylogeny/16S_FR_Run1_Run2_afterQtrim-rooted_tree.qza \
  --i-table dada2_output/16S_FR_Run1_Run2_afterQtrim-merged_dada2_table.qza \
  --p-sampling-depth 1028 \
  --m-metadata-file ../../16S_FR_Run1Run2_manifest.tsv \
  --output-dir merged-core-metrics-results-1028
```

-  No batch effects! Can tell because the plots do not separate "Run1" and "Run2" samples in PCA space

### Let's move forward to the next directory! `02_Run1Run2_combo`

