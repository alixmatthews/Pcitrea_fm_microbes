Now that we know that we don't have batch effects, we can safely move forward with combining Run1 and Run2 data into a single analysis. It is possible to merge these two runs within QIIME2 (using the 'group'/'collapse' method), but it is risky to do so. The reason it is risky is because our DADA2 runs were done separately, so some of the features that DADA2 identified could be repeated between runs but the feature-id could be different (unknowingly, since they are just a bunch of numbers and letters), so it would artifically inflate the feature counts for each sample. The collapse method may account for this, but it is not clear how it would do so (or if it does so). So we'll consider all samples as one big run, analyze them together, and then collapse the samples. This is the safer bet and is a bit cleaner.

#### Data files in this directory and used for this step:
`16S_FR_Run1Run2_manifest.tsv`: concatenated version of Run1 and Run2 manifest files containing metadata and filepaths

---

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
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo
```

## Import and summarize the combined runs
```
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path 16S_FR_Run1Run2_manifest.tsv \
  --output-path 16S_FR_Run1Run2_beforeQtrim-demux.qza

qiime demux summarize \
--i-data 16S_FR_Run1Run2_beforeQtrim-demux.qza \
--o-visualization 16S_FR_Run1Run2_beforeQtrim-demux.qzv
```



## Trim primers and summarize combined trimmed data
Trimming primers will also consequently trim Illumina adapters

These primer seqs were obtained from the sequencing facility. They are the 515F/806R primer set for 16S (V4 region).

This takes less than 5 minutes

```
## trim 
qiime cutadapt trim-paired \
  --p-cores 16 \
  --p-front-f ^GTGCCAGCMGCCGCGGTAA \
  --p-front-r ^GGACTACHVGGGTWTCTAAT \
  --p-error-rate 0.10 \
  --p-minimum-length 100 \
  --p-discard-untrimmed \
  --i-demultiplexed-sequences 16S_FR_Run1Run2_beforeQtrim-demux.qza \
  --o-trimmed-sequences 16S_FR_Run1Run2_afterQtrim-demux.qza

## summarize
qiime demux summarize \
--i-data 16S_FR_Run1Run2_afterQtrim-demux.qza \
--o-visualization 16S_FR_Run1Run2_afterQtrim-demux.qzv
```




## DADA2 on the trimmed (still paired) data for runs combo
This step will merge reads and remove poor quality bases based on parameters I choose.

Parameters in this step were either adjusted incrementally/experimentally or chosen from the `*afterQtrim-demux.qzv` output files. Final parameters were chosen based on examining the results in the stats tables files that are generated in this step. I chose parameters that allowed the most reads to be retained as possible and fewest proportion of chimeras. All parameters in this step should be chosen carefully. 

- `--p-chimera-method`: experiment with this (default: consensus)
- `--p-trunc-len-f`: check `*afterQtrim-demux.qzv` to determine
- `--p-trunc-len-r`: check `*afterQtrim-demux.qzv` to determine
- `--p-min-overlap`: experiment with this (default: 12)
- `--p-min-fold-parent-over-abundance`: experiment with this (default: 1.0)
- `--p-trunc-q`: experiment with this (default: 2)

```
## make new dir for dada2 output
mkdir dada2_output
```

```
## DADA2 - default params (takes about 5 min to run)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S_FR_Run1Run2_afterQtrim-demux.qza \
  --p-n-threads 16 \
  --p-chimera-method consensus \
  --p-hashed-feature-ids \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 130 \
  --p-min-overlap 12 \
  --p-min-fold-parent-over-abundance 1.0 \
  --p-trunc-q 2 \
  --o-table dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_default_table.qza \
  --o-representative-sequences dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_default_rep-seq.qza \
  --o-denoising-stats dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_default_stats.qza
  
## Generate stats tables - default params (takes about 1 min)
qiime metadata tabulate \
  --m-input-file dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_default_stats.qza \
  --o-visualization dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_default_stats.qzv
```


```
## DADA2 - adjusted params (takes about 5 min to run)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 16S_FR_Run1Run2_afterQtrim-demux.qza \
  --p-n-threads 16 \
  --p-chimera-method consensus \
  --p-hashed-feature-ids \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 130 \
  --p-min-overlap 4 \
  --p-min-fold-parent-over-abundance 4.0 \
  --p-trunc-q 15 \
  --o-table dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_table.qza \
  --o-representative-sequences dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
  --o-denoising-stats dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_stats.qza
  
## Generate stats tables - adjusted params
qiime metadata tabulate \
  --m-input-file dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_stats.qza \
  --o-visualization dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_stats.qzv
```




## Group (collapse) Run1 and Run2 data
Put them all in a single table (qza) to use for downstream analyses

Parameters will depend on how you want to group your data, but here my data are grouped by a column called `sample_group_id`, so I am collapsing each sample that has the sample `sample_group_id` identifier into a single data point.


```
# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo

qiime feature-table group \
  --i-table dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_table.qza \
  --p-axis sample \
  --m-metadata-file 16S_FR_Run1Run2_manifest.tsv \
  --m-metadata-column sample_group_id \
  --p-mode sum \
  --o-grouped-table dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped.qza
```


## Summarize the grouped data 
Note that the metadata file used in this step is ONLY metadata (not like the manifest file which has the file-paths)

```
qiime feature-table summarize \
  --i-table dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped.qza \
  --o-visualization dada2_output/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped.qzv \
  --m-sample-metadata-file ../16S_sample_metadata.tsv
```




## Estimate a phylogeny for the diversity stats on grouped data

```
## Get some new directories set
cd dada2_output
mkdir phylogeny

# pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output
```


```
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 16 \
  --i-sequences 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
  --o-alignment phylogeny/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-aligned.qza \
  --o-masked-alignment phylogeny/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-masked.qza \
  --o-tree phylogeny/16S_FR_Run1Run2_afterQtrim-unrooted_tree.qza \
  --o-rooted-tree phylogeny/16S_FR_Run1Run2_afterQtrim-rooted_tree.qza
```





## SILVA
First need to download the newest .qza classifer (Naive Bayes) for our dataset from qiime2 data-resources page: `Silva 138 99% OTUs from 515F/806R region of sequences (MD5: e05afad0fe87542704be96ff483824d4)`

```
Check to see how many cores the node I'm in has to maximize the efficiency

re-ssh into AHPCC in another window
ssh c1421 # change to whatever node
htop
# there are 24 cores in this node, but someone is using 4 right now
```

```
# Then let's make a new directory:
mkdir taxonomy

pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output
```

This step takes ~10 minutes on the AHPCC

```
qiime feature-classifier classify-sklearn \
   --p-n-jobs 20 \
   --i-classifier ../../silva-138-99-515-806-nb-classifier.qza \
   --i-reads 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq.qza \
   --o-classification taxonomy/taxonomy_16S_grouped-silva.qza
```



## Filter out all eukaryota, mitochondrial, chloroplast, and archaeal seqs from my taxonomy

```
mkdir only_bacteria

pwd: /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/02_Run1Run2_combo/dada2_output
```


```
qiime taxa filter-table \
   --i-table 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped.qza \
   --i-taxonomy taxonomy/taxonomy_16S_grouped-silva.qza \
   --p-exclude archaea,chloroplast,mitochondria,eukaryota \
   --o-filtered-table only_bacteria/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered.qza 
```

```
## Summarize the filtered samples data table
qiime feature-table summarize \
  --i-table only_bacteria/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered.qza \
  --o-visualization only_bacteria/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered.qzv \
  --m-sample-metadata-file ../../16S_sample_metadata.tsv
```


Looks like we lost 3 samples between the grouping and the taxonomy here. 3 samples had zero features in the previous dataset (`16S_FR_Run1Run2_afterQtrim-dada2_table-grouped`), so this makes sense. 



### Now we're ready to move on to Decontam! (`03_decontam`)
