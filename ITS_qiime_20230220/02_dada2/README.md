#### Start an interactive srun on the AHPCC

```
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --partition pcon06 --time 6:00:00 --pty /bin/bash
```

#### Load modules and check things

```
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate qiime2-2022.2

cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
```


## DADA2 on the trimmed forward data
This step will remove poor quality bases based on parameters I choose.

Parameters in this step were either adjusted incrementally/experimentally or chosen from the `afterQtrim-demux.qzv` files. Final parameters were chosen based on examining the results in the stats tables files that are generated in this step. I chose parameters that allowed the most reads to be retained as possible and fewest proportion of chimeras. All parameters in this step should be chosen carefully.

- `--p-trunc-len`: check `afterQtrim-demux.qzv` to determine. I also tried p-trunc-q 20, but it had higher % of chimeras than p-trunc-q 15
- `--p-min-fold-parent-over-abundance`: experimented with this (no dif between 1.0 (default) and 4.0 (what I used for 16S), so going with the default (1.0))

```
## DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs 01_ForwardData/ITS_forward_afterQtrim-demux.qza \
  --p-trunc-len 275 \
  --p-chimera-method consensus \
  --p-trunc-q 15 \
  --p-min-fold-parent-over-abundance 1.0 \
  --p-n-threads 16 \
  --o-table 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_table.qza \
  --o-representative-sequences 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qza \
  --o-denoising-stats 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_stats.qza
  
## Generate stats tables for these new runs
qiime metadata tabulate \
  --m-input-file 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_stats.qza \
  --o-visualization 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_stats.qzv
```



```
qiime feature-table tabulate-seqs \
  --i-data 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qza \
  --o-visualization 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qzv
```
  
  


```
qiime feature-table summarize \
  --i-table 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_table.qza \
  --o-visualization 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_table.qzv \
  --m-sample-metadata-file ITS_forward_metadata.tsv
```






---



### Need to infer a phylogeny for diversity stats. *HOWEVER, BE AWARE*...
##### "Although this marker provides more accurate species identification than other fungal marker genes, alignments of ITS sequences are not useful for informing evolutionary distances among distantly related species (but can yield good alignments for closely related sequences). So you have two options for diversity and other potentially phylogenetic analyses: 1) do not use phylogenetic methods with your ITS data. Instead, use non-phylogenetic methods diversity analyses, such as Jaccard (qualitative) or Bray Curtis (quantitative) for estimating beta diversity. 2) Use q2-ghost-tree to infer phylogeny via grafting ITS sequences. This will enable use of phylogenetic methods, including UniFrac."
###### from here: https://forum.qiime2.org/t/fungal-its-analysis-tutorial/7351
##### Issue I see with ghost-tree is that it has not been updated since 2021 with an older version of UNITE, and the developer is no longer working on ghost-tree, so may have to go with the non-phylogenetic methods, and that will have to be that!

##### First set up some dirs and get in the right place

```
cd 02_dada2
mkdir phylogeny
cd ..
# pwd: /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
```

##### Then run it!

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --p-n-threads 16 \
  --i-sequences 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qza \
  --o-alignment 02_dada2/phylogeny/ITS_forward_afterQtrim-pmin1-dada2_repseqs-aligned.qza \
  --o-masked-alignment 02_dada2/phylogeny/ITS_forward_afterQtrim-pmin1-dada2_repseqs-masked.qza \
  --o-tree 02_dada2/phylogeny/ITS_forward_afterQtrim-pmin1-dada2_repseqs-unrooted_tree.qza \
  --o-rooted-tree 02_dada2/phylogeny/ITS_forward_afterQtrim-pmin1-dada2_repseqs-rooted_tree.qza
```

---

### Now can move on to the classifier step.









