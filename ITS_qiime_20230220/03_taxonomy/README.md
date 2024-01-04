### If needed, start an interactive srun on the AHPCC

```
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --partition pcon06 --time 6:00:00 --pty /bin/bash
```

### Load modules and check things

```
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate qiime2-2022.2

# cd to working directory with all files
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
```


First need to do some stuff to the newest UNITE database... Let's do a couple of things in this step, including training the newest classifier (for the step after deblur). I think I might be able to download these that are already pre-trained, but I'm not really sure what is wrong with the code when I do that. So just gonna retrain them myself.



#### Let's get the correctly formatted qza file from Colin Brislawn (https://github.com/colinbrislawn/unite-train/releases). 
1. Go to his source code (zipfile: `unite-train-9.0-qiime2-2022.11-demo.zip`) and navigate to the snakemake file (->workflow->snakefile)
2. Find the correct doi for wget'ing (can double-check with the public DOI link that is in his snakefile: https://api.plutof.ut.ee/v1/public/dois/?format=api&identifier=10.15156/BIO/2483915)
3. Looks like I want the 'normal' one

Download the files and set up some directories

```
cd /scrfs/storage/amatthews/20220621_ITS/UNITE
wget https://files.plutof.ut.ee/public/orig/59/12/591225E8985EFC44B595C79AF5F467421B4D9A95093A0811B13CB4CC13A6DA46.tgz
tar xzf 591225E8985EFC44B595C79AF5F467421B4D9A95093A0811B13CB4CC13A6DA46.tgz

# should have probably made a new directory before doing the wget :) 
mkdir sh_qiime_release_29.11.2022
mv *.fasta sh_qiime_release_29.11.2022
mv *.txt sh_qiime_release_29.11.2022
mv developer sh_qiime_release_29.11.2022

cd sh_qiime_release_29.11.2022/developer
```

```
# import the new UNITE reference seqs

qiime tools import \
  --type FeatureData[Sequence] \
  --input-path sh_refs_qiime_ver9_dynamic_29.11.2022_dev.fasta \
  --output-path unite-ver9-seqs_dynamic_29.11.2022.qza
```

```
# import the taxonomy file

qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path sh_taxonomy_qiime_ver9_dynamic_29.11.2022_dev.txt \
  --output-path unite-ver9-taxonomy_dynamic_29.11.2022.qza \
  --input-format HeaderlessTSVTaxonomyFormat
```

```
# train the classifier; takes 25 min
## this step is doing the taxonomic classification of the ASVs.

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite-ver9-seqs_dynamic_29.11.2022.qza \
  --i-reference-taxonomy unite-ver9-taxonomy_dynamic_29.11.2022.qza \
  --o-classifier unite-ver9-classifier_dynamic_29.11.2022.qza
```

```
# Convert the rep seqs from .fasta to .qza (of type FeatureData[Sequence])
cd ..
# pwd: /scrfs/storage/amatthews/20220621_ITS/UNITE/sh_qiime_release_29.11.2022

qiime tools import \
  --input-path sh_refs_qiime_ver9_dynamic_29.11.2022.fasta \
  --output-path sh_refs_qiime_ver9_dynamic_29.11.2022.qza \
  --type FeatureData[Sequence]
```

---









#### Now go on and run the classifier...

```
pwd: /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

mkdir 03_taxonomy
```

```
qiime feature-classifier classify-sklearn \
  --p-n-jobs 16 \
  --i-classifier ../UNITE/sh_qiime_release_29.11.2022/developer/unite-ver9-classifier_dynamic_29.11.2022.qza \
  --i-reads 02_dada2/ITS_forward_afterQtrim-pmin1-dada2_repseqs.qza \
  --o-classification 03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza
```










---
### Filter out mitochondrial, chloroplast, etc. from taxonomy

```
cd 03_taxonomy

mkdir only_fungi

qiime taxa filter-table \
   --i-table ../02_dada2/ITS_forward_afterQtrim-pmin1-dada2_table.qza \
   --i-taxonomy taxonomy-unite_ver9_dynamic-ITS_forward.qza \
   --p-exclude archaea,chloroplast,mitochondria,prokaryota \
   --o-filtered-table only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qza 

qiime feature-table summarize \
  --i-table only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qza \
  --o-visualization only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qzv \
  --m-sample-metadata-file ../ITS_forward_metadata.tsv
```


### Just check to see how things look; this is not the final taxa barplot!

```
qiime taxa barplot \
  --i-table only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi.qza \
  --i- taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --m-metadata-file ../ITS_forward_metadata.tsv \
  --o-visualization only_fungi/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-barplot.qzv
```

Things look pretty good, nothing seems fishy. Feathers have more diversity on them than mites do (which is probably expected). The controls have very little diversity (good!). So let's move on to decontam and get rid of contaminants, then we'll remake the barplot with just the biological samples.





---

### Move onto decontam!
