# FUNguild

get going on HPC

```
srun -N1 -n1 -c1 -p cloud72 -q cloud -t 72:00:00 --pty /bin/bash 
module load gcc/11.2.1 mkl/20.0.4 python/3.11-anaconda
. /share/apps/bin/conda-3.11.sh
conda activate qiime2-2023.2
```

convert feature table to biom
```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220

qiime tools export \
  --input-path 05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples_500plus.qza \
  --output-path 05_SeparateSamples/exported-feature-table
```

convert taxonomy to biom 
```
qiime tools export \
 --input-path 03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
 --output-path 03_taxonomy/taxonomy

cp 03_taxonomy/taxonomy/taxonomy.tsv 03_taxonomy/taxonomy/biom-taxonomy.tsv

nano 03_taxonomy/taxonomy/biom-taxonomy.tsv

# changed header to:  #OTUID  taxonomy        confidence
```

Add metadata (taxonomy)
```
biom add-metadata \
--input-fp 05_SeparateSamples/exported-feature-table/feature-table.biom \
--output-fp 11_FUNguild/table-biom-with-taxonomy.biom \
--observation-metadata-fp 03_taxonomy/taxonomy/biom-taxonomy.tsv \
--sc-separated taxonomy
```

convert biom to tsv
```
biom convert \
--input-fp 11_FUNguild/table-biom-with-taxonomy.biom \
--output-fp 11_FUNguild/OTU-table-frombiom-with-taxonomy.tsv \
--to-tsv \
--header-key taxonomy
```

Output here: `OTU-table-frombiom-with-taxonomy.tsv`
Converted to .txt here: `OTU-table-frombiom-with-taxonomy.txt`

FunGuild
- copied and pasted the .py script to 11_FUNguild directory
- python 3.11 is already active... so can just run this as is apparently...

```
cd 11_FUNguild
cp OTU-table-frombiom-with-taxonomy.tsv OTU-table-frombiom-with-taxonomy.txt
nano OTU-table-frombiom-with-taxonomy.txt
# removed the `# Constructed from biom file` and changed #OTU ID to OTU_ID
```

Run FUNguilds

Python file: `Guilds_v1.0.py`


```
python Guilds_v1.0.py -otu OTU-table-frombiom-with-taxonomy.txt -db fungi -m -u
```

It worked!



Run FUNguilds version 1.1 instead

Python file: `Guilds_v1.1.py`

```
python Guilds_v1.1.py -otu OTU-table-frombiom-with-taxonomy.txt -db fungi -m -u
```


This also worked!




---
## MaAsLin2 analyses post-FUNguild
#### note that the Maaslin2 version changed (from 1.10.0 to 1.14.1) because I updated R to 4.3 and apparently also updated Maaslin! This update to 1.14.1 required a change in the import step (e.g., `row.names = 1, stringsAsFactors = FALSE), so just keep that in mind if trying to run through Maaslin again...


- Reduce the matrix `OTU-table-frombiom-with-taxonomy.guilds.txt` to only include animal pathogens that are highly probable or probable: `OTU-table-frombiom-with-taxonomy.guilds-animalpathogens.txt` 
- Then transpose it: `OTU-table-frombiom-with-taxonomy.guilds-animalpathogens_transposed.txt` 
- Calculate relative abundance: `RelativeAbundance_funguild_20230624.R`
- MaAsLin2: `MaAsLin2_FUNguild_20230624.R` ... this actually does not work because every ASV is unique and MaAsLin2 filters out all of the features before running the analysis, so sadly cannot do this.
  - but here are the other files if you want to replicate it or try again: `FUNguild_metadata.csv` and `OTU-table-frombiom-with-taxonomy.guilds-animalpathogens_transposed_relabun.csv`
 
  
---

## Just a visual on the different trophic modes for the ms
- `OTU-table-frombiom-with-taxonomy.guilds.txt`
- `FUNguild-20230615.R`

