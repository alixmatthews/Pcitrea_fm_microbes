# Taxa barplots

This step is to create taxa barplots and then use R to make them look better.


#### Data files in this directory and used for this step:

- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2.csv`: exported results from qiime taxa barplot at phylum level

- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5.csv`: exported results from qiime taxa barplot at family level

- `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6.csv`: exported results from qiime taxa barplot at genus level

- `TaxaBarplots_R/500plus`: R scripts for visualization

**NOTE:** for each of the `.csv` files listed above, the first column (`index`) is the sample identifier, while the remaining columns are each phylogenetic unit at specific taxonomic levels (phylum, family, genus, respectively). The values within the matrix are read count values (count data)

---

This should be on un-rarefied data, biological samples only, after removing contaminants. In other words, do not use rarefied data for taxa barplots.

```
cd /scrfs/storage/amatthews/20220621_ITS/01_qiime_20230220
mkdir 07_TaxaBarplots
```

Go for it!

```
qiime taxa barplot \
  --i-table 05_SeparateSamples/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples.qza \
  --i-taxonomy 03_taxonomy/taxonomy-unite_ver9_dynamic-ITS_forward.qza \
  --m-metadata-file ITS_forward_metadata.tsv \
  --o-visualization 07_TaxaBarplots/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-barplot.qzv
```


## Import the .qzv file onto Qiime2View and then download the .csv files for each level desired.

- Phylum: `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2.csv`
- Family: `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5.csv`
- Genus: `ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6.csv`




## Then run through R scripts (only for 500 plus data)
Can be found in `./TaxaBarplots_R/500plus` -- also includes 16S in this same script


R scripts: 
- `RelativeAbundance_ITS_20230421.R`: use this one with the output from the Qiime2View .csv files (above) to calculate relative abundance at each of the taxonomic levels. This also goes through the curation of the files for the taxonomic barplots. First, remove samples that have <500 feature counts. Second, curate features (e.g., combine certain features that should be combined but are not due to the inconsistency with SILVA taxonomy) for taxonomic barplots. Third, curate features for MaAsLin2, an upcoming step (e.g., remove features that are not identified at the level we are interested in, such as only identified to Order, but we are interested in looking at Family).
- `16S_ITS_COMBO_taxabarplots_20230903.R`: this uses the curated .csv files at the second curation step from the `RelativeAbundance_ITS_20230421.R` file and makes the taxonomic barplots. This includes both 16S and ITS data for ease of editing downstream, although this particular file is in the ITS directory.
-`ITS_taxabarplots_20230501.R`: ITS taxabarplots making (same as above, but only ITS)



