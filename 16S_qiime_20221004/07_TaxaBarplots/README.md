# Taxa barplots

This step is to create taxa barplots and then use R to make them look better.


#### Data files in this directory and used for this step:

- `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-2.csv`: exported results from qiime taxa barplot at phylum level

- `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-5.csv`: exported results from qiime taxa barplot at family level

- `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-6.csv`: exported results from qiime taxa barplot at genus level

- `TaxaBarplots_R/500plus`: R scripts for visualization

**NOTE:** for each of the `.csv` files listed above, the first column (`index`) is the sample identifier, while the remaining columns are each phylogenetic unit at specific taxonomic levels (phylum, family, genus, respectively). The values within the matrix are read count values (count data)

---


Not using rarefied data for taxa barplots... all data with seqs >500


```
cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004
mkdir 07_TaxaBarplots
cd 07_TaxaBarplots
```

Go for it!

```
qiime taxa barplot \
  --i-table ../04_SeparateSamples/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples.qza \
  --i-taxonomy ../02_Run1Run2_combo/dada2_output/taxonomy/taxonomy_16S_grouped-silva.qza \
  --m-metadata-file ../16S_sample_metadata.tsv \
  --o-visualization taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva.qzv
```

## Import the .qzv file onto Qiime2View and then download the .csv files for each level desired.

Phylum: `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-2.csv`
Family: `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-5.csv`
Genus: `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-6.csv`


## Then run through R scripts (only for 500 plus data)
Can be found in `./TaxaBarplots_R/500plus` -- also includes ITS in this same script

R scripts: 
- `RelativeAbundance_16S_20230421.R`: use this one with the output from the Qiime2View .csv files (above) to calculate relative abundance at each of the taxonomic levels. This also goes through the curation of the files for the taxonomic barplots. First, remove samples that have <500 feature counts. Second, curate features (e.g., combine certain features that should be combined but are not due to the inconsistency with SILVA taxonomy) for taxonomic barplots. Third, curate features for MaAsLin2, an upcoming step (e.g., remove features that are not identified at the level we are interested in, such as only identified to Order, but we are interested in looking at Family).
- `16S_ITS_COMBO_taxabarplots_20230903.R`: this uses the curated .csv files at the second curation step from the `RelativeAbundance_16S_20230421.R` file and makes the taxonomic barplots. This includes both 16S and ITS data for ease of editing downstream, although this particular file is in the 16S directory.






