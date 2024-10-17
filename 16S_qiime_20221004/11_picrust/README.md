# PiCrust
#### Note: these analyses were not presented in ms

Another note: subdirectory #10 was a trial with q2-longitudinal. Didn't produce anything valuable for these data and wasn't quite what we were going for. So removed it, but have kept the numbering as #11 for picrust in order to stay consistent with my messier notes.

Standalone PiCrust2 worked better than qiime2 plugin at the time of these analyses. More flexibility on output with the standalone version.
 
Inputs and some of the smaller outputs included in the subdirectory here (`picrust2_out_pipeline_20230606`). Files in that subdirectory include:
- `16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.biom`: input for Picrust
- `16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus_out.tre`: tree after seqs had been placed by PiCrust
- `path_abun_unstrat.tsv.gz`: after PiCrust had been run, this was the output of the pathways
- `path_abun_unstrat_descrip.tsv.gz`: output of pathways after PiCrust had been run and includes descriptions of what the pathways are.


## 0. Install PiCrust2 on AHPCC (and my own machine)

```
# on AHPCC, standalone
module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
mamba create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2
mamba activate picrust2

# on own machine, standalone
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh"
bash Mambaforge-pypy3-Linux-x86_64.sh
## close and reopen shell
mamba create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2
mamba activate picrust2

```


---

## 1. Try tutorial with own data (AHPCC standalone) 

### 1a. Convert some files
```
# first need to convert repseqs to .fna/fasta file (=> ASV seqs)... make sure qiime environment is still activated

cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples

qiime tools export \
    --input-path 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.qza \
    --output-path Rep-seqs/

cd Rep-seqs

mv dna-sequences.fasta 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.fasta

# now need to convert table to biom.tsv table (=> ASV by sample)

cd ..

qiime tools export \
    --input-path 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.qza \
    --output-path Table-export/
    
cd Table-export

mv feature-table.biom 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.biom

biom convert --to-tsv -i 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.biom -o 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.tsv

```

### 1b. Take converted files and use for picrust standalone

```
#module load python
#module list
#source /share/apps/bin/conda-3.11.sh

conda deactivate #deactivating qiime
conda activate picrust2

biom summarize-table -i 16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.biom

cd /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/11_picrust

mkdir picrust2_out_pipeline_20230606 #to differentiate from the q2-picrust2 plugin
cd picrust2_out_pipeline_20230606/
```

### 1c. Place sequences 
```
place_seqs.py --study_fasta /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples/Rep-seqs/16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus.fasta \
--placement_tool epa-ng \
--out_tree 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus_out.tre \
--processes 32 \
--intermediate intermediate/place_seqs
```

### 1d. Hidden-state prediction of gene families
```

# First take a look at 16S... i.e., the predicted number of 16S copies per ASVs, followed by the nearest-sequenced taxon index (NSTI) value per ASV. This metric is a rough guide for how similar an ASV is to an existing reference sequence.

hsp.py --in_trait 16S \
--tree 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus_out.tre \
--output marker_predicted_and_nsti.tsv.gz \
--processes 1 \
--calculate_NSTI

# Now do EC... i.e., the predicted copy number of all Enzyme Classification (EC) numbers is shown for each ASV. EC numbers are a type of gene family defined based on the chemical reactions they catalyze. EC numbers can be used to infer MetaCyc pathway levels.

hsp.py --in_trait EC \
--tree 16S_FR_Run1Run2_afterQtrim-dada2_rep_seq-grouped-filtered-decontam-bio_samples_500plus_out.tre \
--output EC_predicted.tsv.gz \
--processes 1

```

### 1e. Generate metagenome predictions

```
# Here, we are interested in the predicted gene families weighted by the relative abundance of ASVs in their community.
# First, the read depth per ASV is divided by the predicted 16S copy numbers. Then, the ASV read depths per sample (after normalizing by 16S copy number) are multiplied by the predicted gene family copy numbers per ASV. --max_nsti of 2 is default, just including for a reminder...


metagenome_pipeline.py --input /scrfs/storage/amatthews/20220523_16S/01_qiime_20221004/04_SeparateSamples/Table-export/16S_FR_Run1Run2_afterQtrim-dada2_table-grouped-filtered-decontam-bio_samples_500plus.biom \
--marker marker_predicted_and_nsti.tsv.gz \
--function EC_predicted.tsv.gz \
--max_nsti 2 \
--out_dir EC_metagenome_out \
--strat_out

```

### 1f. Pathway-level inference

```
# By default this script infers MetaCyc pathway abundances based on EC number abundances, although different gene families and pathways can also be optionally specified... specifically, Regroups EC numbers to MetaCyc reactions; infers which MetaCyc pathways are present based on these reactions with MinPath; Calculates and returns the abundance of pathways identified as present

pathway_pipeline.py --input EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
--out_dir pathways_out \
--processes 2

# Output: include both unstratified and stratified MetaCyc pathway abundances, similar to the EC number tables above. The default stratified pathway abundance table represents how much each ASV is contributing to the community-wide pathway abundance and not what the pathway abundance is predicted to be within the predicted genome of that ASV alone. 
```



### 1g. Add functional descriptions

```
# a description of each functional id in the output abundance tables. The below commands will add these descriptions as new column in gene family and pathway abundance tables

# first on unstratified (metagenome)
add_descriptions.py --input EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
--map_type EC \
--output EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py --input pathways_out/path_abun_unstrat.tsv.gz \
--map_type METACYC \
--output pathways_out/path_abun_unstrat_descrip.tsv.gz

```

### 1h. Downstream analyses

```
# some options `ggpicrust2`, `ALDEx2`, `STAMP`
# shuffle_predictions.py (https://github.com/picrust/picrust2/wiki/Generating-shuffled-predictions)
```
