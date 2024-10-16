This purpose of this step is to start to process the reads for downstream analyses (in this case, trim primers)

#### Data file in this directory and used for this step:
- `ITS_forward_manifest.tsv`: manifest file for QIIME2 input

---

### Start an interactive srun on the AHPCC

```
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=16 --partition pcon06 --time 6:00:00 --pty /bin/bash
```

### Load modules and check things

```
module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate qiime2-2022.2

# cd to working directory with all files and make new dir
cd /scrfs/storage/amatthews/20220621_ITS
mkdir 01_qiime_20230220
cd 01_qiime_20230220
```

---

### Import and summarize the forward reads

```
mkdir 01_ForwardData

qiime tools import \
  --type SampleData[SequencesWithQuality] \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ITS_forward_manifest.tsv \
  --output-path 01_ForwardData/ITS_forward_beforeQtrim-demux.qza

# summarize and view 
qiime demux summarize \
  --i-data 01_ForwardData/ITS_forward_beforeQtrim-demux.qza \
  --o-visualization 01_ForwardData/ITS_forward_beforeQtrim-demux.qzv
```







---

### Now can remove adapters with cutadapt
I want to trim the forward and reverse adapters. To check on this, I viewed the merged .fastq files and the p-front sequence is at the beginning of the reads and the p-anywhere close to the end. The reverse complements are not found, so only need to do this step once. I'm setting the p-minimum-length to 250 because that is when the reads start to get a  little worse. If I set it at 300, there are no reads left after trimming!

```
## run cutadapt
qiime cutadapt trim-single \
  --p-cores 16 \
  --i-demultiplexed-sequences 01_ForwardData/ITS_forward_beforeQtrim-demux.qza \
  --p-front ^CTTGGTCATTTAGAGGAAGTAA \
  --p-anywhere GCATCGATGAAGAACGCAGC \
  --p-match-read-wildcards \
  --p-match-adapter-wildcards \
  --p-minimum-length 250 \
  --p-discard-untrimmed \
  --o-trimmed-sequences 01_ForwardData/ITS_forward_afterQtrim-demux.qza
  
## summarize trimmed ITS forward seqs
qiime demux summarize \
  --i-data 01_ForwardData/ITS_forward_afterQtrim-demux.qza \
  --o-visualization 01_ForwardData/ITS_forward_afterQtrim-demux.qzv
```




Looks better and we will remove poor reads at the next step 



---

### Now can move onto the DADA2/UNITE steps (`02_dada2`)
