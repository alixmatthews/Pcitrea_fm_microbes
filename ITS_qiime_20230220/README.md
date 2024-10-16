Steps for ITS data analyses (cleaned up from private repos)

#### Files in this directory
- `ITS_forward_metadata.tsv` is the original metadata file. Modified in some downstream steps (e.g., for plotting)
- `ReadCounts_ITS.xlsx` shows the change in the number of reads per samples at each step of pre-processing data

NOTE: For ITS analyses, only forward reads were used instead of trying to use PEAR to pair forward and reverse reads and then using deblur. Doing so (i.e., PEAR/deblur) resulted in very poor quality overlap and we lost a lot of reads. Using only forward reads allows us to use DADA2, which more closely matches with 16S analyses.

