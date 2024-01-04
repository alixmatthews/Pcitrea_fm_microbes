Steps for ITS data analyses (cleaned up from private repos)

- `ITS_forward_metadata.tsv` is the original metadata file. Modified in some downstream steps (e.g., for plotting)
- `ReadCounts_ITS.xlsx` shows the change in the number of reads per samples at each step of pre-processing data

Forward reads only instead of trying to use PEAR to pair reads (poor quality overlap lost a lot of reads) and then using deblur. Forward reads only allows us to use DADA2, which more closely matches with 16S analyses.
