### Decontam (ITS seqs, forward reads only)
### Alix Matthews
### 2023 February 20

### R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"

### OS: Ubuntu 20.04.3 LTS


#### Set working directory, install packges, load libraries ####
setwd("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/04_decontam")

# only install these once per machine, comment out otherwise
# install.packages("BiocManager") 
# BiocManager::install("GenomeInfoDbData", force = TRUE)
# BiocManager::install("phyloseq", force = TRUE)
# BiocManager::install("decontam", force = TRUE)

library(BiocManager) # 1.30.18
library(decontam) # 1.16.0
# packageVersion("decontam")
library(phyloseq) # 1.40.0
# packageVersion("phyloseq")
library(ggplot2) # 3.3.6
library(ggthemes) # 4.2.4

### Import biom file and metadata, make phyloseq object ####
biom_otu_tax <- import_biom(BIOMfilename = "/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/03_taxonomy/only_fungi/exported-feature-table/ITS_fungi-all_samples-feature_table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/ITS_forward_metadata.tsv")

mite.phyloseq <- merge_phyloseq(biom_otu_tax, metadata)
mite.phyloseq

head(sample_data(mite.phyloseq))

#### Convert phyloseq to data frame ####
df <- as.data.frame(sample_data(mite.phyloseq))
str(df)

#### Figure of library size x sample type ####
df$LibrarySize <- sample_sums(mite.phyloseq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

ggplot(data=df, aes(x=Index, y=LibrarySize, color=bio_or_control)) + 
  geom_point(size=2.5, alpha = 0.75) +
  theme_clean() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 10)) +
  theme(legend.text = element_text(face = "bold", size = 10)) +
  scale_color_manual(name = "Sample Type", labels = c("Biological", "Control"), values = c("goldenrod","seashell4"))

# saved as A6 landscape


#### Run decontam to identify potential contaminants ####
sample_data(mite.phyloseq)$is.neg <- sample_data(mite.phyloseq)$bio_or_control == "control"


# Make phyloseq object of presence-absence in negative controls and true samples
mite.phyloseq.pa <- transform_sample_counts(mite.phyloseq, function(abund) 1*(abund>0))
mite.phyloseq.pa.neg <- prune_samples(sample_data(mite.phyloseq.pa)$bio_or_control == "control", mite.phyloseq.pa)
mite.phyloseq.pa.pos <- prune_samples(sample_data(mite.phyloseq.pa)$bio_or_control == "biological", mite.phyloseq.pa)


#### + try different thresholds ####
#### +++ 0.1 ####
contamdf.prev01 <- isContaminant(mite.phyloseq, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev01$contaminant)
head(which(contamdf.prev01$contaminant))


# Make data.frame of prevalence in positive and negative samples
df.pa01 <- data.frame(pa.pos=taxa_sums(mite.phyloseq.pa.pos), pa.neg=taxa_sums(mite.phyloseq.pa.neg),
                      contaminant=contamdf.prev01$contaminant)

ggplot(data=df.pa01, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  #geom_point(size=2.5, alpha = 0.75) +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  theme_clean() +
  geom_jitter(size=2.5, alpha = 0.75) +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 10)) +
  theme(legend.text = element_text(face = "bold", size = 10)) +  
  scale_color_manual(name = "Contaminant", labels = c("False", "True"), values = c("orchid","red4")) +
  ggtitle("Prevalence threshold = 0.1")

# saved as A6 landscape


#### +++ 0.5 ####
contamdf.prev05 <- isContaminant(mite.phyloseq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))


# Make data.frame of prevalence in positive and negative samples
df.pa05 <- data.frame(pa.pos=taxa_sums(mite.phyloseq.pa.pos), pa.neg=taxa_sums(mite.phyloseq.pa.neg),
                      contaminant=contamdf.prev05$contaminant)

ggplot(data=df.pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  #geom_point(size=2.5, alpha = 0.75) +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  theme_clean() +
  geom_jitter(size=2.5, alpha = 0.75) +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.title = element_text(face = "bold", size = 10)) +
  theme(legend.text = element_text(face = "bold", size = 10)) +  
  scale_color_manual(name = "Contaminant", labels = c("False", "True"), values = c("orchid","red4")) +
  ggtitle("Prevalence threshold = 0.5")



# saved as A6 landscape


#### Write out the results for final analysis ####
write.csv(df.pa05, "./DecontamResults_prev05.csv")
