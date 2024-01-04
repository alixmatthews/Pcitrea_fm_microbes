#### MaAsLin2 for ITS
#### Alix Matthews
#### 18 July 2023

### R version 4.2.3 (2023-03-15) "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

#### These .csv files were taken from Qiime2 output and ran through the RelativeAbundanceCalculations_$GENE_$DATE.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for MaAsLin2

setwd("~/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/10_Maaslin/20230425")

library(dplyr) # v. 1.0.10
library(Maaslin2) # v. 1.14.1
library(ggplot2) # v. 3.4.0
library(ggpubr) # v. 0.5.0

# 0. Rows/columns to remove for all datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove <- c('PM13', 'PM03', 'PM06', 'PF03', 'PM09', 'PM08', 'PM04', 'PF25', 'PF10', 'PF13', 'PF30')

# 1. Load and Curate Sample Metadata ####

# This file is from the master results database used to create alpha/beta diversity metrics figures.
metadata <- read.csv(file = "/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/09_R/ITS_forward_master_results_20230413.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)
str(metadata)

# Remove control samples and exclude biological samples that were not included in the alpha/beta stats (i.e., had less than 500 seqs). So, in other words, these samples have no missing data in the metadata file...  "nm" = "no missing"

# drop controls
metadata <- filter(metadata, !bio_or_control == "control")

# drop samples that have <500 feature counts
metadata <- filter(metadata, !rownames(metadata) %in% samples_to_remove)

# change variable types

metadata[,c(8:9, 11:31)] <- sapply(metadata[,c(8:9, 11:31)],as.numeric)
# then take care of the chr->factor
metadata[sapply(metadata,is.character)] <- lapply(metadata[sapply(metadata, is.character)], as.factor)

str(metadata)




# 2. Load Data ####

#### These data include only biological samples that were included in the alpha/beta diversity stats (no missing samples = nm) and have undergone some manual curation specifically for MaAsLin (removing taxa that are unidentified at the taxa level of interest, combining taxa, etc... all is outlined in the RelativeAbundanceCalculations_$GENE_$DATE.R script)


## 2.1. Phylum ####

df_rel_abun_phylum_nm <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)
str(df_rel_abun_phylum_nm)


## 2.2. Family ####

df_rel_abun_family_nm <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)
str(df_rel_abun_family_nm)


## 2.3. Genus ####

df_rel_abun_genus_nm <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)
str(df_rel_abun_genus_nm)


# 3. Run MaAsLin2 - Phylum ####
## 3.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_phylum <-
  Maaslin2(
    input_data = df_rel_abun_phylum_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_phylum_20230718", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )

#### Results: 2 fungal phyla are differentially abundant using FDR-correct p-value









# 4. Run MaAsLin2 - Family ####
## 4.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_family <-
  Maaslin2(
    input_data = df_rel_abun_family_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_family_20230718", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )













# 5. Run MaAsLin2 - Genus ####
## 5.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_genus <-
  Maaslin2(
    input_data = df_rel_abun_genus_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_genus_20230718", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )

#### Results: genera that was differential abundant between feathers and mites as significant at the q-value (**** = > 10% and move forward with). None are >10% here...
