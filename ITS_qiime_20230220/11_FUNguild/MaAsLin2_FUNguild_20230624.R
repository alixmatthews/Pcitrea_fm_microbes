#### MaAsLin2 for FUNguild results - THIS HAS ERRORS, EMAILED BIOBAKERY FOR HELP
#### Alix Matthews
#### 24 June 2023

### R version 4.3.0 (2023-04-21) "Already Tomorrow"

### OS: Ubuntu 20.04.3 LTS

#### These .csv files were taken from FUNguild output and ran through the RelativeAbundance_funguild_20230624.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for MaAsLin2

setwd("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/11_FUNguild/")

library(dplyr) # v. 1.1.2
library(Maaslin2) # v. 1.14.1
library(ggplot2) # v. 3.4.0
library(ggpubr) # v. 0.5.0

# 0. Rows/columns to remove for all datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove <- c('PM13', 'PM03', 'PM06', 'PF03', 'PM09', 'PM08', 'PM04', 'PF25', 'PF10', 'PF13', 'PF30')


# 1. Load and Curate Sample Metadata ####

# This file is from the master results database used to create alpha/beta diversity metrics figures.
metadata <- read.csv(file = "/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/09_R/ITS_forward_master_results_20230413.csv", sep = ",", header= TRUE)
str(metadata)

# Remove control samples and exclude biological samples that were not included in the alpha/beta stats (i.e., had less than 500 seqs). So, in other words, these samples have no missing data in the metadata file...  "nm" = "no missing"

# drop controls
metadata <- filter(metadata, !bio_or_control == "control")

# drop samples that have <500 feature counts
metadata <- filter(metadata, !sample.id %in% samples_to_remove)

# change variable types

metadata[,c(9:10, 12:32)] <- sapply(metadata[,c(9:10,12:32)],as.numeric)
# then take care of the chr->factor
metadata[sapply(metadata,is.character)] <- lapply(metadata[sapply(metadata, is.character)], as.factor)

str(metadata)

write.table(metadata , file = "FUNguild_metadata.csv", sep=",", row.names=FALSE)



# error with metadata??????
metadata <- read.csv(file = "/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/11_FUNguild/FUNguild_metadata_ordered.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)
str(metadata)

metadata <- as.data.frame(metadata)

row.names(metadata)


# 2. Load Data ####

#### These data were taken from the FUNguild output, transposed, and then converted to relative abundances in `RelativeAbundance_funguild_20230624.R` 

## 2.1. Animal Pathogens only Data ####


# error with data??????

data_guilds <- read.csv(file = "OTU-table-frombiom-with-taxonomy.guilds-animalpathogens_transposed_relabun.csv", sep = ",", header= TRUE, row.names = 1, stringsAsFactors = FALSE)

str(data_guilds)
data_guilds[sapply(data_guilds,is.integer)] <- lapply(data_guilds[sapply(data_guilds, is.integer)], as.numeric)

data_guilds <- as.data.frame(data_guilds)

row.names(data_guilds)


# 3. Run MaAsLin2 - Animal Pathogen only Data ####

## 3.1. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category <-
  Maaslin2(
    input_data = data_guilds, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )


# ERRORRRRR?????/