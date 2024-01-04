#### Calculate ITS relative abundance at phylum, family, and genus for import for taxa barplots and MaAsLin2
#### Alix Matthews
#### 21 April 2023

#### These .csv files were taken from Qiime2 output `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva.qzv` and downloading the .csv file for the specific level of interest.

setwd("~/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/07_TaxaBarplots/R/500plus")

library(dplyr) # v. 1.0.10

# 0. Rows/columns to remove for all datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove <- c('PM13', 'PM03', 'PM06', 'PF03', 'PM09', 'PM08', 'PM04', 'PF25', 'PF10', 'PF13', 'PF30')

# Unncessary columns to be dropped from every level
drop_cols <- c("bio_or_control", "category_broad", "category_specific", "bird_id", "bird_sex", "bird_age", "capture_date", "num_mites_fm", "num_mites_m", "mite_sp_fm", "mite_sp_m")

# 1. Phylum level (level 2) ####
## 1.1. Load data and make adjustments to data structure as tibble ####
phylum <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2.csv", sep = ",", header= TRUE)
str(phylum)
colnames(phylum)[1] ="SampleID"

# drop unnecessary columns
phylum <- phylum[,!(names(phylum) %in% drop_cols)]

# drop samples that have <500 feature counts
phylum_nomissing <- filter(phylum, !SampleID %in% samples_to_remove)
phylum_nomissing$SampleID <- as.factor(phylum_nomissing$SampleID)

# remove columns that sum to zero
phylum_nomissing <-
  phylum_nomissing %>% 
  select_if(~is.numeric(.) && sum(.) != 0 || !is.numeric(.))


## 1.2. Calculate total reads within a sample ####
colnames(phylum_nomissing[2]) # get first taxa name, k__Fungi.p__Basidiomycota
colnames(phylum_nomissing[ncol(phylum_nomissing)]) # get last taxa name, k__Fungi.p__Chytridiomycota

phylum_nomissing <-
  phylum_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(k__Fungi.p__Basidiomycota:k__Fungi.p__Chytridiomycota), na.rm = T))

# need to ungroup after using rowwise()
phylum_nomissing <- 
  phylum_nomissing %>%
  ungroup()

## 1.3. Calculate relative abundances to new df ####
phylum_rel_abund <-
  phylum_nomissing %>%
  mutate_at(vars(k__Fungi.p__Basidiomycota:k__Fungi.p__Chytridiomycota) , funs(relabun = ./phylum_nomissing$raw_total))

str(phylum_rel_abund)

# make sure the summation of all relabun columns sum to 1
phylum_should_equal_1 <- 
  phylum_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

phylum_should_equal_1$relabun_sum

# all = 1, so things look good!

## 1.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(phylum_rel_abund)) # column indices that have relabun in them [12:20]

phylum_rel_abund <- phylum_rel_abund[-c(2:11)] # remove everything else (but keep SampleID, which is column 1)
str(phylum_rel_abund)
colnames(phylum_rel_abund)<-gsub("_relabun","",colnames(phylum_rel_abund))

write.table(phylum_rel_abund , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun.csv", sep=",", row.names=FALSE)

## 1.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
###### - k__Fungi.__ combined with k__Fungi.p__Fungi_phy_Incertae_sedis into k__Fungi.__

phylum_rel_abund_curated <-
  phylum_rel_abund %>%
  mutate(k__Fungi.__combo = k__Fungi.__ + k__Fungi.p__Fungi_phy_Incertae_sedis)

phylum_rel_abund_curated <-
  select(phylum_rel_abund_curated, -c(k__Fungi.__, k__Fungi.p__Fungi_phy_Incertae_sedis))

## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment)

# write out file for taxa barplots
write.table(phylum_rel_abund_curated , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## No changes needed to be made for MaAsLin2, but I am just renaming it to be consistent with the other taxonomic levels
write.table(phylum_rel_abund_curated , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)































































# 2. family level (level 2) ####
## 2.1. Load data and make adjustments to data structure as tibble ####
family <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5.csv", sep = ",", header= TRUE)
str(family)
colnames(family)[1] ="SampleID"

# drop unnecessary columns
family <- family[,!(names(family) %in% drop_cols)]

# drop samples that have <500 feature counts
family_nomissing <- filter(family, !SampleID %in% samples_to_remove)
family_nomissing$SampleID <- as.factor(family_nomissing$SampleID)

# remove columns that sum to zero
family_nomissing <-
  family_nomissing %>% 
  select_if(~is.numeric(.) && sum(.) != 0 || !is.numeric(.))


## 2.2. Calculate total reads within a sample ####
colnames(family_nomissing[2]) # get first taxa name, k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae
colnames(family_nomissing[ncol(family_nomissing)]) # get last taxa name, k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae

family_nomissing <-
  family_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae), na.rm = T))

# need to ungroup after using rowwise()
family_nomissing <- 
  family_nomissing %>%
  ungroup()

## 2.3. Calculate relative abundances to new df ####
family_rel_abund <-
  family_nomissing %>%
  mutate_at(vars(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae) , funs(relabun = ./family_nomissing$raw_total))

str(family_rel_abund)

# make sure the summation of all relabun columns sum to 1
family_should_equal_1 <- 
  family_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

family_should_equal_1$relabun_sum

# all = 1, so things look good!

## 2.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(family_rel_abund)) # column indices that have relabun in them [321:638]

family_rel_abund <- family_rel_abund[-c(2:320)] # remove everything else (but keep SampleID, which is column 1)
str(family_rel_abund)
colnames(family_rel_abund)<-gsub("_relabun","",colnames(family_rel_abund))

write.table(family_rel_abund , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun.csv", sep=",", row.names=FALSE)


## 2.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
family_rel_abund_curated <-
  family_rel_abund %>%
  mutate(k__Fungi.__.__.__.__combo = k__Fungi.__.__.__.__ + k__Fungi.p__Fungi_phy_Incertae_sedis.c__Fungi_cls_Incertae_sedis.o__Fungi_ord_Incertae_sedis.f__Fungi_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.__.__.__combo = k__Fungi.p__Ascomycota.__.__.__ + k__Fungi.p__Ascomycota.c__Ascomycota_cls_Incertae_sedis.o__Ascomycota_ord_Incertae_sedis.f__Ascomycota_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Dothideomycetes.__.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.__.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideomycetes_ord_Incertae_sedis.f__Dothideomycetes_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__combo = k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__ + k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Chaetothyriales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.__combo = k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.__ + k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Helotiales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.f__Xylariales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariomycetes_ord_Incertae_sedis.f__Sordariomycetes_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Basidiomycota.c__Agaricomycetes.__.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.__.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricomycetes_ord_Incertae_sedis.f__Agaricomycetes_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Mycosphaerellales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__combo = k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__ + k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Lecanoromycetes_ord_Incertae_sedis.f__Lecanoromycetes_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Polyporales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Diaporthales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.f__Trechisporales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Cantharellales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.__combo = k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.__ + k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.f__Pezizales_fam_Incertae_sedis) %>%
  mutate(k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideales_fam_Incertae_sedis)

family_rel_abund_curated <-
  select(family_rel_abund_curated, -c(k__Fungi.__.__.__.__, k__Fungi.p__Fungi_phy_Incertae_sedis.c__Fungi_cls_Incertae_sedis.o__Fungi_ord_Incertae_sedis.f__Fungi_fam_Incertae_sedis, k__Fungi.p__Ascomycota.__.__.__, k__Fungi.p__Ascomycota.c__Ascomycota_cls_Incertae_sedis.o__Ascomycota_ord_Incertae_sedis.f__Ascomycota_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideomycetes_ord_Incertae_sedis.f__Dothideomycetes_fam_Incertae_sedis,k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Chaetothyriales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Helotiales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.f__Xylariales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariomycetes_ord_Incertae_sedis.f__Sordariomycetes_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricomycetes_ord_Incertae_sedis.f__Agaricomycetes_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Mycosphaerellales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__, k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Lecanoromycetes_ord_Incertae_sedis.f__Lecanoromycetes_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Polyporales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Diaporthales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.f__Trechisporales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Cantharellales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.__, k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.f__Pezizales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__,k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideales_fam_Incertae_sedis))

## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment)

# write out file for taxa barplots
write.table(family_rel_abund_curated , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## - FROM THIS FINAL FILE... STILL NEED SOME COLUMNS TO BE REMOVED FOR MAASLIN (ONES THAT ARE NOT IDENTIFIED AT THE LEVEL WE'RE INTERESTED IN). 

family_rel_abund_curated_MaAsLin2 <- 
  select(family_rel_abund_curated, -c(k__Fungi.__.__.__.__combo, k__Fungi.p__Ascomycota.__.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Myriangiales.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Tubeufiales.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.__.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__combo, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Phaeomoniellales.__, k__Fungi.p__Ascomycota.c__Laboulbeniomycetes.o__Pyxidiophorales.f__Pyxidiophorales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__combo, k__Fungi.p__Ascomycota.c__Leotiomycetes.__.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.__combo, k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.__combo, k__Fungi.p__Ascomycota.c__Pezizomycotina_cls_Incertae_sedis.o__Pezizomycotina_ord_Incertae_sedis.f__Pezizomycotina_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Chaetosphaeriales.f__Chaetosphaeriales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Melanosporales.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Savoryellales.f__Savoryellales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariales.f__Sordariales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sporidesmiales.f__Sporidesmiales_fam_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.__combo, k__Fungi.p__Basidiomycota.__.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.__.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Auriculariales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Corticiales.f__Corticiales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Hymenochaetales.f__Hymenochaetales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Russulales.f__Russulales_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Sebacinales.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.__combo, k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.__.__, k__Fungi.p__Basidiomycota.c__Microbotryomycetes.o__Microbotryomycetes_ord_Incertae_sedis.f__Microbotryomycetes_fam_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Tremellomycetes.__.__, k__Fungi.p__Basidiomycota.c__Tremellomycetes.o__Tremellales.__, k__Fungi.p__Blastocladiomycota.c__Blastocladiomycota_cls_Incertae_sedis.o__Blastocladiomycota_ord_Incertae_sedis.f__Blastocladiomycota_fam_Incertae_sedis, k__Fungi.p__Mucoromycota.c__Mucoromycotina_cls_Incertae_sedis.o__Mucoromycotina_ord_Incertae_sedis.f__Mucoromycotina_fam_Incertae_sedis, k__Fungi.p__Rozellomycota.c__Rozellomycota_cls_Incertae_sedis.o__Rozellomycota_ord_Incertae_sedis.f__Rozellomycota_fam_Incertae_sedis, k__Fungi.p__Rozellomycota.c__Rozellomycotina_cls_Incertae_sedis.o__GS11.f__GS11_fam_Incertae_sedis))
 
write.table(family_rel_abund_curated_MaAsLin2, file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)




































# 3. genus level (level 2) ####
## 3.1. Load data and make adjustments to data structure as tibble ####
genus <- read.csv(file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6.csv", sep = ",", header= TRUE)
str(genus)
colnames(genus)[1] ="SampleID"

# drop unnecessary columns
genus <- genus[,!(names(genus) %in% drop_cols)]

# drop samples that have <500 feature counts
genus_nomissing <- filter(genus, !SampleID %in% samples_to_remove)
genus_nomissing$SampleID <- as.factor(genus_nomissing$SampleID)

# remove columns that sum to zero
genus_nomissing <-
  genus_nomissing %>% 
  select_if(~is.numeric(.) && sum(.) != 0 || !is.numeric(.))


## 3.2. Calculate total reads within a sample ####
colnames(genus_nomissing[2]) # get first taxa name, k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium
colnames(genus_nomissing[ncol(genus_nomissing)]) # get last taxa name, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella

genus_nomissing <-
  genus_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella), na.rm = T))

# need to ungroup after using rowwise()
genus_nomissing <- 
  genus_nomissing %>%
  ungroup()

## 3.3. Calculate relative abundances to new df ####
genus_rel_abund <-
  genus_nomissing %>%
  mutate_at(vars(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella) , funs(relabun = ./genus_nomissing$raw_total))

str(genus_rel_abund)

# make sure the summation of all relabun columns sum to 1
genus_should_equal_1 <- 
  genus_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

genus_should_equal_1$relabun_sum

# all = 1, so things look good!

## 3.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(genus_rel_abund)) # column indices that have relabun in them [654:1304]

genus_rel_abund <- genus_rel_abund[-c(2:653)] # remove everything else (but keep SampleID, which is column 1)
str(genus_rel_abund)
colnames(genus_rel_abund)<-gsub("_relabun","",colnames(genus_rel_abund))

write.table(genus_rel_abund , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun.csv", sep=",", row.names=FALSE)



## 3.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
genus_rel_abund_curated <-
  genus_rel_abund %>%
  mutate(k__Fungi.__.__.__.__.___combo = k__Fungi.__.__.__.__.__ + k__Fungi.p__Fungi_phy_Incertae_sedis.c__Fungi_cls_Incertae_sedis.o__Fungi_ord_Incertae_sedis.f__Fungi_fam_Incertae_sedis.g__Fungi_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__.__combo = k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__.__ + k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Chaetothyriales_fam_Incertae_sedis.g__Chaetothyriales_gen_Incertae_sedis ) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis.g__Pleosporales_gen_Incertae_sedis ) %>%
  mutate( k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariomycetes_ord_Incertae_sedis.f__Sordariomycetes_fam_Incertae_sedis.g__Sordariomycetes_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.__ +  k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.g__Gnomoniaceae_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis.g__Capnodiales_gen_Incertae_sedis + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis.__) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideales_fam_Incertae_sedis.g__Dothideales_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Mycosphaerellales_fam_Incertae_sedis.g__Mycosphaerellales_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__.__combo = k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__.__ + k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Lecanoromycetes_ord_Incertae_sedis.f__Lecanoromycetes_fam_Incertae_sedis.g__Lecanoromycetes_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.g__Phaeosphaeriaceae_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__.__ +  k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Polyporales_fam_Incertae_sedis.g__Polyporales_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Cantharellales_fam_Incertae_sedis.g__Cantharellales_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.__combo = k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.__ + k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.g__Valsaceae_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.__combo = k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.__ + k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.g__Dothideaceae_gen_Incertae_sedis) %>%
  mutate( k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.__combo = k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.__ + k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.g__Orbiliaceae_gen_Incertae_sedis) %>%
  mutate(k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.__combo = k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.__ + k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.g__Ceratobasidiaceae_gen_Incertae_sedis)


genus_rel_abund_curated <-
  select(genus_rel_abund_curated, -c(k__Fungi.__.__.__.__.__, k__Fungi.p__Fungi_phy_Incertae_sedis.c__Fungi_cls_Incertae_sedis.o__Fungi_ord_Incertae_sedis.f__Fungi_fam_Incertae_sedis.g__Fungi_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Chaetothyriales_fam_Incertae_sedis.g__Chaetothyriales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Pleosporales_fam_Incertae_sedis.g__Pleosporales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariomycetes_ord_Incertae_sedis.f__Sordariomycetes_fam_Incertae_sedis.g__Sordariomycetes_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.g__Gnomoniaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis.g__Capnodiales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiales_fam_Incertae_sedis.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideales_fam_Incertae_sedis.g__Dothideales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Mycosphaerellales_fam_Incertae_sedis.g__Mycosphaerellales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__.__, k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Lecanoromycetes_ord_Incertae_sedis.f__Lecanoromycetes_fam_Incertae_sedis.g__Lecanoromycetes_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.g__Phaeosphaeriaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Polyporales_fam_Incertae_sedis.g__Polyporales_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Cantharellales_fam_Incertae_sedis.g__Cantharellales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.g__Valsaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.g__Dothideaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.__, k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.g__Orbiliaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.g__Ceratobasidiaceae_gen_Incertae_sedis)) 
                                      
                                      
## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment) *there were none unassigned for fungi

# write out file for taxa barplots
write.table(genus_rel_abund_curated , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## - FROM THIS FINAL FILE... STILL NEED SOME COLUMNS TO BE REMOVED FOR MAASLIN (ONES THAT ARE NOT IDENTIFIED AT THE LEVEL WE'RE INTERESTED IN). 

genus_rel_abund_curated_MaAsLin2 <-
  select(genus_rel_abund_curated, -c(k__Fungi.__.__.__.__.___combo, k__Fungi.p__Ascomycota.__.__.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.__.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Botryosphaeriales.f__Botryosphaeriaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Capnodiaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Schizothyriaceae.g__Schizothyriaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Dothideales.f__Dothideaceae.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Mycosphaerellaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Mycosphaerellales.f__Teratosphaeriaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Myriangiales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.__.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Amorosiaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Cucurbitariaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Didymellaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Didymosphaeriaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Lophiostomataceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Occultibambusaceae.g__Occultibambusaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Phaeosphaeriaceae.__combo, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Pleosporales.f__Teichosporaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Tubeufiales.__.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Tubeufiales.f__Tubeufiaceae.__, k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Venturiales.f__Sympoventuriaceae.g__Sympoventuriaceae_gen_Incertae_sedis,k__Fungi.p__Ascomycota.c__Eurotiomycetes.__.__.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.__.__combo, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Chaetothyriaceae.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Chaetothyriales.f__Herpotrichiellaceae.__, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Mycocaliciales.f__Mycocaliciaceae.g__Mycocaliciaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Eurotiomycetes.o__Phaeomoniellales.__.__, k__Fungi.p__Ascomycota.c__Laboulbeniomycetes.o__Pyxidiophorales.f__Pyxidiophorales_fam_Incertae_sedis.g__Pyxidiophorales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Lecanoromycetes.__.__.__combo, k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Caliciales.f__Physciaceae.__, k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Peltigerales.f__Pannariaceae.g__Pannariaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Lecanoromycetes.o__Teloschistales.f__Teloschistaceae.g__Teloschistaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Leotiomycetes.__.__.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.__.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Discinellaceae.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Helotiaceae.__))

genus_rel_abund_curated_MaAsLin2 <-
  select(genus_rel_abund_curated_MaAsLin2, -c(k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Hyaloscyphaceae.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Mollisiaceae.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Helotiales.f__Sclerotiniaceae.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Rhytismatales.f__Rhytismataceae.__, k__Fungi.p__Ascomycota.c__Leotiomycetes.o__Thelebolales.f__Thelebolaceae.__, k__Fungi.p__Ascomycota.c__Orbiliomycetes.o__Orbiliales.f__Orbiliaceae.__combo, k__Fungi.p__Ascomycota.c__Pezizomycetes.o__Pezizales.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.__.__.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Chaetosphaeriales.f__Chaetosphaeriaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Diaporthaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Gnomoniaceae.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Diaporthales.f__Valsaceae.__combo, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Glomerellales.f__Plectosphaerellaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Nectriaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Stachybotryaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Melanosporales.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Microascales.f__Microascaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Pleurotheciales.f__Pleurotheciaceae.g__Pleurotheciaceae_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Savoryellales.f__Savoryellales_fam_Incertae_sedis.g__Savoryellales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sordariales.f__Chaetomiaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Sporidesmiales.f__Sporidesmiales_fam_Incertae_sedis.g__Sporidesmiales_gen_Incertae_sedis, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.__.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.f__Diatrypaceae.__, k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Xylariales.f__Sporocadaceae.__, k__Fungi.p__Ascomycota.c__Taphrinomycetes.o__Taphrinales.f__Taphrinaceae.g__Taphrinaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.__.__.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.__.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricales.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricales.f__Cyphellaceae.g__Cyphellaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricales.f__Stephanosporaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Agaricales.f__Strophariaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Auriculariales.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Auriculariales.f__Aporpiaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Auriculariales.f__Exidiaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.__.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Cantharellales.f__Ceratobasidiaceae.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Corticiales.f__Corticiales_fam_Incertae_sedis.g__Corticiales_gen_Incertae_sedis))

genus_rel_abund_curated_MaAsLin2 <-
  select(genus_rel_abund_curated_MaAsLin2, -c(k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Hymenochaetales.f__Hymenochaetaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.__.__combo, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Fomitopsidaceae.g__Fomitopsidaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Irpicaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Meruliaceae.g__Meruliaceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Polyporaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Polyporales.f__Steccherinaceae.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Russulales.f__Russulales_fam_Incertae_sedis.g__Russulales_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Sebacinales.__.__, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Thelephorales.f__Thelephoraceae.g__Thelephoraceae_gen_Incertae_sedis, k__Fungi.p__Basidiomycota.c__Agaricomycetes.o__Trechisporales.__.__, k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.__.__.__, k__Fungi.p__Basidiomycota.c__Tremellomycetes.__.__.__, k__Fungi.p__Basidiomycota.c__Tremellomycetes.o__Tremellales.__.__, k__Fungi.p__Blastocladiomycota.c__Blastocladiomycota_cls_Incertae_sedis.o__Blastocladiomycota_ord_Incertae_sedis.f__Blastocladiomycota_fam_Incertae_sedis.g__Blastocladiomycota_gen_Incertae_sedis, k__Fungi.p__Mortierellomycota.c__Mortierellomycetes.o__Mortierellales.f__Mortierellaceae.__, k__Fungi.p__Mucoromycota.c__Mucoromycotina_cls_Incertae_sedis.o__Mucoromycotina_ord_Incertae_sedis.f__Mucoromycotina_fam_Incertae_sedis.g__Bifiguratus, k__Fungi.p__Rozellomycota.c__Rozellomycota_cls_Incertae_sedis.o__Rozellomycota_ord_Incertae_sedis.f__Rozellomycota_fam_Incertae_sedis.g__Rozellomycota_gen_Incertae_sedis, k__Fungi.p__Rozellomycota.c__Rozellomycotina_cls_Incertae_sedis.o__GS11.f__GS11_fam_Incertae_sedis.g__GS11_gen_Incertae_sedis))


write.table(genus_rel_abund_curated_MaAsLin2 , file = "ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)
