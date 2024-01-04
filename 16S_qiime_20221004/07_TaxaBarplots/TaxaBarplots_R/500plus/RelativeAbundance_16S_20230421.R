#### Calculate 16S relative abundance at phylum, family, and genus for import for taxa barplots and MaAsLin2
#### Alix Matthews
#### 21 April 2023

#### These .csv files were taken from Qiime2 output `taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva.qzv` and downloading the .csv file for the specific level of interest.

setwd("~/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/07_TaxaBarplots/R/500plus")

library(dplyr) # v. 1.0.10

# 0. Rows/columns to remove for all datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove <- c('PF05', 'PM24', 'PM29', 'PF14', 'PF15', 'PF04', 'PM19', 'PF19', 'PM06', 'PF28', 'PM08', 'PF25', 'PM05', 'PF20', 'PF08')

# Unncessary columns to be dropped from every level
drop_cols <- c("bio_or_control", "category_broad", "category_specific", "bird_id", "bird_sex", "bird_age", "capture_date", "num_mites_fm", "num_mites_m", "mite_sp_fm", "mite_sp_m")


# 1. Phylum level (level 2) ####
## 1.1. Load data and make adjustments to data structure as tibble ####
phylum <- read.csv(file = "taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-2.csv", sep = ",", header= TRUE)
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
colnames(phylum_nomissing[2]) # get first taxa name, d__Bacteria.p__Armatimonadota
colnames(phylum_nomissing[ncol(phylum_nomissing)]) # get last taxa name, d__Bacteria.p__RCP2.54

phylum_nomissing <-
  phylum_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(d__Bacteria.p__Armatimonadota:d__Bacteria.p__RCP2.54), na.rm = T))

# need to ungroup after using rowwise()
phylum_nomissing <- 
  phylum_nomissing %>%
  ungroup()



## 1.3. Calculate relative abundances to new df ####
phylum_rel_abund <-
  phylum_nomissing %>%
  mutate_at(vars(d__Bacteria.p__Armatimonadota:d__Bacteria.p__RCP2.54) , funs(relabun = ./phylum_nomissing$raw_total))

str(phylum_rel_abund)

# make sure the summation of all relabun columns sum to 1
phylum_should_equal_1 <- 
  phylum_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

phylum_should_equal_1$relabun_sum

# all = 1, so things look good!


## 1.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(phylum_rel_abund)) # column indices that have relabun in them [41:78]

phylum_rel_abund <- phylum_rel_abund[-c(2:40)] # remove everything else (but keep SampleID, which is column 1)
str(phylum_rel_abund)
colnames(phylum_rel_abund)<-gsub("_relabun","",colnames(phylum_rel_abund))

write.table(phylum_rel_abund , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun.csv", sep=",", row.names=FALSE)


## 1.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
###### - d__Bacteria.__ combined with d__Bacteria.p__uncultured into d__Bacteria.__combo

phylum_rel_abund_curated <-
  phylum_rel_abund %>%
  mutate(d__Bacteria.__combo = d__Bacteria.__ + d__Bacteria.p__uncultured)

phylum_rel_abund_curated <-
  select(phylum_rel_abund_curated, -c(d__Bacteria.__, d__Bacteria.p__uncultured))

## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment)

# write out file for taxa barplots
write.table(phylum_rel_abund_curated , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## No changes needed to be made for MaAsLin2, but I am just renaming it to be consistent with the other taxonomic levels
write.table(phylum_rel_abund_curated , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)













# 2. Family level (level 5) ####
## 2.1. Load data and make adjustments to data structure as tibble ####
family <- read.csv(file = "taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-5.csv", sep = ",", header= TRUE)
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
colnames(family_nomissing[2]) # get first taxa name, d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae
colnames(family_nomissing[ncol(family_nomissing)]) # get last taxa name, d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108

family_nomissing <-
  family_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae:d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108), na.rm = T))

# need to ungroup after using rowwise()
family_nomissing <- 
  family_nomissing %>%
  ungroup()

## 2.3. Calculate relative abundances to new df ####
family_rel_abund <-
  family_nomissing %>%
  mutate_at(vars(d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae:d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108) , funs(relabun = ./family_nomissing$raw_total))

str(family_rel_abund)

# make sure the summation of all relabun columns sum to 1
family_should_equal_1 <- 
  family_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

family_should_equal_1$relabun_sum

# all = 1, so things look good!


## 2.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(family_rel_abund)) # column indices that have relabun in them [377:750]

family_rel_abund <- family_rel_abund[-c(2:376)] # remove everything else (but keep SampleID, which is column 1)
str(family_rel_abund)
colnames(family_rel_abund)<-gsub("_relabun","",colnames(family_rel_abund))

write.table(family_rel_abund , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun.csv", sep=",", row.names=FALSE)

## 2.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
family_rel_abund_curated <-
  family_rel_abund %>%
  mutate(d__Bacteria.__.__.__.__combo = d__Bacteria.__.__.__.__ + d__Bacteria.p__uncultured.c__uncultured.o__uncultured.f__uncultured) %>%
  mutate(d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__combo = d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__ + d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__combo = d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__ + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__uncultured.f__uncultured + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria_Incertae_Sedis.f__Unknown_Family) %>%
  mutate(d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__combo = d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__ + d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__uncultured.f__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__uncultured.f__uncultured)

family_rel_abund_curated <-
  select(family_rel_abund_curated, -c(d__Bacteria.__.__.__.__, d__Bacteria.p__uncultured.c__uncultured.o__uncultured.f__uncultured, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__uncultured.f__uncultured,d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria_Incertae_Sedis.f__Unknown_Family,d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__uncultured.f__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__uncultured.f__uncultured))

## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment)

# write out file for taxa barplots
write.table(family_rel_abund_curated , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## - FROM THIS FINAL FILE... STILL NEED SOME COLUMNS TO BE REMOVED FOR MAASLIN (ONES THAT ARE NOT IDENTIFIED AT THE LEVEL WE'RE INTERESTED IN). 

family_rel_abund_curated_MaAsLin2 <- 
  select(family_rel_abund_curated, -c(d__Bacteria.__.__.__.__combo , d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__combo , d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__combo, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__combo , d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__combo, d__Bacteria.p__Proteobacteria.__.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiales_Incertae_Sedis, d__Bacteria.p__Acidobacteriota.c__Acidobacteriae.o__Acidobacteriales.f__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.__, d__Bacteria.p__Acidobacteriota.c__Vicinamibacteria.o__Vicinamibacterales.f__uncultured, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Corynebacteriales.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.__, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.__.__, d__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.Selenomonadales.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales.__, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.__, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Planctomycetales.f__uncultured, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Gaiellales.f__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.__, d__Bacteria.p__Firmicutes.c__Bacilli.__.__, d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Micavibrionales.f__uncultured, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Micrococcales_Incertae_Sedis, d__Bacteria.p__Desulfobacterota.c__uncultured.o__uncultured.f__uncultured, d__Bacteria.p__Patescibacteria.c__Parcubacteria.__.__, d__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.__, d__Bacteria.p__Actinobacteriota.c__Acidimicrobiia.o__Microtrichales.f__uncultured, d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.__.__, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__uncultured.f__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__uncultured, d__Bacteria.p__Bdellovibrionota.c__Oligoflexia.o__Oligoflexales.f__uncultured, d__Bacteria.p__Cyanobacteria.c__Cyanobacteriia.o__Oxyphotobacteria_Incertae_Sedis.f__Unknown_Family, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.__))

# remove the ones that are not taxonomically informative at the family level (i.e., same name at order as family)... in two steps because of how many need to be removed...
family_rel_abund_curated_MaAsLin2 <- 
  select(family_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Acidobacteriota.c__Acidobacteriae.o__Subgroup_13.f__Subgroup_13, d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__11.24.f__11.24, d__Bacteria.p__Acidobacteriota.c__Holophagae.o__Subgroup_7.f__Subgroup_7, d__Bacteria.p__Acidobacteriota.c__Subgroup_18.o__Subgroup_18.f__Subgroup_18, d__Bacteria.p__Acidobacteriota.c__Subgroup_22.o__Subgroup_22.f__Subgroup_22, d__Bacteria.p__Acidobacteriota.c__Vicinamibacteria.o__Subgroup_17.f__Subgroup_17, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__PeM15.f__PeM15, d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108, d__Bacteria.p__Bdellovibrionota.c__Oligoflexia.o__0319.6G20.f__0319.6G20, d__Bacteria.p__Caldatribacteriota.c__JS1.o__JS1.f__JS1, d__Bacteria.p__Chloroflexi.c__Anaerolineae.o__SBR1031.f__SBR1031, d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.o__FW22.f__FW22, d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.o__Napoli.4B.65.f__Napoli.4B.65, d__Bacteria.p__Chloroflexi.c__Gitt.GS.136.o__Gitt.GS.136.f__Gitt.GS.136, d__Bacteria.p__Chloroflexi.c__JG30.KF.CM66.o__JG30.KF.CM66.f__JG30.KF.CM66, d__Bacteria.p__Chloroflexi.c__KD4.96.o__KD4.96.f__KD4.96, d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__C0119.f__C0119, d__Bacteria.p__Chloroflexi.c__TK10.o__TK10.f__TK10, d__Bacteria.p__Elusimicrobiota.c__Elusimicrobia.o__Lineage_IV.f__Lineage_IV, d__Bacteria.p__Elusimicrobiota.c__Lineage_IIb.o__Lineage_IIb.f__Lineage_IIb, d__Bacteria.p__Firmicutes.c__Bacilli.o__RF39.f__RF39, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia_UCG.014.f__Clostridia_UCG.014, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia_vadinBB60_group.f__Clostridia_vadinBB60_group, d__Bacteria.p__Firmicutes.c__Clostridia.o__Peptostreptococcales.Tissierellales.f__Peptostreptococcales.Tissierellales, d__Bacteria.p__Firmicutes.c__TSAC18.o__TSAC18.f__TSAC18, d__Bacteria.p__Myxococcota.c__bacteriap25.o__bacteriap25.f__bacteriap25, d__Bacteria.p__Myxococcota.c__Polyangia.o__Blfdi19.f__Blfdi19, d__Bacteria.p__Myxococcota.c__Polyangia.o__mle1.27.f__mle1.27, d__Bacteria.p__Nitrospirota.c__4.29.1.o__4.29.1.f__4.29.1, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Buchananbacteria.f__Candidatus_Buchananbacteria, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Falkowbacteria.f__Candidatus_Falkowbacteria, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Magasanikbacteria.f__Candidatus_Magasanikbacteria, d__Bacteria.p__Patescibacteria.c__CPR2.o__CPR2.f__CPR2, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Absconditabacteriales_.SR1..f__Absconditabacteriales_.SR1., d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Candidatus_Peregrinibacteria.f__Candidatus_Peregrinibacteria, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Candidatus_Peribacteria.f__Candidatus_Peribacteria, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Gracilibacteria.f__Gracilibacteria, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__JGI_0000069.P22.f__JGI_0000069.P22, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Azambacteria.f__Candidatus_Azambacteria, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Brennerbacteria.f__Candidatus_Brennerbacteria, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Giovannonibacteria.f__Candidatus_Giovannonibacteria, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Kaiserbacteria.f__Candidatus_Kaiserbacteria,  d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Nomurabacteria.f__Candidatus_Nomurabacteria, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Yanofskybacteria.f__Candidatus_Yanofskybacteria, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__GWB1.42.6.f__GWB1.42.6, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Parcubacteria.f__Parcubacteria, d__Bacteria.p__Planctomycetota.c__OM190.o__OM190.f__OM190, d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__CCM11a.f__CCM11a, d__Bacteria.p__Planctomycetota.c__vadinHA49.o__vadinHA49.f__vadinHA49))

# keep removing...
family_rel_abund_curated_MaAsLin2 <- 
  select(family_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__JG36.TzT.191.f__JG36.TzT.191, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__R7C24.f__R7C24, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__WD260.f__WD260, d__Bacteria.p__RCP2.54.c__RCP2.54.o__RCP2.54.f__RCP2.54, d__Bacteria.p__Rs.K70_termite_group.c__Rs.K70_termite_group.o__Rs.K70_termite_group.f__Rs.K70_termite_group, d__Bacteria.p__Verrucomicrobiota.c__Omnitrophia.o__Omnitrophales.f__Omnitrophales, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__S.BQ2.57_soil_group.f__S.BQ2.57_soil_group, d__Bacteria.p__WPS.2.c__WPS.2.o__WPS.2.f__WPS.2))

write.table(family_rel_abund_curated_MaAsLin2, file = "16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)






# 3. Genus level (level 6) ####
## 3.1. Load data and make adjustments to data structure as tibble ####
genus <- read.csv(file = "taxa_barplots-16S_grouped-filtered-decontam-bio_samples-silva-level-6.csv", sep = ",", header= TRUE)
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
colnames(genus_nomissing[2]) # get first taxa name, d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae.g__Fimbriimonadaceae
colnames(genus_nomissing[ncol(genus_nomissing)]) # get last taxa name, d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108.g__MB.A2.108

genus_nomissing <-
  genus_nomissing %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae.g__Fimbriimonadaceae:d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108.g__MB.A2.108), na.rm = T))

# need to ungroup after using rowwise()
genus_nomissing <- 
  genus_nomissing %>%
  ungroup()

## 3.3. Calculate relative abundances to new df ####
genus_rel_abund <-
  genus_nomissing %>%
  mutate_at(vars(d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae.g__Fimbriimonadaceae:d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108.g__MB.A2.108) , funs(relabun = ./genus_nomissing$raw_total))

str(genus_rel_abund)

# make sure the summation of all relabun columns sum to 1
genus_should_equal_1 <- 
  genus_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

genus_should_equal_1$relabun_sum

# all = 1, so things look good!


## 3.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(genus_rel_abund)) # column indices that have relabun in them [763:1522]

genus_rel_abund <- genus_rel_abund[-c(2:762)] # remove everything else (but keep SampleID, which is column 1)
str(genus_rel_abund)
colnames(genus_rel_abund)<-gsub("_relabun","",colnames(genus_rel_abund))

write.table(genus_rel_abund , file = "16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun.csv", sep=",", row.names=FALSE)



## 3.5. Curate relative abundance file for taxa barplots and MaAsLin2 ####
## This file needs to be curated manually prior to taxa barplots and MaAsLin

## - Combine certain columns.
genus_rel_abund_curated <-
  genus_rel_abund %>%
  mutate(d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.__combo = d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.__ + d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.__combo = d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.__ + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.__combo = d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.__ + d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__uncultured) %>%
  mutate(d__Bacteria.__.__.__.__.__combo = d__Bacteria.__.__.__.__.__ + d__Bacteria.p__uncultured.c__uncultured.o__uncultured.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__.__combo = d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__.__ + d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.__combo = d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.__ + d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Patescibacteria.c__Gracilibacteria.__.__.__combo = d__Bacteria.p__Patescibacteria.c__Gracilibacteria.__.__.__ + d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Gracilibacteria.f__Gracilibacteria.g__Gracilibacteria) %>%
  mutate(d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.__combo = d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.__ + d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.__ +  d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.__combo =  d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.__ + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.__combo = d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.__ + d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.__combo = d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.__ + d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__.__combo = d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__.__ + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__uncultured.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.__combo = d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.__ + d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.__combo = d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.__ + d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.__combo = d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.__ + d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.__combo = d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.__ + d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.__combo = d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.__ + d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.__combo = d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.__ + d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__.__combo = d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__.__ + d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__uncultured.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.__combo = d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.__ + d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.__combo = d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.__ + d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.__combo = d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.__ + d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.__combo = d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.__ + d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.__combo = d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.__ + d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__.__combo = d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__.__ + d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__uncultured.f__uncultured.g__uncultured) %>%
  mutate(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.__combo = d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.__ + d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.__combo = d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.__ + d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.g__Incertae_Sedis + d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.g__uncultured) %>%
  mutate(d__Bacteria.p__Patescibacteria.c__Parcubacteria.__.__.__combo = d__Bacteria.p__Patescibacteria.c__Parcubacteria.__.__.__ + d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Parcubacteria.f__Parcubacteria.g__Parcubacteria) %>%
  mutate(d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__uncultured_combo = d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__uncultured + d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__Spirosomaceae)


# need to run this in a few batches
genus_rel_abund_curated <-
  select(genus_rel_abund_curated, -c(d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.__, d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.g__uncultured, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__uncultured, d__Bacteria.__.__.__.__.__, d__Bacteria.p__uncultured.c__uncultured.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__uncultured.g__uncultured, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.__, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.__.__.__, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Gracilibacteria.f__Gracilibacteria.g__Gracilibacteria, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.g__uncultured))

genus_rel_abund_curated <-
  select(genus_rel_abund_curated, -c(d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.__, d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__uncultured, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.__, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.g__uncultured, d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.__, d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__uncultured, d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.__, d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.g__uncultured, d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.__, d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.g__uncultured, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.__, d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.g__uncultured, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.__, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.g__uncultured, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.__, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.g__uncultured, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.__, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.g__uncultured))

genus_rel_abund_curated <-
  select(genus_rel_abund_curated, -c(d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.__, d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.g__uncultured, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.__, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.g__Incertae_Sedis, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__Parcubacteria.__.__.__, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Parcubacteria.f__Parcubacteria.g__Parcubacteria, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__Spirosomaceae))


## - Keep 'unassigned' as is (even though it could be a mix of different taxa within this unassigned assignment)

# write out file for taxa barplots
write.table(genus_rel_abund_curated, file = "16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated.csv", sep=",", row.names=FALSE)

## - FROM THIS FINAL FILE... STILL NEED SOME COLUMNS TO BE REMOVED FOR MAASLIN (ONES THAT ARE NOT IDENTIFIED AT THE LEVEL WE'RE INTERESTED IN). 
#### have to do this in multiple rounds because there are so many...


genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated, -c(d__Bacteria.__.__.__.__.__combo, d__Bacteria.p__Firmicutes.c__Clostridia.o__Lachnospirales.f__Lachnospiraceae.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Salinisphaerales.f__Solimonadaceae.__, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Candidatus_Peribacteria.f__Candidatus_Peribacteria.g__Candidatus_Peribacteria, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Falkowbacteria.f__Candidatus_Falkowbacteria.g__Candidatus_Falkowbacteria, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.__combo, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Steroidobacterales.f__Steroidobacteraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.__combo,  d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.__.__combo, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Oscillospiraceae.__combo, d__Bacteria.p__Bdellovibrionota.c__Bdellovibrionia.o__Bacteriovoracales.f__Bacteriovoracaceae.__, d__Bacteria.p__Acidobacteriota.c__Acidobacteriae.o__Acidobacteriales.f__uncultured.g__uncultured, d__Bacteria.p__Chloroflexi.c__Gitt.GS.136.o__Gitt.GS.136.f__Gitt.GS.136.g__Gitt.GS.136, d__Bacteria.p__Margulisbacteria.c__Margulisbacteria.o__Margulisbacteria.f__Margulisbacteria.g__Margulisbacteria, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia_vadinBB60_group.f__Clostridia_vadinBB60_group.g__Clostridia_vadinBB60_group))


genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Microbacteriaceae.__, d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__11.24.f__11.24.g__11.24, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.__.__.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__SC.I.84.g__SC.I.84, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Magasanikbacteria.f__Candidatus_Magasanikbacteria.g__Candidatus_Magasanikbacteria, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Pedosphaerales.f__Pedosphaeraceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Barnesiellaceae.g__uncultured, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Verrucomicrobiales.f__Verrucomicrobiaceae.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodobacterales.f__Rhodobacteraceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.__combo, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__Sporichthyaceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__Candidatus_Peregrinibacteria.f__Candidatus_Peregrinibacteria.g__Candidatus_Peregrinibacteria, d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__Arcobacteraceae.__combo, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Methylacidiphilales.f__Methylacidiphilaceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Nomurabacteria.f__Candidatus_Nomurabacteria.g__Candidatus_Nomurabacteria, d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Sandaracinaceae.g__uncultured, d__Bacteria.p__Firmicutes.c__Bacilli.o__RF39.f__RF39.g__RF39, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Rhodocyclaceae.__, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__67.14.g__67.14, d__Bacteria.p__Patescibacteria.c__Saccharimonadia.o__Saccharimonadales.f__LWQ8.g__LWQ8, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__GWB1.42.6.f__GWB1.42.6.g__GWB1.42.6, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Streptomycetales.f__Streptomycetaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Sphingobacteriales.f__NS11.12_marine_group.g__NS11.12_marine_group, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Saprospiraceae.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__uncultured, d__Bacteria.p__Acidobacteriota.c__Vicinamibacteria.o__Vicinamibacterales.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Sutterellaceae.g__uncultured, d__Bacteria.p__Chloroflexi.c__JG30.KF.CM66.o__JG30.KF.CM66.f__JG30.KF.CM66.g__JG30.KF.CM66, d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__Tepidisphaerales.f__WD2101_soil_group.g__WD2101_soil_group, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Xanthobacteraceae.__combo, d__Bacteria.p__Nitrospirota.c__4.29.1.o__4.29.1.f__4.29.1.g__4.29.1, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Nocardioidaceae.__, d__Bacteria.p__Bdellovibrionota.c__Oligoflexia.o__0319.6G20.f__0319.6G20.g__0319.6G20, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Corynebacteriales.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Hyphomonadaceae.__combo))

genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Synergistota.c__Synergistia.o__Synergistales.f__Synergistaceae.__, d__Bacteria.p__Patescibacteria.c__Gracilibacteria.o__JGI_0000069.P22.f__JGI_0000069.P22.g__JGI_0000069.P22, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__BD2.7.g__BD2.7, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae.__combo, d__Bacteria.p__Verrucomicrobiota.c__Chlamydiae.o__Chlamydiales.f__Parachlamydiaceae.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.__combo, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Paracaedibacterales.f__Paracaedibacteraceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__ABY1.o__Candidatus_Buchananbacteria.f__Candidatus_Buchananbacteria.g__Candidatus_Buchananbacteria, d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__CCM11a.f__CCM11a.g__CCM11a, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.__.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylomonadaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Oxalobacteraceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__WD260.f__WD260.g__WD260, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Erwiniaceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Cellvibrionales.f__Cellvibrionaceae.g__uncultured, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Clostridiaceae.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Neisseriaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__NS9_marine_group.g__NS9_marine_group, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Azambacteria.f__Candidatus_Azambacteria.g__Candidatus_Azambacteria, d__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.Selenomonadales.__.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__PeM15.f__PeM15.g__PeM15, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__.Eubacterium._coprostanoligenes_group.g__.Eubacterium._coprostanoligenes_group, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Gemmatales.f__Gemmataceae.__combo, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales.__.__, d__Bacteria.p__Cyanobacteria.c__Vampirivibrionia.o__Gastranaerophilales.f__Gastranaerophilales.g__Gastranaerophilales, d__Bacteria.p__Chloroflexi.c__KD4.96.o__KD4.96.f__KD4.96.g__KD4.96, d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__Polyangiaceae.__combo, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.__.__, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__A21b.g__A21b, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Planctomycetales.f__uncultured.g__uncultured, d__Bacteria.p__Chloroflexi.c__Chloroflexia.o__Chloroflexales.f__Roseiflexaceae.g__uncultured, d__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Bacillaceae.__))



genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Intrasporangiaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Sphingobacteriales.f__Sphingobacteriaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.__combo, d__Bacteria.p__Planctomycetota.c__vadinHA49.o__vadinHA49.f__vadinHA49.g__vadinHA49, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Gaiellales.f__uncultured.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__uncultured.g__uncultured, d__Bacteria.p__Cyanobacteria.c__Vampirivibrionia.o__Vampirovibrionales.f__Vampirovibrionaceae.g__Vampirovibrionaceae, d__Bacteria.p__Caldatribacteriota.c__JS1.o__JS1.f__JS1.g__JS1, d__Bacteria.p__Acidobacteriota.c__Blastocatellia.o__Blastocatellales.f__Blastocatellaceae.__combo, d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Geobacterales.f__Geobacteraceae.__combo, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.__.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Cryomorphaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__JG36.TzT.191.f__JG36.TzT.191.g__JG36.TzT.191, d__Bacteria.p__Firmicutes.c__Bacilli.__.__.__, d__Bacteria.p__Myxococcota.c__Polyangia.o__Polyangiales.f__BIrii41.g__BIrii41, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__SM2D12.g__SM2D12, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiales_Incertae_Sedis.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.__.__.__combo, d__Bacteria.p__Campilobacterota.c__Campylobacteria.o__Campylobacterales.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__R7C24.f__R7C24.g__R7C24, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Chthoniobacterales.f__Chthoniobacteraceae.__, d__Bacteria.p__Myxococcota.c__Myxococcia.o__Myxococcales.f__Myxococcaceae.__combo, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae.__combo, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__Solirubrobacterales.f__Solirubrobacteraceae.__combo, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia.f__Hungateiclostridiaceae.__combo, d__Bacteria.p__Chloroflexi.c__Anaerolineae.o__SBR1031.f__A4b.g__A4b, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Planctomycetales.f__Rubinisphaeraceae.g__uncultured, d__Bacteria.p__Planctomycetota.c__OM190.o__OM190.f__OM190.g__OM190, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__UCG.010.g__UCG.010, d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__Ktedonobacterales.f__Ktedonobacteraceae.__combo, d__Bacteria.p__Firmicutes.c__Clostridia.o__Peptococcales.f__Peptococcaceae.g__uncultured))

genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.__.__combo, d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__Phycisphaerales.f__Phycisphaeraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Micavibrionales.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Micropepsales.f__Micropepsaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Amoebophilaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.__.__.__combo, d__Bacteria.p__Rs.K70_termite_group.c__Rs.K70_termite_group.o__Rs.K70_termite_group.f__Rs.K70_termite_group.g__Rs.K70_termite_group, d__Bacteria.p__Acidobacteriota.c__Subgroup_18.o__Subgroup_18.f__Subgroup_18.g__Subgroup_18, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__TRA3.20.g__TRA3.20, d__Bacteria.p__Dependentiae.c__Babeliae.o__Babeliales.f__Babeliales.g__Babeliales, d__Bacteria.p__Acidobacteriota.c__Acidobacteriae.o__Subgroup_13.f__Subgroup_13.g__Subgroup_13, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Weeksellaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Sphingobacteriales.f__env.OPS_17.g__env.OPS_17, d__Bacteria.p__Acidobacteriota.c__Subgroup_22.o__Subgroup_22.f__Subgroup_22.g__Subgroup_22, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Diplorickettsiales.f__Diplorickettsiaceae.__combo, d__Bacteria.p__Chloroflexi.c__Chloroflexia.o__Thermomicrobiales.f__JG30.KF.CM45.g__JG30.KF.CM45, d__Bacteria.p__Firmicutes.c__Clostridia.o__Oscillospirales.f__Ruminococcaceae.__combo, d__Bacteria.p__Chloroflexi.c__Anaerolineae.o__SBR1031.f__SBR1031.g__SBR1031, d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridia_UCG.014.f__Clostridia_UCG.014.g__Clostridia_UCG.014, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Brennerbacteria.f__Candidatus_Brennerbacteria.g__Candidatus_Brennerbacteria, d__Bacteria.p__Firmicutes.c__Syntrophomonadia.o__Syntrophomonadales.f__Syntrophomonadaceae.__, d__Bacteria.p__Deinococcota.c__Deinococci.o__Deinococcales.f__Deinococcaceae.g__uncultured, d__Bacteria.p__Desulfobacterota.c__uncultured.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Planococcaceae.__, d__Bacteria.p__Elusimicrobiota.c__Lineage_IIb.o__Lineage_IIb.f__Lineage_IIb.g__Lineage_IIb, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Saccharospirillaceae.g__uncultured, d__Bacteria.p__Patescibacteria.c__Parcubacteria.__.__.__combo, d__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.__.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.__, d__Bacteria.p__Gemmatimonadota.c__Gemmatimonadetes.o__Gemmatimonadales.f__Gemmatimonadaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__SAR11_clade.f__Clade_II.g__Clade_II, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Legionellales.f__Legionellaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__Rickettsiaceae.g__uncultured, d__Bacteria.p__Chloroflexi.c__TK10.o__TK10.f__TK10.g__TK10, d__Bacteria.p__Actinobacteriota.c__Acidimicrobiia.o__Microtrichales.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__T34.g__T34, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Microscillaceae.g__uncultured))


genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.__.__.__, d__Bacteria.p__Acidobacteriota.c__Holophagae.o__Subgroup_7.f__Subgroup_7.g__Subgroup_7, d__Bacteria.p__Firmicutes.c__Clostridia.o__Christensenellales.f__Christensenellaceae.g__uncultured, d__Bacteria.p__Actinobacteriota.c__Acidimicrobiia.o__Microtrichales.f__Microtrichaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae.g__uncultured_combo, d__Bacteria.p__Patescibacteria.c__CPR2.o__CPR2.f__CPR2.g__CPR2, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Giovannonibacteria.f__Candidatus_Giovannonibacteria.g__Candidatus_Giovannonibacteria, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Rhodanobacteraceae.__, d__Bacteria.p__Bdellovibrionota.c__Oligoflexia.o__Silvanigrellales.f__Silvanigrellaceae.g__uncultured, d__Bacteria.p__Elusimicrobiota.c__Elusimicrobia.o__Lineage_IV.f__Lineage_IV.g__Lineage_IV, d__Bacteria.p__Actinobacteriota.c__Thermoleophilia.o__uncultured.f__uncultured.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Morganellaceae.g__endosymbionts, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Devosiaceae.__, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Yanofskybacteria.f__Candidatus_Yanofskybacteria.g__Candidatus_Yanofskybacteria, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae.g__uncultured, d__Bacteria.p__Chloroflexi.c__Ktedonobacteria.o__C0119.f__C0119.g__C0119, d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Pirellulales.f__Pirellulaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__uncultured.g__uncultured, d__Bacteria.p__Bdellovibrionota.c__Oligoflexia.o__Oligoflexales.f__uncultured.g__uncultured, d__Bacteria.p__WPS.2.c__WPS.2.o__WPS.2.f__WPS.2.g__WPS.2, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Methyloligellaceae.__, d__Bacteria.p__Patescibacteria.c__Parcubacteria.o__Candidatus_Kaiserbacteria.f__Candidatus_Kaiserbacteria.g__Candidatus_Kaiserbacteria, d__Bacteria.p__Myxococcota.c__Polyangia.o__mle1.27.f__mle1.27.g__mle1.27))


genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Myxococcota.c__bacteriap25.o__bacteriap25.f__bacteriap25.g__bacteriap25, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Opitutales.f__Puniceicoccaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.__.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Opitutales.f__Opitutaceae.__, d__Bacteria.p__Acidobacteriota.c__Acidobacteriae.o__Acidobacteriales.f__Acidobacteriaceae_.Subgroup_1..__, d__Bacteria.p__Planctomycetota.c__Phycisphaerae.o__Tepidisphaerales.f__CPla.3_termite_group.g__CPla.3_termite_group, d__Bacteria.p__Chloroflexi.c__Anaerolineae.o__Caldilineales.f__Caldilineaceae.g__uncultured, d__Bacteria.p__Acidobacteriota.c__Vicinamibacteria.o__Subgroup_17.f__Subgroup_17.g__Subgroup_17, d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.o__FW22.f__FW22.g__FW22, d__Bacteria.p__Myxococcota.c__Polyangia.o__Blfdi19.f__Blfdi19.g__Blfdi19, d__Bacteria.p__Desulfobacterota.c__Desulfovibrionia.o__Desulfovibrionales.f__Desulfovibrionaceae.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Dermacoccaceae.__, d__Bacteria.p__Firmicutes.c__TSAC18.o__TSAC18.f__TSAC18.g__TSAC18, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.__, d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__S.BQ2.57_soil_group.f__S.BQ2.57_soil_group.g__S.BQ2.57_soil_group, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__Rhodospirillaceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Hydrogenophilaceae.g__uncultured, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Rs.E47_termite_group.g__Rs.E47_termite_group, d__Bacteria.p__Fibrobacterota.c__Fibrobacteria.o__Fibrobacterales.f__Fibrobacteraceae.g__uncultured, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__Magnetospirillaceae.g__uncultured, d__Bacteria.p__RCP2.54.c__RCP2.54.o__RCP2.54.f__RCP2.54.g__RCP2.54, d__Bacteria.p__Actinobacteriota.c__MB.A2.108.o__MB.A2.108.f__MB.A2.108.g__MB.A2.108, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Cellulomonadaceae.__, d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Dysgonomonadaceae.__, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micromonosporales.f__Micromonosporaceae.__))

# remove the ones that are not taxonomically informative at the family level (i.e., same name at order as family)... in two steps because of how many need to be removed...
genus_rel_abund_curated_MaAsLin2 <- 
  select(genus_rel_abund_curated_MaAsLin2, -c(d__Bacteria.p__Bacteroidota.c__Kapabacteria.o__Kapabacteriales.f__Kapabacteriales.g__Kapabacteriales, d__Bacteria.p__Chloroflexi.c__Dehalococcoidia.o__Napoli.4B.65.f__Napoli.4B.65.g__Napoli.4B.65, d__Bacteria.p__Cyanobacteria.c__Sericytochromatia.o__Sericytochromatia.f__Sericytochromatia.g__Sericytochromatia, d__Bacteria.p__Dependentiae.c__Babeliae.o__Babeliales.f__UBA12409.g__UBA12409, d__Bacteria.p__Dependentiae.c__Babeliae.o__Babeliales.f__Vermiphilaceae.g__Vermiphilaceae, d__Bacteria.p__Desulfobacterota.c__Desulfuromonadia.o__Bradymonadales.f__Bradymonadales.g__Bradymonadales, d__Bacteria.p__Patescibacteria.c__Saccharimonadia.o__Saccharimonadales.f__Saccharimonadales.g__Saccharimonadales, d__Bacteria.p__Proteobacteria.__.__.__.__, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rickettsiales.f__AB1.g__AB1, d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Salinisphaerales.f__Solimonadaceae.g__uncultured, d__Bacteria.p__Verrucomicrobiota.c__Omnitrophia.o__Omnitrophales.f__Omnitrophales.g__Omnitrophales, d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Dermacoccaceae.g__Dermacoccaceae, d__Bacteria.p__Acidobacteriota.c__Vicinamibacteria.o__Vicinamibacterales.f__Vicinamibacteraceae.g__Vicinamibacteraceae, d__Bacteria.p__Armatimonadota.c__Fimbriimonadia.o__Fimbriimonadales.f__Fimbriimonadaceae.g__Fimbriimonadaceae,d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Muribaculaceae.g__Muribaculaceae,d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Rikenellaceae, d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.g__Beijerinckiaceae))

write.table(genus_rel_abund_curated_MaAsLin2, file = "16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep=",", row.names=FALSE)

