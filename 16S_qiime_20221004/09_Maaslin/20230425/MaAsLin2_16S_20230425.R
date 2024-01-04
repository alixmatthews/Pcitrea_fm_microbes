#### MaAsLin2 for 16S
#### Alix Matthews
#### 25 April 2023

### R version 4.2.3 (2023-03-15) "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

#### These .csv files were taken from Qiime2 output and ran through the RelativeAbundanceCalculations_$GENE_$DATE.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for MaAsLin2

setwd("~/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/09_Maaslin/20230425")

library(dplyr) # v. 1.0.10
library(Maaslin2) # v. 1.10.1
library(ggplot2) # v. 3.4.0
library(ggpubr) # v. 0.5.0

# 0. Rows/columns to remove for all datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove <- c('PF05', 'PM24', 'PM29', 'PF14', 'PF15', 'PF04', 'PM19', 'PF19', 'PM06', 'PF28', 'PM08', 'PF25', 'PM05', 'PF20', 'PF08')

# 1. Load and Curate Sample Metadata ####

# This file is from the master results database used to create alpha/beta diversity metrics figures.
metadata <- read.csv(file = "/home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/08_R/Master_results_20230413.csv", sep = ",", header= TRUE)
str(metadata)

# Remove control samples and exclude biological samples that were not included in the alpha/beta stats (i.e., had less than 500 seqs). So, in other words, these samples have no missing data in the metadata file...  "nm" = "no missing"

# drop controls
metadata <- filter(metadata, !bio_or_control == "control")

# drop samples that have <500 feature counts
metadata <- filter(metadata, !sample.id %in% samples_to_remove)

# also remove 3 more samples that had zero counts at the very beginning and never made it to the alpha/beta stats
more_samples_to_remove <- c('PM09', 'PF10', 'PM20')
metadata <- filter(metadata, !sample.id %in% more_samples_to_remove)


# change variable types

metadata[,c(9:10, 12:45)] <- sapply(metadata[,c(9:10,12:45)],as.numeric)
# then take care of the chr->factor
metadata[sapply(metadata,is.character)] <- lapply(metadata[sapply(metadata, is.character)], as.factor)

str(metadata)




# 2. Load Data ####

#### These data include only biological samples that were included in the alpha/beta diversity stats (no missing samples = nm) and have undergone some manual curation specifically for MaAsLin (removing taxa that are unidentified at the taxa level of interest, combining taxa, etc... all is outlined in the RelativeAbundanceCalculations_$GENE_$DATE.R script)


## 2.1. Phylum ####

df_rel_abun_phylum_nm <- read.csv(file = "16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylum_nm)


## 2.2. Family ####

df_rel_abun_family_nm <- read.csv(file = "16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_family_nm)


## 2.3. Genus ####

df_rel_abun_genus_nm <- read.csv(file = "16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_genus_nm)






# 3. Run MaAsLin2 - Phylum ####

## 3.1. Age and Sex ####
#### Explanation: Try testing for effect of bird age and bird sex on differential microbial abundance. If these are NS, then we can exclude them from the downstream models.

fit_data_cplm_none_agesex_phylum <-
  Maaslin2(
    input_data = df_rel_abun_phylum_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_agesex_phylum", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("bird_sex,bird_age"),
    random_effects = c("bird_id")
  )

#### Results: Age and sex appear to be NS. Now let's just test for differential microbial abundance between the groups (feathers versus mites) while controlling for bird_id (we will still need to control for this as mite/feather samples from the same bird are not independent)


## 3.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_phylum <-
  Maaslin2(
    input_data = df_rel_abun_phylum_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_phylum", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )

#### Results: There are some taxa that are differential abundant between feathers and mites. (**** = > 10% and move forward with)
# d__Bacteria.p__Actinobacteriota ****
# d__Bacteria.p__Cyanobacteria
# d__Bacteria.p__Deinococcota
# d__Bacteria.p__Armatimonadota
# d__Bacteria.p__Proteobacteria ****
# d__Bacteria.p__Bdellovibrionota
# d__Bacteria.p__Gemmatimonadota
# d__Bacteria.p__Bacteroidota ****
# d__Bacteria.__combo



## 3.3 Plotting significant results ####

# let's cbind the metadata df with the rel_abun df

df_rel_abun_phylum_nm_metadata <- cbind(df_rel_abun_phylum_nm, metadata)
names(df_rel_abun_phylum_nm_metadata)
df_rel_abun_phylum_nm_metadata<-df_rel_abun_phylum_nm_metadata[,-39] # remove duplicate sample.id column


#### 3.3.1 Actinobacteriota ####
relabun_Actinobacteriota<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Actinobacteriota, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Actinobacteriota") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

relabun_Actinobacteriota

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Actinobacteriota.pdf", device=cairo_pdf)


relabun_Actinobacteriota_jitter<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Actinobacteriota, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Actinobacteriota") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

relabun_Actinobacteriota_jitter

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Actinobacteriota_jitter.pdf", device=cairo_pdf)


#### 3.3.2 Proteobacteria ####
relabun_Proteobacteria<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Proteobacteria") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Proteobacteria

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Proteobacteria.pdf", device=cairo_pdf)



relabun_Proteobacteria_jitter<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Proteobacteria") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Proteobacteria_jitter

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Proteobacteria_jitter.pdf", device=cairo_pdf)




#### 3.3.3 Bacteroidota ####
relabun_Bacteroidota<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Bacteroidota, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Bacteroidota") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Bacteroidota

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Bacteroidota.pdf", device=cairo_pdf)



relabun_Bacteroidota_jitter<-
  ggplot(df_rel_abun_phylum_nm_metadata, aes(y=d__Bacteria.p__Bacteroidota, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Bacteroidota") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Bacteroidota_jitter

ggsave("./fit_data_cplm_none_category_phylum/figures/relabun_Bacteroidota_jitter.pdf", device=cairo_pdf)










#### 3.3.4. Put figures all together ####
ggarrange(
  relabun_Actinobacteriota_jitter,
  relabun_Bacteroidota_jitter,
  relabun_Proteobacteria_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 1)

# saved as 3.5 x 10 landscape




























# 4. Run MaAsLin2 - Family ####

## 4.1. Age and Sex ####
#### Explanation: Try testing for effect of bird age and bird sex on differential microbial abundance. If these are NS, then we can exclude them from the downstream models.

fit_data_cplm_none_agesex_family <-
  Maaslin2(
    input_data = df_rel_abun_family_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_agesex_family", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("bird_sex,bird_age"),
    random_effects = c("bird_id")
  )

#### Results: Age and sex appear to be NS. Now let's just test for differential microbial abundance between the groups (feathers versus mites) while controlling for bird_id (we will still need to control for this as mite/feather samples from the same bird are not independent)




## 4.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_family <-
  Maaslin2(
    input_data = df_rel_abun_family_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_family", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )

#### Results: There are some taxa that are differential abundant between feathers and mites (**** = > 10% and move forward with)

# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Hymenobacteraceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae
# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Chitinophagales.f__Chitinophagaceae
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae **** 
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__Geodermatophilaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Nocardioidaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Micrococcales.f__Intrasporangiaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Pseudonocardiales.f__Pseudonocardiaceae
# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Spirosomaceae
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Propionibacteriaceae **** 
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Steroidobacterales.f__Steroidobacteraceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Oxalobacteraceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Salinisphaerales.f__Solimonadaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Kineosporiales.f__Kineosporiaceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae **** 
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__Frankiaceae
# d__Bacteria.p__Verrucomicrobiota.c__Verrucomicrobiae.o__Opitutales.f__Puniceicoccaceae
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Acetobacterales.f__Acetobacteraceae
# d__Bacteria.p__Firmicutes.c__Bacilli.o__Paenibacillales.f__Paenibacillaceae
# d__Bacteria.p__Gemmatimonadota.c__Gemmatimonadetes.o__Gemmatimonadales.f__Gemmatimonadaceae
# d__Bacteria.p__Firmicutes.c__Clostridia.o__Peptostreptococcales.Tissierellales.f__Peptostreptococcaceae
# d__Bacteria.p__Planctomycetota.c__Planctomycetes.o__Isosphaerales.f__Isosphaeraceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales.f__Vibrionaceae


## 4.3 Plotting significant results ####

# let's cbind the metadata df with the rel_abun df

df_rel_abun_family_nm_metadata <- cbind(df_rel_abun_family_nm, metadata)
names(df_rel_abun_family_nm_metadata)
df_rel_abun_family_nm_metadata<-df_rel_abun_family_nm_metadata[,-276] # remove duplicate sample.id column

# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Propionibacteriaceae
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae


#### 4.3.1 Rhizobiaceae ####
relabun_Rhizobiaceae<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Rhizobiaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

relabun_Rhizobiaceae

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Rhizobiaceae.pdf", device=cairo_pdf)



relabun_Rhizobiaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Rhizobiaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)


relabun_Rhizobiaceae_jitter

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Rhizobiaceae_jitter.pdf", device=cairo_pdf)




#### 4.3.2 Propionibacteriaceae ####
relabun_Propionibacteriaceae<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Propionibacteriaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Propionibacteriaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

relabun_Propionibacteriaceae

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Propionibacteriaceae.pdf", device=cairo_pdf)



relabun_Propionibacteriaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Propionibacteriaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Propionibacteriaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

relabun_Propionibacteriaceae_jitter

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Propionibacteriaceae_jitter.pdf", device=cairo_pdf)




#### 4.3.3 Burkholderiaceae ####
# not > 10%...
relabun_Burkholderiaceae<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Burkholderiaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Burkholderiaceae

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Burkholderiaceae.pdf", device=cairo_pdf)



relabun_Burkholderiaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Burkholderiaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Burkholderiaceae_jitter

ggsave("./fit_data_cplm_none_category_family/figures/relabun_Burkholderiaceae_jitter.pdf", device=cairo_pdf)







#### 4.3.4. Put figures all together ####
ggarrange(
  relabun_Propionibacteriaceae_jitter,
  relabun_Rhizobiaceae_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 2,
  nrow = 1)

# saved as 3.75 x 7 landscape






















# 5. Run MaAsLin2 - Genus ####

## 5.1. Age and Sex ####
#### Explanation: Try testing for effect of bird age and bird sex on differential microbial abundance. If these are NS, then we can exclude them from the downstream models.

fit_data_cplm_none_agesex_genus <-
  Maaslin2(
    input_data = df_rel_abun_genus_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_agesex_genus", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("bird_sex,bird_age"),
    random_effects = c("bird_id")
  )

#### Results: Age and sex appear to be NS. Now let's just test for differential microbial abundance between the groups (feathers versus mites) while controlling for bird_id (we will still need to control for this as mite/feather samples from the same bird are not independent)



## 5.2. Feathers vs. Mites ####
#### Explanation: Try testing for effect of sample type (feathers vs. mites) on differential microbial abundance

fit_data_cplm_none_category_genus <-
  Maaslin2(
    input_data = df_rel_abun_genus_nm, 
    input_metadata = metadata, 
    output = "fit_data_cplm_none_category_genus", 
    analysis_method = "CPLM",
    normalization = "NONE",
    transform = "NONE",
    fixed_effects = c("category_broad"),
    random_effects = c("bird_id")
  )

#### Results: genera that was differential abundant between feathers and mites as significant at the q-value (**** = > 10% and move forward with)


# d__Bacteria.p__Gemmatimonadota.c__Gemmatimonadetes.o__Gemmatimonadales.f__Gemmatimonadaceae.g__Roseisolibacter
# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Cytophagales.f__Hymenobacteraceae.g__Hymenobacter
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Methylophilaceae.g__Methylotenera
# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Sphingobacteriales.f__Sphingobacteriaceae.g__Solitalea
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Novosphingobium
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Propionibacteriales.f__Nocardioidaceae.g__Nocardioides
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Pseudonocardiales.f__Pseudonocardiaceae.g__Actinomycetospora
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Bartonella ****
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Delftia
# d__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Flavobacteriaceae.g__Flavobacterium
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.g__Methylobacterium.Methylorubrum
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Beijerinckiaceae.g__1174.901.12
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Salinisphaerales.f__Solimonadaceae.g__Nevskia
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Oxalobacteraceae.g__Massilia
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Caulobacterales.f__Caulobacteraceae.g__Phenylobacterium
# d__Bacteria.p__Actinobacteriota.c__Actinobacteria.o__Frankiales.f__Frankiaceae.g__Jatrophihabitans
# d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Sphingomonas ****
# d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae.g__Limnobacter



## 5.3 Plotting significant results ####

# let's cbind the metadata df with the rel_abun df

df_rel_abun_genus_nm_metadata <- cbind(df_rel_abun_genus_nm, metadata)
which(colnames(df_rel_abun_genus_nm_metadata) == "sample.id")
df_rel_abun_genus_nm_metadata<-df_rel_abun_genus_nm_metadata[,-501] # remove duplicate sample.id column





#### 5.3.1 Bartonella ####
relabun_Bartonella<-
  ggplot(df_rel_abun_genus_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Bartonella, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab(expression(paste("Relative Abundance of ", italic("Bartonella")))) +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "***")

relabun_Bartonella

ggsave("./fit_data_cplm_none_category_genus/figures/relabun_Bartonella.pdf", device=cairo_pdf)



relabun_Bartonella_jitter<-
  ggplot(df_rel_abun_genus_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Bartonella, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab(expression(paste("Relative Abundance of ", italic("Bartonella")))) +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "***")

relabun_Bartonella_jitter

ggsave("./fit_data_cplm_none_category_genus/figures/relabun_Bartonella_jitter.pdf", device=cairo_pdf)





#### 5.3.2 Sphingomonas ####
relabun_Sphingomonas<-
  ggplot(df_rel_abun_genus_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Sphingomonas, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  # geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Sphingomonas") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Sphingomonas

ggsave("./fit_data_cplm_none_category_genus/figures/relabun_Sphingomonas.pdf", device=cairo_pdf)



relabun_Sphingomonas_jitter<-
  ggplot(df_rel_abun_genus_nm_metadata, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae.g__Sphingomonas, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Relative Abundance of Sphingomonas") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "*")

relabun_Sphingomonas_jitter

ggsave("./fit_data_cplm_none_category_genus/figures/relabun_Sphingomonas_jitter.pdf", device=cairo_pdf)






#### 5.3.3. Put figures all together ####
ggarrange(
  relabun_Bartonella_jitter,
  relabun_Sphingomonas_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 2,
  nrow = 1)

# saved as 3.75 x 7 landscape


















# 6. Bartonella taxonomy figures for ms ####

relabun_Proteobacteria_jitter
relabun_Rhizobiaceae_jitter
relabun_Bartonella_jitter


ggarrange(
  relabun_Proteobacteria_jitter,
  relabun_Rhizobiaceae_jitter,
  relabun_Bartonella_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 1)

# 4x9 landscape