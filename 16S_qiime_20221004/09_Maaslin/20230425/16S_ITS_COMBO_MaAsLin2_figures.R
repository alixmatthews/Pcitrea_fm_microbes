# Differential abundance figures for manuscript
# MaAsLin2 16S and ITS
# 20230719
# updated 20230903

setwd("~/Desktop/Alix/")

library(dplyr) # v. 1.0.10
library(Maaslin2) # v. 1.14.1
library(ggplot2) # v. 3.4.0
library(ggpubr) # v. 0.5.0

# 16S PATH: 20220523_16S/AHPCC/01_qiime_20221004/09_Maaslin/20230425
# ITS PATH: 20220621_ITS/AHPCC/01_qiime_20230220/10_Maaslin/20230425

# 0. Rows/columns to remove for 16S datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove_16S <- c('PF05', 'PM24', 'PM29', 'PF14', 'PF15', 'PF04', 'PM19', 'PF19', 'PM06', 'PF28', 'PM08', 'PF25', 'PM05', 'PF20', 'PF08')

# 1. Load and Curate Sample metadata_16S ####

# This file is from the master results database used to create alpha/beta diversity metrics figures.
metadata_16S <- read.csv(file = "/home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/08_R/Master_results_20230413.csv", sep = ",", header= TRUE)
str(metadata_16S)

# Remove control samples and exclude biological samples that were not included in the alpha/beta stats (i.e., had less than 500 seqs). So, in other words, these samples have no missing data in the metadata_16S file...  "nm" = "no missing"

# drop controls
metadata_16S <- filter(metadata_16S, !bio_or_control == "control")

# drop samples that have <500 feature counts
metadata_16S <- filter(metadata_16S, !sample.id %in% samples_to_remove_16S)

# also remove 3 more samples that had zero counts at the very beginning and never made it to the alpha/beta stats
more_samples_to_remove_16S <- c('PM09', 'PF10', 'PM20')
metadata_16S <- filter(metadata_16S, !sample.id %in% more_samples_to_remove_16S)


# change variable types

metadata_16S[,c(9:10, 12:45)] <- sapply(metadata_16S[,c(9:10,12:45)],as.numeric)
# then take care of the chr->factor
metadata_16S[sapply(metadata_16S,is.character)] <- lapply(metadata_16S[sapply(metadata_16S, is.character)], as.factor)

str(metadata_16S)

# check order of sample id
metadata_16S$sample.id










# 0. Rows/columns to remove for ITS datasets ####

# Samples that have less than 500 sequences that we will drop from every level
samples_to_remove_ITS <- c('PM13', 'PM03', 'PM06', 'PF03', 'PM09', 'PM08', 'PM04', 'PF25', 'PF10', 'PF13', 'PF30')

# 1. Load and Curate Sample metadata_ITS ####

# This file is from the master results database used to create alpha/beta diversity metrics figures, but I have rearranged it to match the same order as the taxonomy file.
metadata_ITS <- read.csv(file = "/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/09_R/ITS_forward_master_results_20230413_reordered.csv", sep = ",", header= TRUE)
str(metadata_ITS)

# Remove control samples and exclude biological samples that were not included in the alpha/beta stats (i.e., had less than 500 seqs). So, in other words, these samples have no missing data in the metadata_ITS file...  "nm" = "no missing"

# drop controls
metadata_ITS <- filter(metadata_ITS, !bio_or_control == "control")

# drop samples that have <500 feature counts
metadata_ITS <- filter(metadata_ITS, !sample.id %in% samples_to_remove_ITS)

# change variable types

metadata_ITS[,c(9:10, 12:32)] <- sapply(metadata_ITS[,c(9:10, 12:32)],as.numeric)
# then take care of the chr->factor
metadata_ITS[sapply(metadata_ITS,is.character)] <- lapply(metadata_ITS[sapply(metadata_ITS, is.character)], as.factor)

str(metadata_ITS)
























# Taxonomy - 16S ####

df_rel_abun_phylum_nm_16S <- read.csv(file = "./20220523_16S/AHPCC/01_qiime_20221004/09_Maaslin/20230425/16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylum_nm_16S)

df_rel_abun_family_nm_16S <- read.csv(file = "./20220523_16S/AHPCC/01_qiime_20221004/09_Maaslin/20230425/16S_grouped-filtered-decontam-bio_samples-silva-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_family_nm_16S)

df_rel_abun_genus_nm_16S <- read.csv(file = "./20220523_16S/AHPCC/01_qiime_20221004/09_Maaslin/20230425/16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_genus_nm_16S)

# check that the sample id order is the same with metadata for cbinding
metadata_16S$sample.id
df_rel_abun_phylum_nm_16S$SampleID
df_rel_abun_family_nm_16S$SampleID
df_rel_abun_genus_nm_16S$SampleID


# let's cbind the metadata df with the rel_abun df

df_rel_abun_phylum_nm_metadata_16S <- cbind(df_rel_abun_phylum_nm_16S, metadata_16S)
which(colnames(df_rel_abun_phylum_nm_metadata_16S) == "sample.id")
df_rel_abun_phylum_nm_metadata_16S<-df_rel_abun_phylum_nm_metadata_16S[,-39] # remove duplicate sample.id column

df_rel_abun_family_nm_metadata_16S <- cbind(df_rel_abun_family_nm_16S, metadata_16S)
which(colnames(df_rel_abun_family_nm_metadata_16S) == "sample.id")
df_rel_abun_family_nm_metadata_16S<-df_rel_abun_family_nm_metadata_16S[,-276] # remove duplicate sample.id column

df_rel_abun_genus_nm_metadata_16S <- cbind(df_rel_abun_genus_nm_16S, metadata_16S)
which(colnames(df_rel_abun_genus_nm_metadata_16S) == "sample.id")
df_rel_abun_genus_nm_metadata_16S<-df_rel_abun_genus_nm_metadata_16S[,-501] # remove duplicate sample.id column




# Taxonomy - ITS ####

df_rel_abun_phylum_nm_ITS <- read.csv(file = "./20220621_ITS/AHPCC/01_qiime_20230220/10_Maaslin/20230425/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylum_nm_ITS)

df_rel_abun_family_nm_ITS <- read.csv(file = "./20220621_ITS/AHPCC/01_qiime_20230220/10_Maaslin/20230425/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_family_nm_ITS)

df_rel_abun_genus_nm_ITS <- read.csv(file = "./20220621_ITS/AHPCC/01_qiime_20230220/10_Maaslin/20230425/ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2.csv", sep = ",", header= TRUE)
str(df_rel_abun_genus_nm_ITS)


# check that the sample id order is the same with metadata for cbinding
metadata_ITS$sample.id
df_rel_abun_phylum_nm_ITS$SampleID
df_rel_abun_family_nm_ITS$SampleID
df_rel_abun_genus_nm_ITS$SampleID


# let's cbind the metadata df with the rel_abun df

df_rel_abun_phylum_nm_metadata_ITS <- cbind(df_rel_abun_phylum_nm_ITS, metadata_ITS)
which(colnames(df_rel_abun_phylum_nm_metadata_ITS) == "sample.id")
df_rel_abun_phylum_nm_metadata_ITS<-df_rel_abun_phylum_nm_metadata_ITS[,-10] # remove duplicate sample.id column

df_rel_abun_family_nm_metadata_ITS <- cbind(df_rel_abun_family_nm_ITS, metadata_ITS)
which(colnames(df_rel_abun_family_nm_metadata_ITS) == "sample.id")
df_rel_abun_family_nm_metadata_ITS<-df_rel_abun_family_nm_metadata_ITS[,-254] # remove duplicate sample.id column

df_rel_abun_genus_nm_metadata_ITS <- cbind(df_rel_abun_genus_nm_ITS, metadata_ITS)
which(colnames(df_rel_abun_genus_nm_metadata_ITS) == "sample.id")
df_rel_abun_genus_nm_metadata_ITS<-df_rel_abun_genus_nm_metadata_ITS[,-537] # remove duplicate sample.id column





















# Proteobacteria ####

relabun_Proteobacteria_jitter<-
  ggplot(df_rel_abun_phylum_nm_metadata_16S, aes(y=d__Bacteria.p__Proteobacteria, x = category_broad, fill = category_broad)) +
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


# Rhizobiaceae ####

relabun_Rhizobiaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata_16S, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae, x = category_broad, fill = category_broad)) +
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



# Rhizobiaceae (Bartonellaceae) ####

relabun_RhizobiaceaeBartonellaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata_16S, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae, x = category_broad, fill = category_broad)) +
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
  ylab("Relative Abundance of Rhizobiaceae (Bartonellaceae)") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)


relabun_RhizobiaceaeBartonellaceae_jitter



# Bartonella ####

relabun_Bartonella_jitter<-
  ggplot(df_rel_abun_genus_nm_metadata_16S, aes(y=d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Rhizobiaceae.g__Bartonella, x = category_broad, fill = category_broad)) +
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















# Basidiomycota ####

relabun_Basidiomycota_jitter<-
  ggplot(df_rel_abun_phylum_nm_metadata_ITS, aes(y=k__Fungi.p__Basidiomycota, x = category_broad, fill = category_broad)) +
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
  ylab("Relative Abundance of Basidiomycota") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "***")

relabun_Basidiomycota_jitter







# Filobasidiaceae ####

relabun_Filobasidiaceae_jitter<-
  ggplot(df_rel_abun_family_nm_metadata_ITS, aes(y=k__Fungi.p__Basidiomycota.c__Tremellomycetes.o__Filobasidiales.f__Filobasidiaceae, x = category_broad, fill = category_broad)) +
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
  ylab("Relative Abundance of Filobasidiaceae") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "***")

relabun_Filobasidiaceae_jitter





# Naganishia ####



relabun_Naganishia_jitter<-
  ggplot(df_rel_abun_genus_nm_metadata_ITS, aes(y=k__Fungi.p__Basidiomycota.c__Tremellomycetes.o__Filobasidiales.f__Filobasidiaceae.g__Naganishia, x = category_broad, fill = category_broad)) +
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
  ylab(expression(paste("Relative Abundance of ", italic("Naganishia")))) +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "***")

relabun_Naganishia_jitter




# Cladosporium for talk ####
relabun_Cladosporium_jitter<-
  ggplot(df_rel_abun_genus_nm_metadata_ITS, aes(y=k__Fungi.p__Ascomycota.c__Dothideomycetes.o__Capnodiales.f__Cladosporiaceae.g__Cladosporium, x = category_broad, fill = category_broad)) +
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
  ylab(expression(paste("Relative Abundance of ", italic("Cladosporium")))) +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=FALSE,
              textsize = 5,
              annotations = "**")

relabun_Cladosporium_jitter














# Final taxonomy figures for ms ####

relabun_Proteobacteria_jitter
relabun_Rhizobiaceae_jitter
relabun_RhizobiaceaeBartonellaceae_jitter
relabun_Bartonella_jitter
relabun_Basidiomycota_jitter
relabun_Filobasidiaceae_jitter
relabun_Naganishia_jitter


ggarrange(
  relabun_Proteobacteria_jitter,
  relabun_Rhizobiaceae_jitter,
  relabun_Bartonella_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 1)

# 4x9 landscape


ggarrange(
  relabun_Proteobacteria_jitter,
  relabun_RhizobiaceaeBartonellaceae_jitter,
  relabun_Bartonella_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 1)

# 4.5x10 landscape








ggarrange(
  relabun_Proteobacteria_jitter,
  relabun_Rhizobiaceae_jitter,
  relabun_Bartonella_jitter,
  relabun_Basidiomycota_jitter,
  relabun_Filobasidiaceae_jitter,
  relabun_Naganishia_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 2)


# 7x9 landscape





ggarrange(
  relabun_Proteobacteria_jitter,
  relabun_RhizobiaceaeBartonellaceae_jitter,
  relabun_Bartonella_jitter,
  relabun_Basidiomycota_jitter,
  relabun_Filobasidiaceae_jitter,
  relabun_Naganishia_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 2)

# 8.5x9 landscape




relabun_Proteobacteria_jitter_noxlabels <- relabun_Proteobacteria_jitter +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

relabun_Rhizobiaceae_jitter_noxlabels <- relabun_Rhizobiaceae_jitter +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

relabun_Bartonella_jitter_noxlabels <- relabun_Bartonella_jitter +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

relabun_RhizobiaceaeBartonellaceae_jitter_noxlabels <- relabun_RhizobiaceaeBartonellaceae_jitter +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())




ggarrange(
  relabun_Proteobacteria_jitter_noxlabels,
  relabun_Rhizobiaceae_jitter_noxlabels,
  relabun_Bartonella_jitter_noxlabels,
  relabun_Basidiomycota_jitter,
  relabun_Filobasidiaceae_jitter,
  relabun_Naganishia_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 2)


# 7x9 landscape




ggarrange(
  relabun_Proteobacteria_jitter_noxlabels,
  relabun_RhizobiaceaeBartonellaceae_jitter_noxlabels,
  relabun_Bartonella_jitter_noxlabels,
  relabun_Basidiomycota_jitter,
  relabun_Filobasidiaceae_jitter,
  relabun_Naganishia_jitter,
  align = "hv",
  labels = "AUTO",
  ncol = 3,
  nrow = 2)

# 8.5x9 landscape








relabun_Proteobacteria_jitter_noxlabels_yaxis <-
  relabun_Proteobacteria_jitter_noxlabels +
  theme(
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=12)
    )

relabun_Bartonella_jitter_noxlabels_yaxis <-
  relabun_Bartonella_jitter_noxlabels +
  theme(
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=12)
  )

relabun_Basidiomycota_jitter_yaxis <-
  relabun_Basidiomycota_jitter +
  theme(
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=12)
  )


relabun_Naganishia_jitter_yaxis <-
  relabun_Naganishia_jitter +
  theme(
    axis.title.y = element_text(size=12),
    axis.text.y = element_text(size=12)
  )




ggarrange(
  relabun_Proteobacteria_jitter_noxlabels_yaxis,
  relabun_Bartonella_jitter_noxlabels_yaxis,
  relabun_Basidiomycota_jitter_yaxis,
  relabun_Naganishia_jitter_yaxis,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 2,
  nrow = 2)

# 8x8 portrait
