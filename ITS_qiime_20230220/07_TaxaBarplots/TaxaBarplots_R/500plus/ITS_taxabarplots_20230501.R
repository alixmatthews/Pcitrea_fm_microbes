#### Taxa barplots for ITS
#### Alix Matthews
#### 2023 May 1

#### These .csv files were taken from Qiime2 output and ran through the RelativeAbundanceCalculations_$GENE_$DATE.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for making these taxa barplots

### R version 4.2.3 (2023-03-15) "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

# 1. Set working directory, install packges, load libraries ####
setwd("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/07_TaxaBarplots/R/500plus")

library(dplyr) # 1.0.10
library(ggplot2) # 3.4.0 
library(data.table) # 1.14.6
library(randomcoloR) # 1.1.0.1
library(RColorBrewer) # 1.1-3
library(plotrix) # 3.8-2

# 2. Phylum level ####

df_rel_abun_phylum <- read.csv(file = "./ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylum)

df_rel_abun_phylum_long <- melt(setDT(df_rel_abun_phylum), id.vars = c("SampleID"), variable.name = c("OTU_ID_Phylum"))

ggplot(df_rel_abun_phylum_long, aes(x = SampleID, y = value, fill = OTU_ID_Phylum)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))


## 2.1 Collapse phyla with very low relative abundances ####

# Well this is fun and all, but now I want to get collapse all those ASVs that are in really tiny proportions; they're not really telling us what's going on

df_rel_abun_phylum_long$OTU_ID_Phylum_condensed_1per <- df_rel_abun_phylum_long$OTU_ID_Phylum # copy OTU_ID column to a new column so to not lose the original, just in case


df_rel_abun_phylum_long_OTU_ID_condensed <- 
  df_rel_abun_phylum_long %>%
  mutate(OTU_ID_Phylum_condensed_1per = ifelse(value <= 0.01, '< 1%', as.character(OTU_ID_Phylum_condensed_1per)))




## Take a quick look
ggplot(df_rel_abun_phylum_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))




## 2.2 Change colors etc. ####

# Let's change colors and factor names
# First let's get a sorted comma separated list of the names we need, with double quotes
sort(unique(df_rel_abun_phylum_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per <- sort(unique(df_rel_abun_phylum_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per_rs <- paste("\"",as.character(list_phyla_1per),"\"",collapse=", ",sep="")
cat(list_phyla_1per_rs)

# new names; find and replace magic of output from `cat(list..)`
phylum_new_labels <- c("< 1%", "Unidentifed fungi", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Mortierellomycota")

phylum_color <- c("gray40","ivory3",
                          "plum4",
                          "slategray2",
                          "violetred4",
                          "lightpink")
                          
ggplot(df_rel_abun_phylum_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color)


# Let's try to facet this plot
# First make a new column with a factor as either mite or feather
# View(df_rel_abun_phylum_long_OTU_ID_condensed)

df_rel_abun_phylum_long_OTU_ID_condensed$Category_broad <- df_rel_abun_phylum_long_OTU_ID_condensed$SampleID # copy sampleID column to a new column so to not lose the original, just in case

# convert category broad to feather or mites for plotting purposes
df_rel_abun_phylum_long_OTU_ID_condensed <- df_rel_abun_phylum_long_OTU_ID_condensed %>%  
  mutate(Category_broad = case_when(grepl("PM", Category_broad) ~ "Mites",
                                    grepl("PF", Category_broad, ignore.case = TRUE) ~"Feather"))

## 2.3. Final plots ####

ggplot(df_rel_abun_phylum_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_wrap(~Category_broad, scales = "free") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  theme(axis.text.x=element_text(angle=75,hjust=1)) +
  theme(strip.text.x = element_text(size = 18)) +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.y = element_text(size=12))

# saved as 8 x 13 landscape


taxa_ITS_final<-
  ggplot(df_rel_abun_phylum_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  theme(legend.title=element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  theme(strip.text.x = element_text(size = 26)) +
  theme(axis.title.y = element_text(size=22)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=16)) +
  guides(fill=guide_legend(ncol=1))

# saved as 8 x 13 landscape




### 2.3.1. Do some summarizing of the data for manuscript #####
str(df_rel_abun_phylum_long_OTU_ID_condensed)

# look at everything
print(n = 100, df_rel_abun_phylum_long_OTU_ID_condensed %>%
        group_by(Category_broad) %>%
        filter(value != 0) %>%
        distinct(OTU_ID_Phylum))

# summarize # of phyla per category while filtering out certain 'OTUs' that aren't really phyla
df_rel_abun_phylum_long_OTU_ID_condensed %>%
  group_by(Category_broad) %>%
  filter(value != 0) %>%
  filter(OTU_ID_Phylum != "k__Fungi.__combo") %>%
  summarise(OTU_ID_Phylum=n_distinct(OTU_ID_Phylum))

# get the average and SE of the values per category to determine the most common phyla OTUs
print(n = 100, df_rel_abun_phylum_long_OTU_ID_condensed %>%
        group_by(Category_broad,OTU_ID_Phylum) %>%
        filter(value != 0) %>%
        filter(OTU_ID_Phylum != "Unassigned.__") %>%
        filter(OTU_ID_Phylum != "d__Bacteria.__") %>%
        summarize(Average = mean(value), SE = std.error(value)) %>%
        ungroup() %>%
        arrange(desc(Average)))




