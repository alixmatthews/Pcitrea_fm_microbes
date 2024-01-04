#### Taxa barplots for 16S and ITS combos
#### Alix Matthews
#### 2023 Sept 3
#### updated 2023 Dec 6 (reorder bars by prevalence)

#### These .csv files were taken from Qiime2 output and ran through the RelativeAbundanceCalculations_$GENE_$DATE.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for making these taxa barplots

### R version 4.3.1 (2023-06-16)

### OS: Ubuntu 20.04.3 LTS

# 1. Set working directory, install packges, load libraries ####
setwd("/home/lse305/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/07_TaxaBarplots/R/500plus")

library(dplyr) # 1.1.2
library(ggplot2) # 3.4.2 
library(data.table) # 1.14.8
library(randomcoloR) # 1.1.0.1
library(RColorBrewer) # 1.1-3
library(plotrix) # 3.8-2
library(ggpubr) # 0.6.0


# 2. Phylum level ####

df_rel_abun_phylum16s <- read.csv(file = "./16S_grouped-filtered-decontam-bio_samples-silva-level-2_nomissing_relabun_curated.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylum16s)

df_rel_abun_phylum16s_long <- melt(setDT(df_rel_abun_phylum16s), id.vars = c("SampleID"), variable.name = c("OTU_ID_Phylum"))

ggplot(df_rel_abun_phylum16s_long, aes(x = SampleID, y = value, fill = OTU_ID_Phylum)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))


## 2.1 Collapse phyla with very low relative abundances ####

# Well this is fun and all, but now I want to get collapse all those ASVs that are in really tiny proportions; they're not really telling us what's going on

df_rel_abun_phylum16s_long$OTU_ID_Phylum_condensed_1per <- df_rel_abun_phylum16s_long$OTU_ID_Phylum # copy OTU_ID column to a new column so to not lose the original, just in case

df_rel_abun_phylum16s_long_OTU_ID_condensed <- 
  df_rel_abun_phylum16s_long %>%
  mutate(OTU_ID_Phylum_condensed_1per = ifelse(value <= 0.01, '< 1%', as.character(OTU_ID_Phylum_condensed_1per)))

## Take a quick look
ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))


## 2.2 Change colors etc. ####

# Let's change colors and factor names
# First let's get a sorted comma separated list of the names we need, with double quotes
sort(unique(df_rel_abun_phylum16s_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per <- sort(unique(df_rel_abun_phylum16s_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per_rs <- paste("\"",as.character(list_phyla_1per),"\"",collapse=", ",sep="")
cat(list_phyla_1per_rs)

phylum_new_labels <- c("< 1%", "Unidentified bacteria", "Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota", "Caldatribacteriota", "Campilobacterota", "Chloroflexi", "Cloacimonadota", "Cyanobacteria", "Deferribacterota", "Deinococcota", "Desulfobacterota", "Elusimicrobiota", "Firmicutes", "Fusobacteriota", "Myxococcota", "Patescibacteria", "Planctomycetota", "Proteobacteria", "Synergistota", "Thermotogota", "Verrucomicrobiota", "Unassigned")

phylum_color <- c("gray40","lightcyan2",
                  "brown3","navajowhite3","khaki3","olivedrab4","skyblue4","violetred4",
                  "gray60",
                  "tomato3","sandybrown","lightgoldenrod2","palegreen3","lightblue3","plum4",
                  "gray80",
                  "orangered4","chocolate2","darkgoldenrod1","forestgreen","darkslateblue","paleturquoise4",
                  "ivory3", 
                  "hotpink3", "lightpink", "bisque3")



ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color)


# Let's try to facet this plot
# First make a new column with a factor as either mite or feather
# View(df_rel_abun_phylum16s_long_OTU_ID_condensed)

df_rel_abun_phylum16s_long_OTU_ID_condensed$Category_broad <- df_rel_abun_phylum16s_long_OTU_ID_condensed$SampleID # copy sampleID column to a new column so to not lose the original, just in case

# convert category broad to feather or mites for plotting purposes
df_rel_abun_phylum16s_long_OTU_ID_condensed <- df_rel_abun_phylum16s_long_OTU_ID_condensed %>%  
  mutate(Category_broad = case_when(grepl("PM", Category_broad) ~ "Mites",
                                    grepl("PF", Category_broad, ignore.case = TRUE) ~"Feather"))

## 2.3. Final plots ####

ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
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


taxa_16S_final<-
  ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scales = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_16S_final

# saved as 8 x 13 landscape


## order by decreasing abundance of Proteobacteria (most common)

taxa_16S_feather_ordered <-
  df_rel_abun_phylum16s_long_OTU_ID_condensed %>%
  filter(Category_broad == "Feather") %>%
  filter(OTU_ID_Phylum == "d__Bacteria.p__Proteobacteria") %>%
  arrange(desc(value)) %>%
  mutate(order = 1:nrow(.))

taxa_16S_mites_ordered <-
  df_rel_abun_phylum16s_long_OTU_ID_condensed %>%
  filter(Category_broad == "Mites") %>%
  filter(OTU_ID_Phylum == "d__Bacteria.p__Proteobacteria") %>%
  arrange(desc(value)) %>%
  mutate(order = 1:nrow(.))


# make copy of original dataset because that is a good idea

df_rel_abun_phylum16s_long_OTU_ID_condensedOrdered <- df_rel_abun_phylum16s_long_OTU_ID_condensed

df_rel_abun_phylum16s_long_OTU_ID_condensedOrdered$SampleIDOrdered <- factor(df_rel_abun_phylum16s_long_OTU_ID_condensedOrdered$SampleID, levels = c("PF11", "PF13", "PF18", "PF26", "PF03", "PF22","PF21","PF01","PF17","PF12","PF23","PF16","PF09","PF07","PF30","PF29","PF24","PF27","PF02","PF06","PM02","PM16","PM27","PM11","PM01","PM18","PM10","PM07","PM26","PM14","PM13","PM22","PM12","PM21","PM03","PM23","PM17","PM30","PM28","PM25","PM04"))


taxa_16S_finalOrdered<-
  ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensedOrdered, aes(x = SampleIDOrdered, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scale = "free_x", space = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_16S_finalOrdered


# saved as 8 x 13 landscape




taxa_16S_finalOrdered2<-
  ggplot(df_rel_abun_phylum16s_long_OTU_ID_condensedOrdered, aes(x = SampleIDOrdered, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scale = "free_x", space = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_16S_finalOrdered2


# saved as 8 x 13 landscape















### 2.3.1. Do some summarizing of the data for manuscript #####
str(df_rel_abun_phylum16s_long_OTU_ID_condensed)

# look at everything
print(n = 100, df_rel_abun_phylum16s_long_OTU_ID_condensed %>%
        group_by(Category_broad) %>%
        filter(value != 0) %>%
        distinct(OTU_ID_Phylum))

# summarize # of phyla per category while filtering out certain 'OTUs' that aren't really phyla
df_rel_abun_phylum16s_long_OTU_ID_condensed %>%
  group_by(Category_broad) %>%
  filter(value != 0) %>%
  filter(OTU_ID_Phylum != "Unassigned.__") %>%
  filter(OTU_ID_Phylum != "d__Bacteria.__combo") %>%
  summarise(OTU_ID_Phylum=n_distinct(OTU_ID_Phylum))

# get the average and SE of the values per category to determine the most common phyla OTUs
print(n = 100, df_rel_abun_phylum16s_long_OTU_ID_condensed %>%
        group_by(Category_broad,OTU_ID_Phylum) %>%
        filter(value != 0) %>%
        filter(OTU_ID_Phylum != "Unassigned.__") %>%
        filter(OTU_ID_Phylum != "d__Bacteria.__") %>%
        summarize(Average = mean(value), SE = std.error(value)) %>%
        ungroup() %>%
        arrange(desc(Average)))



































#### Taxa barplots for ITS

#### These .csv files were taken from Qiime2 output and ran through the RelativeAbundanceCalculations_$GENE_$DATE.R script, further curated the dataset (as is outlined in that script), and now they're officially ready for import for making these taxa barplots


# 1. Set working directory, install packges, load libraries ####
setwd("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/07_TaxaBarplots/R/500plus")

# 2. Phylum level ####

df_rel_abun_phylumITS <- read.csv(file = "./ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated.csv", sep = ",", header= TRUE)
str(df_rel_abun_phylumITS)

df_rel_abun_phylumITS_long <- melt(setDT(df_rel_abun_phylumITS), id.vars = c("SampleID"), variable.name = c("OTU_ID_Phylum"))

ggplot(df_rel_abun_phylumITS_long, aes(x = SampleID, y = value, fill = OTU_ID_Phylum)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))


## 2.1 Collapse phyla with very low relative abundances ####

# Well this is fun and all, but now I want to get collapse all those ASVs that are in really tiny proportions; they're not really telling us what's going on

df_rel_abun_phylumITS_long$OTU_ID_Phylum_condensed_1per <- df_rel_abun_phylumITS_long$OTU_ID_Phylum # copy OTU_ID column to a new column so to not lose the original, just in case


df_rel_abun_phylumITS_long_OTU_ID_condensed <- 
  df_rel_abun_phylumITS_long %>%
  mutate(OTU_ID_Phylum_condensed_1per = ifelse(value <= 0.01, '< 1%', as.character(OTU_ID_Phylum_condensed_1per)))




## Take a quick look
ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))




## 2.2 Change colors etc. ####

# Let's change colors and factor names
# First let's get a sorted comma separated list of the names we need, with double quotes
sort(unique(df_rel_abun_phylumITS_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per <- sort(unique(df_rel_abun_phylumITS_long_OTU_ID_condensed$OTU_ID_Phylum_condensed_1per))
list_phyla_1per_rs <- paste("\"",as.character(list_phyla_1per),"\"",collapse=", ",sep="")
cat(list_phyla_1per_rs)

# new names; find and replace magic of output from `cat(list..)`
phylum_new_labels <- c("< 1%", "Unidentifed fungi", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Mortierellomycota")

phylum_color <- c("gray40","ivory3",
                  "plum4",
                  "slategray2",
                  "violetred4",
                  "lightpink")

ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color)


# Let's try to facet this plot
# First make a new column with a factor as either mite or feather
# View(df_rel_abun_phylumITS_long_OTU_ID_condensed)

df_rel_abun_phylumITS_long_OTU_ID_condensed$Category_broad <- df_rel_abun_phylumITS_long_OTU_ID_condensed$SampleID # copy sampleID column to a new column so to not lose the original, just in case

# convert category broad to feather or mites for plotting purposes
df_rel_abun_phylumITS_long_OTU_ID_condensed <- df_rel_abun_phylumITS_long_OTU_ID_condensed %>%  
  mutate(Category_broad = case_when(grepl("PM", Category_broad) ~ "Mites",
                                    grepl("PF", Category_broad, ignore.case = TRUE) ~"Feather"))

## 2.3. Final plots ####

ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
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
  ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensed, aes(x = SampleID, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scales = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_ITS_final

# saved as 8 x 13 landscape



## order by decreasing abundance of Proteobacteria (most common)

taxa_ITS_feather_ordered <-
  df_rel_abun_phylumITS_long_OTU_ID_condensed %>%
  filter(Category_broad == "Feather") %>%
  filter(OTU_ID_Phylum == "k__Fungi.p__Ascomycota") %>%
  arrange(desc(value)) %>%
  mutate(order = 1:nrow(.))

taxa_ITS_mites_ordered <-
  df_rel_abun_phylumITS_long_OTU_ID_condensed %>%
  filter(Category_broad == "Mites") %>%
  filter(OTU_ID_Phylum == "k__Fungi.p__Ascomycota") %>%
  arrange(desc(value)) %>%
  mutate(order = 1:nrow(.))


# make copy of original dataset because that is a good idea

df_rel_abun_phylumITS_long_OTU_ID_condensedOrdered <- df_rel_abun_phylumITS_long_OTU_ID_condensed

df_rel_abun_phylumITS_long_OTU_ID_condensedOrdered$SampleIDOrdered <- factor(df_rel_abun_phylumITS_long_OTU_ID_condensedOrdered$SampleID, levels = c("PF15", "PF27", "PF09", "PF16", "PF14", "PF19","PF21","PF24","PF20","PF02","PF29","PF23","PF04","PF26","PF11","PF28","PF06","PF12","PF05","PF07","PF08", "PF18", "PF22", "PF01", "PF17", "PM30","PM25","PM21","PM17","PM19","PM29","PM16","PM14","PM20","PM05","PM10","PM12","PM23","PM02","PM22","PM27","PM18","PM26","PM01","PM11","PM24","PM07","PM28"))


taxa_ITS_finalOrdered<-
  ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensedOrdered, aes(x = SampleIDOrdered, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scale = "free_x", space = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_ITS_finalOrdered


# saved as 8 x 13 landscape




taxa_ITS_finalOrdered2<-
  ggplot(df_rel_abun_phylumITS_long_OTU_ID_condensedOrdered, aes(x = SampleIDOrdered, y = value, fill = OTU_ID_Phylum_condensed_1per)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(labels = phylum_new_labels, values = phylum_color) +
  ylab("Relative Abundance") +
  xlab("Sample ID") +
  facet_grid(~Category_broad, scale = "free_x", space = "free", switch = "both") +
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
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
  theme(strip.background = element_rect(fill = "gray90"),
        strip.text.x = element_text(margin = margin(0,0,0,0, "pt")),
        strip.text.x.bottom = element_text(size = 20))

taxa_ITS_finalOrdered2


# saved as 8 x 13 landscape







### 2.3.1. Do some summarizing of the data for manuscript #####
str(df_rel_abun_phylumITS_long_OTU_ID_condensed)

# look at everything
print(n = 100, df_rel_abun_phylumITS_long_OTU_ID_condensed %>%
        group_by(Category_broad) %>%
        filter(value != 0) %>%
        distinct(OTU_ID_Phylum))

# summarize # of phyla per category while filtering out certain 'OTUs' that aren't really phyla
df_rel_abun_phylumITS_long_OTU_ID_condensed %>%
  group_by(Category_broad) %>%
  filter(value != 0) %>%
  filter(OTU_ID_Phylum != "k__Fungi.__combo") %>%
  summarise(OTU_ID_Phylum=n_distinct(OTU_ID_Phylum))

# get the average and SE of the values per category to determine the most common phyla OTUs
print(n = 100, df_rel_abun_phylumITS_long_OTU_ID_condensed %>%
        group_by(Category_broad,OTU_ID_Phylum) %>%
        filter(value != 0) %>%
        filter(OTU_ID_Phylum != "Unassigned.__") %>%
        filter(OTU_ID_Phylum != "d__Bacteria.__") %>%
        summarize(Average = mean(value), SE = std.error(value)) %>%
        ungroup() %>%
        arrange(desc(Average)))



























# Combine taxa barplots for 16S and ITS ####


taxa_16S_final_combo <- taxa_16S_final
taxa_ITS_final_combo <- taxa_ITS_final



ggarrange(
  taxa_16S_final_combo,
  taxa_ITS_final_combo,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 1,
  nrow = 2)

# save as 15x10 portrait




# ordered barcharts

taxa_16S_final_comboOrdered <- taxa_16S_finalOrdered
taxa_ITS_final_comboOrdered <- taxa_ITS_finalOrdered



ggarrange(
  taxa_16S_final_comboOrdered,
  taxa_ITS_final_comboOrdered,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 1,
  nrow = 2)

# save as 15x10 portrait





# ordered barcharts with white space

taxa_16S_final_comboOrdered2 <- taxa_16S_finalOrdered2
taxa_ITS_final_comboOrdered2 <- taxa_ITS_finalOrdered2



ggarrange(
  taxa_16S_final_comboOrdered2,
  taxa_ITS_final_comboOrdered2,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 1,
  nrow = 2)

# save as 15x10 portrait