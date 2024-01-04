### Density-dependent selection 16S data
### Alix Matthews
### 2023 May 22

### R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

# 1. Set working directory and load libraries ####
setwd("~/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/08_R/dietr/20230501/electivity_20230522")

library(dplyr) # 1.0.10
library(ggplot2) # 3.3.6
library(ggthemes) # 4.2.4


# 2. Genus V-S selection ####
## 2.1. Load data ####

# Curation steps of this file: This file consists of the top 10 most abundant ASVs that are available. From this file used for dietR analyses (`16S_grouped-filtered-decontam-bio_samples-silva-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`), I calculated the sum of each column (ASV). I then copied/paste-transposed these values and their associated ASV into a new datasheet. I sorted the output to find the top 10 with the highest overall abundance value. I took those associated ASV identities and search for them in the `genus_vs_long.csv` file I created in the `dietR_16S_20230502.R` script. Then I copied over the associated V-S value for each of those top 10 ASVs into a reduced metadata file which is what is loaded here.

electivity_db_top10avail <- read.csv("16S_electivity_top10avail_20230522.csv", header = TRUE)
str(electivity_db_top10avail)
electivity_db_top10avail[,c(4:5, 7:17)] <- sapply(electivity_db_top10avail[,c(4:5, 7:17)],as.numeric)

# View(electivity_db_top10avail)

# count number of NAs per sample and output the ones with more than 5 NAs (not including 5)
which(rowSums(is.na(electivity_db_top10avail)) > 5 )

# find ID of those samples, out of curiosity
electivity_db_top10avail$sample.id[12] # PM22
electivity_db_top10avail$sample.id[13] # PM23

# remove those outliers
electivity_db_top10avail_nooutliers <- electivity_db_top10avail[-c(12,13),]

## 2.2 Stats ####

### 2.2.1. Check normality of residuals on the models ####

m1 <- lm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail)
err <- resid(m1)
shapiro.test(err) # W = 0.98588, p-value = 0.9935, normal

m2 <- lm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail)
err2 <- resid(m2)
shapiro.test(err2) # W = 0.93758, p-value = 0.3205, normal

m3 <- lm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail_nooutliers)
err3 <- resid(m3)
shapiro.test(err3) # W = 0.98234, p-value = 0.9861, normal

m4 <- lm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail_nooutliers)
err4 <- resid(m4)
shapiro.test(err4) # W = 0.9428, p-value = 0.4554, normal

### 2.2.2. GLMs ####
mod1 <- glm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail)
summary(mod1) # t14 = 1.250, p = 0.232

mod1.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod1.sex.age.int) # all NS, so remove interaction

mod1.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod1.sex.age.add) # mite effect: t11 = 1.242, p = 0.23992



mod2 <- glm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail)
summary(mod2) # t14 = 1.928, p = 0.0744

mod2.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod2.sex.age.int) # all NS, so remove interaction

mod2.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod2.sex.age.add) # mite effect: t11 = 1.519, p = 0.15694



mod3 <- glm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail_nooutliers)
summary(mod3) # t12 = 1.152, p = 0.27166

mod3.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod3.sex.age.int) # all NS, so remove interaction

mod3.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod3.sex.age.int) # mite effect: t9 = 1.077, p = 0.30938



mod4 <- glm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail_nooutliers)
summary(mod4) # t14 = 1.879, p = 0.084734

mod4.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod4.sex.age.int) # all NS, so remove interaction

mod4.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod4.sex.age.add) # mite effect: t9 = 1.341, p = 0.21263


## 2.3 Plotting ####

### 2.3.1. VS_electivity_AVERAGE~num_mites_feather (mod1) ####
ggplot(electivity_db_top10avail, aes(y = VS_electivity_AVERAGE, x = num_mites_feather)) +
  geom_point(size = 3, col = "black") +
  geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE) +
  xlab("Number of Mites (Feather)") +
  ylab("Average Electivity Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold")) +
  geom_hline(yintercept = 0,linetype = "dashed", color="darkred", size = 1)

# save as A6  


### 2.3.2. VS_electivity_AVERAGE~num_mites_total (mod2) ####
ggplot(electivity_db_top10avail, aes(y = VS_electivity_AVERAGE, x = num_mites_total)) +
  geom_point(size = 3, col = "black") +
  geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE) +
  xlab("Number of Mites (Total)") +
  ylab("Average Electivity Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold")) +
  geom_hline(yintercept = 0,linetype = "dashed", color="darkred", size = 1)

# save as A6  


### 2.3.3. VS_electivity_AVERAGE~num_mites_feather, no outliers (mod3) ####
ggplot(electivity_db_top10avail_nooutliers, aes(y = VS_electivity_AVERAGE, x = num_mites_feather)) +
  geom_point(size = 3, col = "black") +
  geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE) +
  xlab("Number of Mites (Feather)") +
  ylab("Average Electivity Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold")) +
  geom_hline(yintercept = 0,linetype = "dashed", color="darkred", size = 1)

# save as A6  



### 2.3.4. VS_electivity_AVERAGE~num_mites_total, no outliers (mod4) ####
ggplot(electivity_db_top10avail_nooutliers, aes(y = VS_electivity_AVERAGE, x = num_mites_total)) +
  geom_point(size = 3, col = "black") +
  geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE) +
  xlab("Number of Mites (Total)") +
  ylab("Average Electivity Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold")) +
  geom_hline(yintercept = 0,linetype = "dashed", color="darkred", size = 1)
