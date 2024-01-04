### Density-dependent selection ITS data
### Alix Matthews
### 2023 May 22

### R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

# 1. Set working directory and load libraries ####
setwd("~/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/09_R/dietr/20230501/electivity_20230522")

library(dplyr) # 1.0.10
library(ggplot2) # 3.3.6
library(ggthemes) # 4.2.4


# 2. Genus V-S selection ####
## 2.1. Load data ####

# Curation steps of this file: This file consists of the top 10 most abundant ASVs that are available. From this file used for dietR analyses (`ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv`), I calculated the sum of each column (ASV). I then copied/paste-transposed these values and their associated ASV into a new datasheet. I sorted the output to find the top 10 with the highest overall abundance value. I took those associated ASV identities and search for them in the `genus_vs_long.csv` file I created in the `dietR_ITS_20230514.R` script. Then I copied over the associated V-S value for each of those top 10 ASVs into a reduced metadata file which is what is loaded here.

electivity_db_top10avail <- read.csv("ITS_electivity_top10avail_20230522.csv", header = TRUE)
str(electivity_db_top10avail)
electivity_db_top10avail[,c(4:5, 7:17)] <- sapply(electivity_db_top10avail[,c(4:5, 7:17)],as.numeric)

# View(electivity_db_top10avail)

# count number of NAs per sample and output the ones with more than 5 NAs (not including 5)
which(rowSums(is.na(electivity_db_top10avail)) > 5 )

# none based on this criteria, but P14 seems to be a numerical outlier (super high number of mites on the feather which might be driving some of the significance in age later on...)


# remove the outlier
electivity_db_top10avail$sample.id[7] # P14
electivity_db_top10avail_nooutliers <- electivity_db_top10avail[-7,]


## 2.2 Stats ####

### 2.2.1. Check normality of residuals on the models ####

m1 <- lm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail)
err <- resid(m1)
shapiro.test(err) # W = 0.96677, p-value = 0.6858, normal

m2 <- lm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail)
err2 <- resid(m2)
shapiro.test(err2) # W = 0.95446, p-value = 0.4399, normal

m3 <- lm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail_nooutliers)
err3 <- resid(m3)
shapiro.test(err3) # W = 0.97932, p-value = 0.9341, normal

m4 <- lm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail_nooutliers)
err4 <- resid(m4)
shapiro.test(err4) # W = 0.96055, p-value = 0.5832, normal


### 2.2.2. GLMs ####
mod1 <- glm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail)
summary(mod1) # t18 = 0.500, p = 0.623

mod1.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod1.sex.age.int) # all NS, so remove interaction

mod1.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod1.sex.age.add) # mite effect: t14 = 0.493, p = 0.629421... age SY is significant: t14 = -2.338, p = 0.034718



mod2 <- glm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail)
summary(mod2) # t18 = 0.926, p = 0.367

mod2.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod2.sex.age.int) # all NS, so remove interaction

mod2.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail)
summary(mod2.sex.age.add) # mite effect: t14 = 0.489, p = 0.632155, age SY is significant: t14 = -2.288, p = 0.038219




mod3 <- glm(VS_electivity_AVERAGE~num_mites_feather, data = electivity_db_top10avail_nooutliers)
summary(mod3) # t17 = 0.656, p = 0.521

mod3.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod3.sex.age.int) # all NS, so remove interaction

mod3.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_feather + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod3.sex.age.int) # mite effect: t9 = 1.077, p = 0.30938, age SY is significant: t13 = -2.708, p = 0.01791



mod4 <- glm(VS_electivity_AVERAGE~num_mites_total, data = electivity_db_top10avail_nooutliers)
summary(mod4) # t14 = 0.945, p = 0.358

mod4.sex.age.int <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod4.sex.age.int) # all NS, so remove interaction

mod4.sex.age.add <- glm(VS_electivity_AVERAGE~num_mites_total + bird_sex*bird_age, data = electivity_db_top10avail_nooutliers)
summary(mod4.sex.age.add) # mite effect: t13 = -0.158, p = 0.87711, age SY is significant: t13 = -2.416, p = 0.03114




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















## 2.4 Plotting with age ####


### 2.4.1. VS_electivity_AVERAGE~num_mites_feather (mod1) ####
ggplot(electivity_db_top10avail, aes(y = VS_electivity_AVERAGE, x = num_mites_feather)) +
  geom_point(size = 3, aes(col = bird_age)) +
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
  geom_point(size = 3, aes(col = bird_age)) +
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
  geom_point(size = 3, aes(col = bird_age)) +
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
  geom_point(size = 3, aes(col = bird_age)) +
  geom_smooth(method = glm, colour = "black", alpha = 0.3, se = TRUE) +
  xlab("Number of Mites (Total)") +
  ylab("Average Electivity Index") +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold")) +
  geom_hline(yintercept = 0,linetype = "dashed", color="darkred", size = 1)
