### dietR ITS data
### Alix Matthews
### 2023 May 14

### R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

### OS: Ubuntu 20.04.3 LTS

# 1. Set working directory and load libraries ####
setwd("~/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/09_R/dietr/20230501")

# library(devtools)
# install_github("sborstein/dietr")
library(dietr) # 1.1.4
library(dplyr) # 1.0.10
library(ggplot2) # 3.3.6
library(RColorBrewer) # 1.1-3
library(tidyverse) # 1.3.2
library(tidyr) # 1.2.1


# 2. Phyla level ####
## 2.1. Load phyla data ####
phyla_consumed <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv", header = TRUE)
ncol(phyla_consumed)
phyla_consumed[3:10] <- lapply(phyla_consumed[3:10], as.numeric)
str(phyla_consumed)

phyla_available <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-2_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv", header = TRUE)
phyla_available[2:9] <- lapply(phyla_available[2:9], as.numeric)
str(phyla_available)

## 2.2. Calculate phyla electivity ####

phyla.indices <- Electivity(phyla_consumed, phyla_available, Indices = c("ForageRatio", "Ivlev", "Strauss", "JacobsQ", "JacobsD", "Chesson", "VanderploegScavia"), LogQ = TRUE, CalcAbundance = FALSE, Depleting = FALSE)

# plot example one...
PlotElectivity(phyla.indices, Indices = "Chesson")


### 2.2.1. Get means for electivity indices ####

phyla.indices$ForageRatio %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota:k__Fungi.__combo), na.rm = TRUE))

phyla.indices$VanderploegScavia %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota:k__Fungi.__combo), na.rm = TRUE))

phyla.indices$Chesson %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota:k__Fungi.__combo), na.rm = TRUE))



## 2.3. Plotting phyla electivity ####


### 2.3.1. Chesson's ####
phyla_chesson <- as.data.frame(phyla.indices$Chesson)

colnames(phyla_chesson)

phyla_chesson_long <-
  phyla_chesson %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota:k__Fungi.__combo,
    names_to = "ASV",
    values_to = "Chesson"
  )

# View(phyla_chesson_long_plotting)

write.csv(phyla_chesson_long_plotting, file = "phyla_chesson_long_plotting.csv")

phyla_chesson_long_plotting <- 
  phyla_chesson_long %>%
  group_by(ASV) %>%
  summarize(mean.Chesson = mean(Chesson, na.rm = TRUE),
            sd.Chesson = sd(Chesson, na.rm = TRUE),
            n.Chesson = n()) %>%
  mutate(se.Chesson = sd.Chesson / sqrt(n.Chesson),
         lower.ci.Chesson = mean.Chesson - qt(1-(0.05/2), n.Chesson-1) * se.Chesson, 
         upper.ci.Chesson = mean.Chesson + qt(1-(0.05/2), n.Chesson-1) * se.Chesson)


ggplot(data = phyla_chesson_long_plotting, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(phyla_available)-1)), lty = 2) +
  xlab("Funfi (Phyla)") +
  ylab("Chesson's Electivity Index (95% CI)") +
  theme_bw() +
  scale_x_discrete(labels = c("Unidentified Fungi", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Rozellomycota")) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))



# save as landscape 6x4




# filter out NAs and zeros for plot
phyla_chesson_long_plotting_nozero <- 
  phyla_chesson_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != 0))


ggplot(data = phyla_chesson_long_plotting_nozero, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 5, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(phyla_available)-1)), lty = 2) +
  xlab("Fungi (Phyla)") +
  ylab("Chesson Electivity Index (95% CI)") +
  theme_bw() +
  scale_x_discrete(labels = c("Unidentified Fungi", "Ascomycota", "Basidiomycota")) +
  theme(axis.title = element_text(size = 14))









### 2.3.2. Vanderploeg and Scavia's ####
phyla_vs <- as.data.frame(phyla.indices$VanderploegScavia)

phyla_vs_long <-
  phyla_vs %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota:k__Fungi.__combo,
    names_to = "ASV",
    values_to = "VanderploegScavia"
  )

# View(phyla_vs_long)

write.csv(phyla_vs_long, file = "phyla_vs_long.csv")


phyla_vs_long_plotting <- 
  phyla_vs_long %>%
  group_by(ASV) %>%
  summarize(mean.VanderploegScavia = mean(VanderploegScavia, na.rm = TRUE),
            sd.VanderploegScavia = sd(VanderploegScavia, na.rm = TRUE),
            n.VanderploegScavia = n()) %>%
  mutate(se.VanderploegScavia = sd.VanderploegScavia / sqrt(n.VanderploegScavia),
         lower.ci.VanderploegScavia = mean.VanderploegScavia - qt(1-(0.05/2), n.VanderploegScavia-1) * se.VanderploegScavia, 
         upper.ci.VanderploegScavia = mean.VanderploegScavia + qt(1-(0.05/2), n.VanderploegScavia-1) * se.VanderploegScavia)


write.csv(phyla_vs_long_plotting, file = "phyla_vs_long_plotting.csv")



ggplot(data = phyla_vs_long_plotting, aes(x = ASV, y = mean.VanderploegScavia, ymin = lower.ci.VanderploegScavia, ymax = upper.ci.VanderploegScavia)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (Phyla)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  scale_x_discrete(labels = c("Unidentified Fungi", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Rozellomycota")) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# save as landscape 4x7










### 2.3.3. Forage Ratio ####
phyla_forage <- as.data.frame(phyla.indices$ForageRatio)

phyla_forage_long <-
  phyla_forage %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota:k__Fungi.__combo,
    names_to = "ASV",
    values_to = "ForageRatio"
  )

# View(phyla_forage_long)

phyla_forage_long_plotting <- 
  phyla_forage_long %>%
  group_by(ASV) %>%
  summarize(mean.ForageRatio = mean(ForageRatio, na.rm = TRUE),
            sd.ForageRatio = sd(ForageRatio, na.rm = TRUE),
            n.ForageRatio = n()) %>%
  mutate(se.ForageRatio = sd.ForageRatio / sqrt(n.ForageRatio),
         lower.ci.ForageRatio = mean.ForageRatio - qt(1-(0.05/2), n.ForageRatio-1) * se.ForageRatio, 
         upper.ci.ForageRatio = mean.ForageRatio + qt(1-(0.05/2), n.ForageRatio-1) * se.ForageRatio)

write.csv(phyla_forage_long_plotting, file = "phyla_forage_long_plotting.csv")


ggplot(data = phyla_forage_long_plotting, aes(x = ASV, y = mean.ForageRatio, ymin = lower.ci.ForageRatio, ymax = upper.ci.ForageRatio)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Fungi (Phyla)") +
  ylab("Forage Ratio (95% CI)") +
  theme_bw() +
  scale_x_discrete(labels = c("Unidentified Fungi", "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", "Rozellomycota")) +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# save as portrait 6x6







































# 3. Family level ####
## 3.1. Load family data ####
family_consumed <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv", header = TRUE)
ncol(family_consumed)
family_consumed[3:254] <- lapply(family_consumed[3:254], as.numeric)

family_available <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-5_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv", header = TRUE)
ncol(family_available)
family_available[2:253] <- lapply(family_available[2:253], as.numeric)

## 3.2. Calculate family electivity ####

family.indices <- Electivity(family_consumed, family_available, Indices = c("ForageRatio", "Ivlev", "Strauss", "JacobsQ", "JacobsD", "Chesson", "VanderploegScavia"), LogQ = TRUE, CalcAbundance = FALSE, Depleting = FALSE)

# plot example one...
PlotElectivity(family.indices, Indices = "Chesson")


### 3.2.1. Get means for electivity indices ####

colnames(family_available)

family.indices$ForageRatio %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae), na.rm = TRUE))

family.indices$VanderploegScavia %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae), na.rm = TRUE))

family.indices$Chesson %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae), na.rm = TRUE))






## 3.3. Plotting family electivity ####


### 3.3.1. Chesson's ####
family_chesson <- as.data.frame(family.indices$Chesson)

colnames(family_chesson)

family_chesson_long <-
  family_chesson %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae,
    names_to = "ASV",
    values_to = "Chesson"
  )

family_chesson_long_plotting <- 
  family_chesson_long %>%
  group_by(ASV) %>%
  summarize(mean.Chesson = mean(Chesson, na.rm = TRUE),
            sd.Chesson = sd(Chesson, na.rm = TRUE),
            n.Chesson = n()) %>%
  mutate(se.Chesson = sd.Chesson / sqrt(n.Chesson),
         lower.ci.Chesson = mean.Chesson - qt(1-(0.05/2), n.Chesson-1) * se.Chesson, 
         upper.ci.Chesson = mean.Chesson + qt(1-(0.05/2), n.Chesson-1) * se.Chesson)

# View(family_chesson_long_plotting)

write.csv(family_chesson_long_plotting, file = "family_chesson_long_plotting.csv")


ggplot(data = family_chesson_long_plotting, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(family_available)-1)), lty = 2) +
  xlab("Fungi (family)") +
  ylab("Chesson's Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 16x10 portrait


# filter out NAs and zeros for plot
family_chesson_long_plotting_nozero <- 
  family_chesson_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != 0))


ggplot(data = family_chesson_long_plotting_nozero, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(family_available)-1)), lty = 2) +
  xlab("Fungi (family)") +
  ylab("Chesson's Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x10 landscape















### 3.3.2. Vanderploeg and Scavia's ####

family_vs <- as.data.frame(family.indices$VanderploegScavia)

colnames(family_vs)

family_vs_long <-
  family_vs %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae,
    names_to = "ASV",
    values_to = "vs"
  )


write.csv(family_vs_long, file = "family_vs_long.csv")


family_vs_long_plotting <- 
  family_vs_long %>%
  group_by(ASV) %>%
  summarize(mean.vs = mean(vs, na.rm = TRUE),
            sd.vs = sd(vs, na.rm = TRUE),
            n.vs = n()) %>%
  mutate(se.vs = sd.vs / sqrt(n.vs),
         lower.ci.vs = mean.vs - qt(1-(0.05/2), n.vs-1) * se.vs, 
         upper.ci.vs = mean.vs + qt(1-(0.05/2), n.vs-1) * se.vs)

# View(family_vs_long_plotting)

write.csv(family_vs_long_plotting, file = "family_vs_long_plotting.csv")



ggplot(data = family_vs_long_plotting, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (family)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 16x11 portrait


# filter out NAs and zeros (aka -1 for VS) for plot
family_vs_long_plotting_nozero <- 
  family_vs_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != -1))


ggplot(data = family_vs_long_plotting_nozero, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (family)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x10 landscape
























### 3.3.3. Forage Ratio ####
family_forage <- as.data.frame(family.indices$ForageRatio)


colnames(family_forage)

family_forage_long <-
  family_forage %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae:k__Fungi.p__Ascomycota.c__Sareomycetes.o__Sareales.f__Sareaceae,
    names_to = "ASV",
    values_to = "forage"
  )

family_forage_long_plotting <- 
  family_forage_long %>%
  group_by(ASV) %>%
  summarize(mean.forage = mean(forage, na.rm = TRUE),
            sd.forage = sd(forage, na.rm = TRUE),
            n.forage = n()) %>%
  mutate(se.forage = sd.forage / sqrt(n.forage),
         lower.ci.forage = mean.forage - qt(1-(0.05/2), n.forage-1) * se.forage, 
         upper.ci.forage = mean.forage + qt(1-(0.05/2), n.forage-1) * se.forage)

# View(family_forage_long_plotting)

write.csv(family_forage_long_plotting, file = "family_forage_long_plotting.csv")


ggplot(data = family_forage_long_plotting, aes(x = ASV, y = mean.forage, ymin = lower.ci.forage, ymax = upper.ci.forage)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Fungi (family)") +
  ylab("Forage Ratio (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 16x10 portrait


# filter out NAs and zeros for plot
family_forage_long_plotting_nozero <- 
  family_forage_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != 0))


ggplot(data = family_forage_long_plotting_nozero, aes(x = ASV, y = mean.forage, ymin = lower.ci.forage, ymax = upper.ci.forage)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Fungi (family)") +
  ylab("Forage Ratio (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x10 landscape






































# 4. Genus level ####
## 4.1. Load genus data ####
genus_consumed <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_consumed.csv", header = TRUE)
ncol(genus_consumed)
genus_consumed[3:537] <- lapply(genus_consumed[3:537], as.numeric)

genus_available <- read.csv("ITS_forward_afterQtrim-pmin1-dada2_table-fungi-decontam-bio_samples-unite-level-6_nomissing_relabun_curated_MaAsLin2_dietr_paired_available.csv", header = TRUE)
ncol(genus_available)
genus_available[2:536] <- lapply(genus_available[2:536], as.numeric)



## 4.2. Calculate genus electivity ####

genus.indices <- Electivity(genus_consumed, genus_available, Indices = c("ForageRatio", "Ivlev", "Strauss", "JacobsQ", "JacobsD", "Chesson", "VanderploegScavia"), LogQ = TRUE, CalcAbundance = FALSE, Depleting = FALSE)

# plot example one...
PlotElectivity(genus.indices, Indices = "Chesson")


### 4.2.1. Get means for electivity indices ####

colnames(genus_available)

genus.indices$ForageRatio %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella), na.rm = TRUE))

genus.indices$VanderploegScavia %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella), na.rm = TRUE))

genus.indices$Chesson %>%
  rowwise(Record) %>%
  summarise(m = mean(c_across(k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella), na.rm = TRUE))






## 4.3. Plotting genus electivity ####


### 4.3.1. Chesson's ####
genus_chesson <- as.data.frame(genus.indices$Chesson)

colnames(genus_chesson)

genus_chesson_long <-
  genus_chesson %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella,
    names_to = "ASV",
    values_to = "Chesson"
  )


genus_chesson_long_plotting <- 
  genus_chesson_long %>%
  group_by(ASV) %>%
  summarize(mean.Chesson = mean(Chesson, na.rm = TRUE),
            sd.Chesson = sd(Chesson, na.rm = TRUE),
            n.Chesson = n()) %>%
  mutate(se.Chesson = sd.Chesson / sqrt(n.Chesson),
         lower.ci.Chesson = mean.Chesson - qt(1-(0.05/2), n.Chesson-1) * se.Chesson, 
         upper.ci.Chesson = mean.Chesson + qt(1-(0.05/2), n.Chesson-1) * se.Chesson)

# View(genus_chesson_long_plotting)

write.csv(genus_chesson_long_plotting, file = "genus_chesson_long_plotting.csv")


ggplot(data = genus_chesson_long_plotting, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(genus_available)-1)), lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Chesson's Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 26x10 portrait





# filter out NAs and zeros for plot
genus_chesson_long_plotting_nozero <- 
  genus_chesson_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != 0))


ggplot(data = genus_chesson_long_plotting_nozero, aes(x = ASV, y = mean.Chesson, ymin = lower.ci.Chesson, ymax = upper.ci.Chesson)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = (1/(ncol(genus_available)-1)), lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Chesson's Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x10 landscape




### 4.3.2. Vanderploeg and Scavia's ####

genus_vs <- as.data.frame(genus.indices$VanderploegScavia)

colnames(genus_vs)

genus_vs_long <-
  genus_vs %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella,
    names_to = "ASV",
    values_to = "vs"
  )


write.csv(genus_vs_long, file = "genus_vs_long.csv")

genus_vs_long_plotting <- 
  genus_vs_long %>%
  group_by(ASV) %>%
  summarize(mean.vs = mean(vs, na.rm = TRUE),
            sd.vs = sd(vs, na.rm = TRUE),
            n.vs = n()) %>%
  mutate(se.vs = sd.vs / sqrt(n.vs),
         lower.ci.vs = mean.vs - qt(1-(0.05/2), n.vs-1) * se.vs, 
         upper.ci.vs = mean.vs + qt(1-(0.05/2), n.vs-1) * se.vs)

# View(genus_vs_long_plotting)

write.csv(genus_vs_long_plotting, file = "genus_vs_long_plotting.csv")

# Now go through this file and add a color column value by selectivity (manually...)
genus_vs_long_plotting_colors <- read.csv("genus_vs_long_plotting_20230522.csv", header = TRUE)


# all the same color


ggplot(data = genus_vs_long_plotting, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 26x12 portrait




# colored by selectivity

ggplot(data = genus_vs_long_plotting_colors, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs, col = electivity_category)) +
  geom_pointrange(fatten = 3) +
  scale_color_manual(values = c("darkred", "darkgreen", "gray50")) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6),
        legend.position = "none")

# saved 26x12 portrait







# filter out NAs and zeros (aka -1 for VS) for plot
genus_vs_long_plotting_nozero <- 
  genus_vs_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != -1))


ggplot(data = genus_vs_long_plotting_nozero, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x11 landscape





# filter out NAs and zeros (aka -1 for VS) for plot (colors by electivity)
genus_vs_long_plotting_colors_nozero <- 
  genus_vs_long_plotting_colors %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != -1))



ggplot(data = genus_vs_long_plotting_colors_nozero, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs, col = electivity_category)) +
  geom_pointrange(fatten = 3) +
  scale_color_manual(values = c("darkred", "darkgreen", "gray50")) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6),
        legend.position = "none")


# saved as 6x11 landscape




# change labels

ggplot(data = genus_vs_long_plotting_colors_nozero, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs, col = electivity_category)) +
  geom_pointrange(fatten = 3) +
  scale_color_manual(values = c("darkred", "darkgreen", "gray50")) +
  scale_x_discrete(labels = c("Botryosphaeria", "Cladosporium", "Rachicladosporium", "Neodevriesia", "Aureobasidium", "Ramularia", "Corylicola", "Didymella", "Neoascochyta", "Nothophoma", "Kalmusia", "Paracamarosporium", "Paraconiothyrium", "Nigrograna", "Leptospora", "Neosetophoma", "Sclerostagonospora", "Setophaeosphaeria", "Alternaria", "Neoroussoella", "Roussoella", "Xenoroussoella", "Magnibotryascoma", "Cyphellophora", "Rhinocladiella", "Trichomerium", "Penicillium", "Sporopachydermia", "Kluyveromyces", "Candida", "Pestalotiopsis", "Diaporthe", "Phomopsis", "Dendrostoma", "Tubakia", "Beauveria", "Trichoderma", "Sarocladium", "Fusarium", "Alfaria", "Thyridium", "Taphrina", "Pluteus", "Coprinellus", "Fuscoporia", "Hydnoporia", "Bjerkandera", "Trametes", "Peniophora", "Stereum", "Trechispora", "Cystobasidium", "Filobasidium", "Naganishia")) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6),
        legend.position = "none")


# saved as 6x8 landscape









### 4.3.3. Forage Ratio ####
genus_forage <- as.data.frame(genus.indices$ForageRatio)


colnames(genus_forage)

genus_forage_long <-
  genus_forage %>%
  pivot_longer(
    cols = k__Fungi.p__Basidiomycota.c__Cystobasidiomycetes.o__Erythrobasidiales.f__Erythrobasidiaceae.g__Erythrobasidium:k__Fungi.p__Ascomycota.c__Sordariomycetes.o__Hypocreales.f__Hypocreaceae.g__Sphaerostilbella,
    names_to = "ASV",
    values_to = "forage"
  )




genus_forage_long_plotting <- 
  genus_forage_long %>%
  group_by(ASV) %>%
  summarize(mean.forage = mean(forage, na.rm = TRUE),
            sd.forage = sd(forage, na.rm = TRUE),
            n.forage = n()) %>%
  mutate(se.forage = sd.forage / sqrt(n.forage),
         lower.ci.forage = mean.forage - qt(1-(0.05/2), n.forage-1) * se.forage, 
         upper.ci.forage = mean.forage + qt(1-(0.05/2), n.forage-1) * se.forage)

# View(genus_forage_long_plotting)

write.csv(genus_forage_long_plotting, file = "genus_forage_long_plotting.csv")



ggplot(data = genus_forage_long_plotting, aes(x = ASV, y = mean.forage, ymin = lower.ci.forage, ymax = upper.ci.forage)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Forage Ratio (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved 26x10 portrait




# filter out NAs and zeros for plot
genus_forage_long_plotting_nozero <- 
  genus_forage_long_plotting %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != 0))


ggplot(data = genus_forage_long_plotting_nozero, aes(x = ASV, y = mean.forage, ymin = lower.ci.forage, ymax = upper.ci.forage)) +
  geom_pointrange(fatten = 3, col = "darkgreen") +
  coord_flip() +
  geom_hline(yintercept = 1, lty = 2) +
  xlab("Fungi (genus)") +
  ylab("Forage Ratio (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6))

# saved as 6x10 landscape
