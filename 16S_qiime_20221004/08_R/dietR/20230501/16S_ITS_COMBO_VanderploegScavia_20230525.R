# 16S and ITS V-S electivity by genus figures
# updated with new (corrected) ITS Maaslin2 results
setwd("~/Desktop/Alix")


# 16S PATH: 20220523_16S/AHPCC/01_qiime_20221004/08_R/dietr/20230501
# ITS PATH: 20220621_ITS/AHPCC/01_qiime_20230220/09_R/dietr/20230501

# library(devtools)
# install_github("sborstein/dietr")
library(dietr) # 1.1.4
library(dplyr) # 1.0.10
library(ggplot2) # 3.3.6
library(RColorBrewer) # 1.1-3
library(tidyverse) # 1.3.2
library(tidyr) # 1.2.1





# 16S data
genus_vs_long_plotting_colors_16S <- read.csv("./20220523_16S/AHPCC/01_qiime_20221004/08_R/dietr/20230501/genus_vs_long_plotting_20230522.csv", header = TRUE)

# filter out -1 and samples with no CI (only 1 'available' feather sample... filter out NAs and zeros (aka -1 for VS) for plot (colors by electivity)
genus_vs_long_plotting_colors_nozero_16S <- 
  genus_vs_long_plotting_colors_16S %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != -1))





# ITS data
genus_vs_long_plotting_colors_ITS <- read.csv("./20220621_ITS/AHPCC/01_qiime_20230220/09_R/dietr/20230501/genus_vs_long_plotting_20230522.csv", header = TRUE)

# filter out -1 and samples with no CI (only 1 'available' feather sample... filter out NAs and zeros (aka -1 for VS) for plot (colors by electivity)
genus_vs_long_plotting_colors_nozero_ITS <- 
  genus_vs_long_plotting_colors_ITS %>%
  na.omit() %>%
  filter(if_all(everything(.), ~. != -1))










# 16S plots

p16s <- ggplot(data = genus_vs_long_plotting_colors_nozero_16S, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs, col = electivity_category)) +
  geom_pointrange(fatten = 3) +
  scale_color_manual(values = c("darkred", "darkgreen", "gray50")) +
  scale_x_discrete(labels = c("Corynebacterium", "Mycobacterium", "Jatrophihabitans*", "hgcI_clade", "Candidatus Limnoluna", "Marmoricola", "Nocardioides*", "Actinomycetospora*", "Hymenobacter*", "Flavobacterium*", "Solitalea*", "Allobaculum", "Enterococcus", "Lactobacillus", "Streptococcus", "Roseisolibacter*", "Craurococcus-Caldovatus", "Roseomonas", "Brevundimonas", "Phenylobacterium*", "Methylobacterium-Methylorubrum*", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Bartonella*", "Sphingomonas*", "Limnobacter*", "Polynucleobacter", "Delftia*", "Methylotenera*", "Massilia*", "Escherichia.Shigella", "Pseudomonas", "Nevskia*", "Candidatus Omnitrophus")) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Bacteria") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6),
        legend.position = "none")

p16s




p16stop <- p16s +
  scale_color_manual(values = c("darkred", "skyblue2", "gray65")) +
  ylim(-1.20, 1.20) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 22)
        #axis.title.y=element_blank(),
        # axis.text.x=element_blank()
  )






# ITS plots

its <- ggplot(data = genus_vs_long_plotting_colors_nozero_ITS, aes(x = ASV, y = mean.vs, ymin = lower.ci.vs, ymax = upper.ci.vs, col = electivity_category)) +
  geom_pointrange(fatten = 3) +
  scale_color_manual(values = c("darkred", "darkgreen", "gray50")) +
  scale_x_discrete(labels = c("Botryosphaeria*", "Cladosporium*", "Rachicladosporium*", "Neodevriesia*", "Aureobasidium", "Ramularia*", "Corylicola", "Didymella*", "Neoascochyta", "Nothophoma*", "Kalmusia", "Paracamarosporium", "Paraconiothyrium", "Nigrograna", "Leptospora*", "Neosetophoma", "Sclerostagonospora", "Setophaeosphaeria", "Alternaria*", "Neoroussoella", "Roussoella", "Xenoroussoella*", "Magnibotryascoma", "Cyphellophora*", "Rhinocladiella", "Trichomerium*", "Penicillium", "Sporopachydermia*", "Kluyveromyces", "Candida*", "Pestalotiopsis", "Diaporthe*", "Phomopsis", "Dendrostoma", "Tubakia", "Beauveria", "Trichoderma*", "Sarocladium*", "Fusarium*", "Alfaria", "Thyridium*", "Taphrina*", "Pluteus", "Coprinellus", "Fuscoporia", "Hydnoporia*", "Bjerkandera", "Trametes", "Peniophora", "Stereum*", "Trechispora", "Cystobasidium", "Filobasidium", "Naganishia*")) +
  coord_flip() +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Fungi") +
  ylab("Vanderploeg and Scavia's Relativized Electivity Index (95% CI)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 6),
        legend.position = "none")

its









ylim_its <- its +
  scale_color_manual(values = c("darkred", "skyblue2", "gray65")) +
  theme(axis.title.y=element_text(size = 22)) +
  ylim(-1.20, 1.20)



# alpha - arrange 'em
ggarrange(
  p16stop,
  ylim_its,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 1,
  nrow = 2)        


# export as 8 x 11 portrait