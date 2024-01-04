# With rarefaction value set at 1451...
# Combine plots from `Plotting_16S.R` and `Plotting_ITS.R`

setwd("~/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/08_R/")


## 1.1 Read in data and view, subset as needed #####
rarefied_db <- read.csv(file = "Plotting_16S_ITS_combos_rarefied_data.csv", sep = ",")
# str(rarefied_db)

## 1.2 reformat columns #####

# first take care of the chr->numeric
rarefied_db[,c(6:17)] <- sapply(rarefied_db[,c(6:17)],as.numeric)
# then take care of the chr->factor
rarefied_db[sapply(rarefied_db,is.character)] <- lapply(rarefied_db[sapply(rarefied_db, is.character)], as.factor)

## 1.3 Subsetting #####
### 1.3.1 biological samples ######

biological_db<-subset(rarefied_db, bio_or_control !="control")
# str(biological_db)
levels(biological_db$bio_or_control) # controls still there
biological_db$bio_or_control <- factor(biological_db$bio_or_control) # factor it again
levels(biological_db$bio_or_control) # controls gone












# 16S Shannon
shannon1451_jitter<-
  ggplot(biological_db, aes(y=shannon_1451, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Shannon Index") +
  #ggtitle("Sampling read depth = 1451") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

shannon1451_jitter

ggsave("./figures/20230627/shannon1451_jitter.pdf", device=cairo_pdf)





# 16S features
features1451_jitter<-
  ggplot(biological_db, aes(y=observed_features_1451, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values=c("goldenrod1", "darkslategray4")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  ylab("ASV Richness") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  #ggtitle("Sampling read depth = 1451") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

features1451_jitter

ggsave("./figures/20230627/features1451_jitter.pdf", device=cairo_pdf)



# 16S RPCA 0
rpca0_pc1_pc2_16S<-
  ggplot(biological_db, aes(y=RPCA_pc2_0_16S, x = RPCA_pc1_0_16S, color = category_broad)) +
  geom_point(size = 3.5, alpha = 0.75) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=c("goldenrod1", "darkslategray4"), name = "Category", labels=c("Feather", "Mites")) +
  ylab("PC2 (26.96%)") +
  xlab("PC1 (66.05%)") +
  #ggtitle("Sampling read depth = 0\nRPCA") +
  stat_ellipse(size = 1)

rpca0_pc1_pc2_16S

ggsave("./figures/20230627/rpca0_pc1_pc2_16S.pdf", device=cairo_pdf)



# 16S RPCA 1451
rpca1451_pc1_pc2_16S<-
  ggplot(biological_db, aes(y=RPCA_pc2_1451, x = RPCA_pc1_1451, color = category_broad)) +
  geom_point(size = 3.5, alpha = 0.75) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=c("goldenrod1", "darkslategray4"), name = "Category", labels=c("Feather", "Mites")) +
  ylab("PC2 (29.85%)") +
  xlab("PC1 (66.76%)") +
  #ggtitle("Sampling read depth = 1451\nRPCA") +
  stat_ellipse(size = 1)

rpca1451_pc1_pc2_16S

ggsave("./figures/20230627/rpca1451_pc1_pc2_16S.pdf", device=cairo_pdf)





# ITS Shannon
shannon995_jitter<-
  ggplot(biological_db, aes(y=shannon_995, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14)) +
  scale_fill_manual(values=c("wheat3", "midnightblue")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("Shannon Diversity") +
  #ggtitle("Sampling read depth = 995") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

shannon995_jitter

ggsave("./figures/20230627/shannon995_jitter.pdf", device=cairo_pdf)






# ITS features
features995_jitter<-
  ggplot(biological_db, aes(y=observed_features_995, x = category_broad, fill = category_broad)) +
  geom_boxplot(outlier.shape=NA, aes(alpha = 0.6)) +
  geom_line(aes(group = bird_id), alpha = 0.6, linetype = "longdash", position = position_dodge(0.3)) +
  geom_point(aes(fill = category_broad, group = bird_id), size = 2, shape = 21, position = position_dodge(0.3)) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14)) +
  scale_fill_manual(values=c("wheat3", "midnightblue")) +
  scale_x_discrete(labels = c("Feather", "Mites")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  ylab("ASV Richness") +
  #ggtitle("Sampling read depth = 995") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

features995_jitter

ggsave("./figures/20230627/features995_jitter.pdf", device=cairo_pdf)






# ITS RPCA 0
rpca0_pc1_pc2_ITS<-
  ggplot(biological_db, aes(y=RPCA_pc2_0_ITS, x = RPCA_pc1_0_ITS, color = category_broad)) +
  geom_point(size = 3.5, alpha = 0.75) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=c("wheat3", "midnightblue"), name = "Category", labels=c("Feather", "Mites")) +
  ylab("PC2 (36.64%)") +
  xlab("PC1 (45.16%)") +
  #ggtitle("Sampling read depth = 0\nRPCA") +
  stat_ellipse(size = 1)

rpca0_pc1_pc2_ITS

ggsave("./figures/20230627/rpca0_pc1_pc2_ITS.pdf", device=cairo_pdf)



# ITS RPCA 995
rpca995_pc1_pc2_ITS<-
  ggplot(biological_db, aes(y=RPCA_pc2_995, x = RPCA_pc1_995, color = category_broad)) +
  geom_point(size = 3.5, alpha = 0.75) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values=c("wheat3", "midnightblue"), name = "Category", labels=c("Feather", "Mites")) +
  ylab("PC2 (28.34%)") +
  xlab("PC1 (52.85%)") +
  #ggtitle("Sampling read depth = 995\nRPCA") +
  stat_ellipse(size = 1)

rpca995_pc1_pc2_ITS

ggsave("./figures/20230627/rpca995_pc1_pc2.pdf", device=cairo_pdf)









# fix 'em up a little


####  observed features
features1451_jitter_combo <- features1451_jitter +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=18))


features995_jitter_combo <- features995_jitter +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())


#### shannon

shannon1451_jitter_combo <- shannon1451_jitter +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=18))+
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))




shannon995_jitter_combo <- shannon995_jitter +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))













#### rpca ####
rpca1451_pc1_pc2_16S_combo <- rpca1451_pc1_pc2_16S +
  annotate("label", x=-0.28, y = -0.45, label = "PERMANOVA \nF = 6.98; p = 0.001 \n\nPERMDISP \nF = 20.56; p = 0.001", size = 3, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_blank())


rpca0_pc1_pc2_16S_combo <- rpca0_pc1_pc2_16S +
  annotate("label", x=-0.31, y = 0.35, label = "PERMANOVA \nF = 8.29; p = 0.001 \n\nPERMDISP \nF = 11.57; p = 0.001", size = 3, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_blank())




rpca995_pc1_pc2_ITS_combo <- rpca995_pc1_pc2_ITS +
  annotate("label", x=0.50, y = -0.25, label = "PERMANOVA \nF = 1.92; p = 0.06 \n\nPERMDISP \nF = 3.18; p = 0.06", size = 3, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_blank())


rpca0_pc1_pc2_ITS_combo <- rpca0_pc1_pc2_ITS +
  annotate("label", x=0.50, y = 0.20, label = "PERMANOVA \nF = 5.64; p = 0.001 \n\nPERMDISP \nF = 10.47; p = 0.001", size = 3, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_blank())




# alpha - shannon and obs fts only
ggarrange(
  shannon1451_jitter_combo,
  shannon995_jitter_combo,
  features1451_jitter_combo,
  features995_jitter_combo,
  align = "hv",
  labels = "AUTO",
  ncol = 2,
  nrow = 2)

# save as 10x10 portrait







# save as 8 x 10 portrait

# RPCA 0 - arrange 'em
ggarrange(
  rpca0_pc1_pc2_16S_combo,
  rpca0_pc1_pc2_ITS_combo,
  align = "hv",
  labels = "AUTO",
  ncol = 2,
  nrow = 1)

# save as 5 x 10 landscape



# RPCA rarefaction level - arrange 'em
ggarrange(
  rpca1451_pc1_pc2_16S_combo,
  rpca995_pc1_pc2_ITS_combo,
  align = "hv",
  labels = "AUTO",
  ncol = 2,
  nrow = 1)

# save as 5 x 10 landscape
