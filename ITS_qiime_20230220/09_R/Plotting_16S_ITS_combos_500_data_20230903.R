# Combine plots from `Plotting_16S.R` and `Plotting_ITS.R` using 500 cutoff (main figures)


library(ggplot2) # 3.4.2
library(RColorBrewer) # 1.1-3
library(ggthemes) # 4.2.4
library(ggpubr) # 0.6.0
library(ggsignif) # 0.6.4
library(dplyr) # 1.1.2
library(tidyr) # 1.3.0
library(lmerTest) # 3.1-3

setwd("~/Desktop/Alix/20220523_16S/AHPCC/01_qiime_20221004/08_R/")


## 1.1 Read in data and view, subset as needed #####
db_500 <- read.csv(file = "Plotting_16S_ITS_combos_500_data.csv", sep = ",")
# str(db_500)

## 1.2 reformat columns #####

# first take care of the chr->numeric
db_500[,c(6:13)] <- sapply(db_500[,c(6:13)],as.numeric)
# then take care of the chr->factor
db_500[sapply(db_500,is.character)] <- lapply(db_500[sapply(db_500, is.character)], as.factor)

## 1.3 Subsetting #####
### 1.3.1 biological samples ######

biological_db<-subset(db_500, bio_or_control !="control")
# str(biological_db)
levels(biological_db$bio_or_control) # controls still there
biological_db$bio_or_control <- factor(biological_db$bio_or_control) # factor it again
levels(biological_db$bio_or_control) # controls gone












# 16S Shannon
shannon_500_jitter_16s<-
  ggplot(biological_db, aes(y=shannon_500_16s, x = category_broad, fill = category_broad)) +
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
  scale_y_continuous(limits = c(1,8.5),
                     breaks=seq(2,8, by = 2)) +
  ylab("Shannon Index") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

shannon_500_jitter_16s

ggsave("./figures/20230903/shannon_500_jitter_16s.pdf", device=cairo_pdf)





# 16S chao
chao_500_16s_jitter<-
  ggplot(biological_db, aes(y=chao_500_16s, x = category_broad, fill = category_broad)) +
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
  ylab("Chao1 Index") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

chao_500_16s_jitter

ggsave("./figures/20230903/chao_500_16s_jitter.pdf", device=cairo_pdf)



# 16S RPCA 500
rpca500_pc1_pc2_16S<-
  ggplot(biological_db, aes(y=RPCA_pc2_500_16S, x = RPCA_pc1_500_16S, color = category_broad)) +
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
  ylab("PC2 (28.20%)") +
  xlab("PC1 (67.08%)") +
  stat_ellipse(linewidth = 1)

rpca500_pc1_pc2_16S

ggsave("./figures/20230903/rpca500_pc1_pc2_16S.pdf", device=cairo_pdf)





# ITS Shannon
shannon_500_its_jitter<-
  ggplot(biological_db, aes(y=shannon_500_its, x = category_broad, fill = category_broad)) +
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
  scale_y_continuous(limits = c(1,8.5),
                     breaks=seq(2,8, by = 2)) +
  ylab("Shannon Index") +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

shannon_500_its_jitter

ggsave("./figures/20230903/shannon_500_its_jitter.pdf", device=cairo_pdf)






# ITS chao1
chao_500_its_jitter<-
  ggplot(biological_db, aes(y=chao_500_its, x = category_broad, fill = category_broad)) +
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
  ylab("Chao1 Index") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  geom_signif(comparisons = list(c("feather","mites")), 
              map_signif_level=TRUE,
              textsize = 5)

chao_500_its_jitter

ggsave("./figures/20230903/chao_500_its_jitter.pdf", device=cairo_pdf)






# ITS RPCA 500
rpca500_pc1_pc2_its<-
  ggplot(biological_db, aes(y=RPCA_pc2_500_its, x = RPCA_pc1_500_its, color = category_broad)) +
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
  ylab("PC2 (32.42%)") +
  xlab("PC1 (47.16%)") +
  stat_ellipse(linewidth = 1)

rpca500_pc1_pc2_its

ggsave("./figures/20230903/rpca500_pc1_pc2_its.pdf", device=cairo_pdf)






# fix 'em up a little


####  chao1
chao_500_16s_jitter_combo <- chao_500_16s_jitter +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=18))


chao_500_its_jitter_combo <- chao_500_its_jitter +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())


#### shannon

shannon_500_16s_jitter_combo <- shannon_500_jitter_16s +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=18))+
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))




shannon_500_its_jitter_combo <- shannon_500_its_jitter +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank()) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))













#### rpca ####
rpca500_pc1_pc2_16S_combo <- rpca500_pc1_pc2_16S +
  annotate("label", x=-0.25, y = -0.40, label = "PERMANOVA \nF = 7.72; p = 0.001 \n\nPERMDISP \nF = 22.06; p = 0.001", size = 4, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.background = element_blank()) +
  ggtitle("Bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))

rpca500_pc1_pc2_16S_combo




rpca500_pc1_pc2_its_combo <- rpca500_pc1_pc2_its +
  annotate("label", x=0.50, y = -0.2, label = "PERMANOVA \nF = 5.43; p = 0.001 \n\nPERMDISP \nF = 5.74; p = 0.001", size = 4, col = "black", fill = "gray95", alpha = 0.2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_blank()) +
  ggtitle("Fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 26))

rpca500_pc1_pc2_its_combo





# alpha - shannon and chao only
ggarrange(
  shannon_500_16s_jitter_combo,
  shannon_500_its_jitter_combo,
  chao_500_16s_jitter_combo,
  chao_500_its_jitter_combo,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 2,
  nrow = 2)

# save as 10x10 portrait






# RPCA 500 - arrange 'em
ggarrange(
  rpca500_pc1_pc2_16S_combo,
  rpca500_pc1_pc2_its_combo,
  align = "hv",
  labels = "AUTO",
  font.label = list(size = 22),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE, legend = 'bottom')

# save as 5 x 10 landscape


