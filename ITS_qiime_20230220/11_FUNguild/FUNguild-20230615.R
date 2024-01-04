library(dplyr)
library(ggplot2)
library(reshape2)


setwd("/Users/alix22792/Desktop/tmp")

# read in full dataset
funguild<-read.delim("OTU-table-frombiom-with-taxonomy.guilds.txt", header=T)

# filter high confidence assignments - Probable and Highly Probable
funguild.conf <-
    funguild %>%
    filter(Confidence.Ranking %in% c("Probable", "Highly Probable"))

# reshape data
funguild.conf_m<-melt(funguild.conf)

# add category variable
funguild.conf_m <-
    funguild.conf_m %>% mutate(category = substr(variable, 1, 2))

funguild.conf_m$category <- as.factor(funguild.conf_m$category)



#calculate abundance

funguild_sum <-
    funguild.conf_m %>%
    group_by(Trophic.Mode, category) %>%
    summarise(total.reads = sum(value))

sum(funguild_sum$total.reads)
#114412 TOTAL reads



aggregate(funguild_sum$total.reads, list(funguild_sum$category), FUN=sum)
#106302 reads for feathers
#8110 reads for mites


funguild_sum_perc<-
    funguild_sum %>%
    group_by(category) %>%
    mutate(percent = total.reads/sum(total.reads))



ggplot(funguild_sum_perc, aes(category, percent, fill = Trophic.Mode))+
    geom_bar(stat='identity')+
    theme_classic()+
    ylab("Relative Abundance")+
    xlab("")+
    labs(fill="Trophic Mode") +
    scale_fill_brewer(palette="Dark2") +
    scale_x_discrete(labels=c("PF" = "Feather", "PM" = "Mites")) +
    theme(legend.title = element_text(size=16),
          axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size=16),
          legend.text=element_text(size=14))

# save as 6x8 landscape