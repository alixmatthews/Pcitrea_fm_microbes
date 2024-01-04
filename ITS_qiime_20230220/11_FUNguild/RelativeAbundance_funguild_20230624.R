#### Calculate relative abundance of FUNguild (animal pathogens only) for import for MaAsLin2 (post FUNguild analyses)
#### Alix Matthews
#### 24 June 2023

#### This .txt file was taken from FUNguild output `OTU-table-frombiom-with-taxonomy.guilds.txt`, selecting only animal pathogens (guild column) with 'highly probable' or 'probable' confidence (another column), and then reformatted by transposing the rows/columns

setwd("/home/lse305/Desktop/Alix/20220621_ITS/AHPCC/01_qiime_20230220/11_FUNguild/")

library(dplyr) # v. 1.1.10


# 1. Full dataset ####
## 1.1. Load data and make adjustments to data structure as tibble ####
data <- read.csv(file = "OTU-table-frombiom-with-taxonomy.guilds-animalpathogens_transposed.txt", sep = "\t", header= TRUE)
str(data)


## 1.2. Calculate total reads within a sample ####
colnames(data[2]) # get first OTU name, "X6e5c88724ce4ad5457ab9e3e47170818"
colnames(data[ncol(data)]) # get last pathway name, "b333f8fa0261da8011ef5b9884e27d60"

data <-
  data %>%
  rowwise() %>%
  mutate(raw_total = sum(c_across(X6e5c88724ce4ad5457ab9e3e47170818:b333f8fa0261da8011ef5b9884e27d60), na.rm = T))

# need to ungroup after using rowwise()
data <-
  data %>%
  ungroup()



## 1.3. Calculate relative abundances to new df ####
data_rel_abund <-
  data %>%
  mutate_at(vars(X6e5c88724ce4ad5457ab9e3e47170818:b333f8fa0261da8011ef5b9884e27d60) , funs(relabun = ./data$raw_total))

str(data_rel_abund)

# make sure the summation of all relabun columns sum to 1
data_rel_abund_should_equal_1 <-
  data_rel_abund %>%
  rowwise() %>%
  mutate(relabun_sum = sum(across(ends_with("_relabun")), na.rm = T))

data_rel_abund_should_equal_1$relabun_sum

# all that should =1 do = 1, so things look good!


## 1.4.  Now drop the raw counts and only keep relative abundance ####
grep("relabun", colnames(data_rel_abund)) # column indices that have relabun in them [28:92]

data_rel_abund <- data_rel_abund[-c(2:47)] # remove everything else (but keep SampleID, which is column 1)
str(data_rel_abund)
colnames(data_rel_abund)<-gsub("_relabun","",colnames(data_rel_abund))

write.table(data_rel_abund , file = "OTU-table-frombiom-with-taxonomy.guilds-animalpathogens_transposed_relabun.csv", sep=",", row.names=FALSE)



