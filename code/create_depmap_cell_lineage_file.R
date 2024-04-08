#In some notebooks I load in the full 22q1_sample_info.csv file, in other I load a dataset with only the two columns of interest.
#It's not clear why I did this, but I did.
#I think the easiest path is to just write a small script to extract the columns of interest and re-make depmap_cell_lineage.csv

#Set up the environment
library(tidyverse)

#Read the file and subset
setwd("../data")
sample.info = read.table("22q1_sample_info.csv", sep = ",", header = T, check.names = FALSE, fill = TRUE, row.names = NULL, quote="\"") %>%
select(DepMap_ID, primary_disease)

#Write the file
setwd("../data")
write.table(sample.info, "depmap_cell_lineage.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


