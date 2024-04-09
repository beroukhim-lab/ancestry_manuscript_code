#Set up the environment
library(tidyverse)
library(reshape2)


#Load in the top snp data
setwd("../data")
top.snp.for.extraction = read.table("top.snp.fdr.txt", sep = "\t", header = T) %>%
select(snp) %>%
data.frame()

#write the bed file into the working directory for this project
setwd("../data")
write.table(top.snp.for.extraction, "top_snp_for_extraction.txt", sep ="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

