#Set up the environment
library(tidyverse)
library(reshape2)



#Load in random guide data
setwd("../data")
random.sites = read.table("final_randoms.txt", sep = "\t", header = T)



#Calculate the number of affected guides in each individual
affected.guides.per.person = random.sites %>%
select(-seq, -chr, -posL, -posR) %>%
apply(2, sum)



#Convert it to a data frame
#Write it, so we don't need to do all of that math again
affected.guides.df = cbind(names(affected.guides.per.person), affected.guides.per.person) %>%
data.frame() %>%
rename("sample" = 1, "num_affected" = 2)

#Now write it
setwd("../data")
write.table(affected.guides.df, "affected_guides_per_gnomad_sample.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)