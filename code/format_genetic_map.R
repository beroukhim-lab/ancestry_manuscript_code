#Set up the environment
library(tidyverse)
library(reshape2)

#Format the genetic map
#Load in the data and format it nicely
setwd("../data")
genetic.map = read.table("genetic_map_hg38_withX.txt", sep = " ", header = T) %>%
rename("chr" = 1, "pos" = 2, "combined" = 3, "genetic_pos" = 4) %>%
select(chr, pos, genetic_pos) %>%
mutate(chr = paste("chr", chr, sep = ""))

#Get a list of all of the chromosomes
unique.chromosomes = genetic.map %>% pull(chr) %>% unique()

#Loop through all of the chromosomes, format the data and then write the output file.
for(chrom in unique.chromosomes){
    
    isolated.chromosome = genetic.map %>%
    filter(chr %in% chrom) %>%
    mutate(genetic_pos = format(genetic_pos, scientific=F)) %>%
    mutate(genetic_pos = as.numeric(genetic_pos)) %>%
    arrange(genetic_pos)
    
    file.name = paste(chrom, "_genetic_map.txt", sep = "")
    
    write.table(isolated.chromosome, file.name, sep = "\t", col.names = T, row.names = F, quote = F)

}

