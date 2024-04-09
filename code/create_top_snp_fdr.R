#Set up the environment
library(tidyverse)
library(reshape2)
library(data.table)


setwd("../data")
file.list = list.files(path = ".", pattern = "glm.linear")

#extract the SNP IDs from the dataset
all.data = read.table(file.list[1], sep = "\t") %>%
pull(3) %>%
data.frame()

for(file in file.list){

    #Load in the file and select only the p-value column
    temp.file = read.table(file, sep = "\t") %>%
    pull(12) %>%
    data.frame()
    
    #Add the pvalues to the big file
    all.data = cbind(all.data, temp.file)
}





#Set the names correctly
#First, get a vector of names
col.name.vector = gsub("output.", "", file.list) %>%
gsub(".glm.linear", "", .) %>%
c("snp", .)
colnames(all.data) = col.name.vector #Assign the names
head(all.data)


#Write out the file so that we don't need to re-load and merge everything again
setwd("../data")
write.table(all.data, "merged.pvals.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





#FDR correct everything
all.data.fdr = all.data %>%
select(-snp) %>%
apply(2, function(x){ p.adjust(x, method = "BH", n = length(x)) }) %>%
data.frame() %>%
mutate("snp" = all.data$snp, .before = 1)

#Write the FDR data so that we don't need to do all of those calculations again.
setwd("../data")
write.table(all.data.fdr, "merged.fdr.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





#Assemble a dataframe that has the gene, lowest fdr, and the corresponding SNP
top.snp = all.data.fdr %>%
select(-snp) %>%
apply(2, function(x){
    
    #Get the index of the smallest fdr
    min.index = which.min(x)
    
    #Then get the SNP with that pval
    min.snp = all.data.fdr$snp[min.index]
    
    #Then get the SNP with the smallest fdr
    min.fdr = min(x)
    
    #paste them together and return it outside of the vector
    c(min.snp, min.fdr) %>% return()

}) %>%
data.frame() %>%
t() %>%
data.frame() %>%
rename("snp" = 1, "fdr" = 2) %>%
mutate("gene" = colnames(all.data.fdr[-1]), .before = 1) %>%
mutate(fdr = as.numeric(fdr)) %>%
arrange(fdr)

#Write it, so that if I do need to go back and re-make this figure, I don't need to do a whole bunch of caluclations
setwd("../data")
write.table(top.snp, "top.snp.fdr.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)