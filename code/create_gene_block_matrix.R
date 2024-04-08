#Set up the environment
library(tidyverse)
library(reshape2)
library(janitor)


#Load in the hg38 refGene gene annotations. These were downloaded from ucsc genome browser on 6.22.22.
#https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
#Then process them so that we have a chr and position for each gene 
setwd("../data")
refgene = read.table("hg38.refGene.gtf", header = F, sep = "\t") %>%
select(V1, V3, V4, V5, 7, V9) %>%
rename("chr" = V1, "type" = V3, "start" = V4, "end" = V5, "strand" = V7, "gene" = V9) %>%
filter(type %in% "transcript") %>%
mutate("position" = ifelse(strand %in% "+", start, end)) %>%
mutate(gene = gsub(";..*", "", gene)) %>%
mutate(gene = gsub("gene_id ", "", gene)) %>%
select(chr, position, gene)




#Load the RFMix data into one big df.
setwd("../data")
rfmix.file.names = list.files(pattern = "msp.tsv")
rfmix.data = NULL
for(rfmix.file in rfmix.file.names){
    rfmix.data.temp = read.table(rfmix.file, sep = "\t", header = T, comment.char = "^", skip = 1)
    rfmix.data = rbind(rfmix.data, rfmix.data.temp)
}
rfmix.data = rfmix.data %>%
select(X.chm, spos, epos, contains(".1")) 




#Annotate rfmix data so that it each row has an (arbitrary) block ID
block_ids = paste("block", seq(1:nrow(rfmix.data)), sep = "") #build a vector of block IDs
rfmix.data = rfmix.data %>%
mutate(block_id = block_ids, .before = X.chm)
head(rfmix.data)




#Loop through all of the genes in 'refgene' and calculate which block each gene falls in
block.assignment = refgene %>% apply(1, function(x){
    #First filter the rfmix data by chromosome
    rfmix.chr.filtered = rfmix.data %>%
    filter(X.chm %in% x[1])
    #Calculate if the gene of interest is between the start and end position of the rfmix block
    between.vector = rfmix.chr.filtered %>% apply(1, function(y){
        between(y[4], y[5], x[2])
    })
    #Now get the block for each gene
    is.true = between.vector[between.vector == TRUE] %>% length()
    if(is.true >= 1){
        block = rfmix.chr.filtered %>%
        pull(block_id) %>%
        .[is.true]
    } else {
        block = "NA"
    }
    return(block)  
})




#Add the block assignment to the refgene object
refgene = refgene %>%
mutate("block.assignment" = block.assignment)




#Now assemble the matrix
#In this code I will go through each row in 
gene.block.matrix = refgene %>%
filter(block.assignment != "NA") %>%
distinct(gene, .keep_all = TRUE) %>%
apply(1, function(x){
    #Get the block of interest for this loop
    block.of.interest = x[4] #This is the block for this iteration of the loop
    gene.of.interest = x[3] #This is the gene for this interation of the loop
    #Filter the row of rfmix.data that has that block
    rfmix.block = rfmix.data %>%
    filter(block_id %in% block.of.interest)
    #Add the gene of interest
    final.output = c(gene.of.interest, rfmix.block) #Add the rfmix gene to that row
    #Export it from the apply loop
    return(final.output)
}) %>%
map_dfr(as.list) #Convert it all to a matrix




#Now export the gene block matrix so that I can use it in other downstream analysis
setwd("../data")
write.table(gene.block.matrix, "gene.block.matrix.txt", sep = "\t", col.names = T, quote = F)