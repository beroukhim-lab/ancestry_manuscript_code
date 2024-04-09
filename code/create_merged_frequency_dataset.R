#Set up the environment
library(tidyverse)
library(reshape2)
library(data.table)


#Read the data back in if needed
setwd("../data")
top.snp.df = read.table("ancestry_top_snp_df.txt", sep = "\t", header = TRUE)


#Get a vector of all of the significant SNPs
significant.snp.vector = top.snp.df %>%
filter(fdr < 0.05) %>%
pull(snp)


#First, create a bed file so that we can subset the ccle genotyping dataset
ancestry.snp.bed.for.subsetting = top.snp.df %>%
select(-gene, -pval, -fdr) %>%
separate(snp, sep = ":", into = c("chr", "pos", "ref", "alt")) %>%
select(-ref, -alt) %>%
mutate(pos = as.numeric(pos)) %>%
mutate("low_pos" = pos - 10) %>%
mutate("high_pos" = pos + 10) %>%
select(-pos)

#Write the bed file
setwd("../data")
write.table(ancestry.snp.bed.for.subsetting, "ancestry_subsetting.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#Subset the genotyping matrix
system(glue::glue("
cd ../data
bcftools filter ../data/ccle_all_called.vcf.gz -R ../data/ancestry_subsetting.bed -o ../data/ccle.ancestry.snps.vcf.gz -Oz
"))


#Load in the vcf sample names
setwd("../data")
ccle.vcf.sample.names = read.table("ccle.vcf.sample.names.txt", sep = "\t") %>%
pull(V1)


#Load in the genotyping matrix and format it
setwd("../data")
formatted.ancestry.snp.genotyping.matrix = read.table("ccle.ancestry.snps.vcf.gz", sep = "\t") %>%
select(-V1, -V2, -V4, -V5, -V6, -V7, -V8, -V9) %>%
rename("snp" = V3) %>%
filter(snp %in% all_of(top.snp.df$snp))


#Do the final steps of formatting to adjust the genotype call format
final.ancestry.snp.genotyping.matrix = formatted.ancestry.snp.genotyping.matrix %>%
select(-snp) %>%
mutate_all(funs(gsub(":..*", "", .))) %>%
mutate_all(funs(gsub("0\\|0", "0", .))) %>%
mutate_all(funs(gsub("0\\|1", "1", .))) %>%
mutate_all(funs(gsub("1\\|0", "1", .))) %>%
mutate_all(funs(gsub("1\\|1", "2", .))) %>%
mutate("snp" = formatted.ancestry.snp.genotyping.matrix$snp, .before = 1)


#Set the column names
colnames(final.ancestry.snp.genotyping.matrix) = c("snp", ccle.vcf.sample.names)


#Get a list of samples to keep that intersect between the genotyping matrix and the chronos matrix
samples.to.keep = intersect(colnames(final.ancestry.snp.genotyping.matrix), chronos.22q1$sample)


#Filter the dataset so that it only includes the samples that have chronos scores.
final.ancestry.snp.genotyping.matrix = final.ancestry.snp.genotyping.matrix %>%
select(snp, all_of(samples.to.keep))


#Load in the data that lists the ancestry group for all of the cell lines
setwd('../data')
rfmix.data = read.table("collapsed.ancestry.information.txt", sep = "\t", header = TRUE) %>%
mutate("max_ancestry" = ifelse(AFR >= 0.8, "afr", 
                              ifelse(AMR >= 0.8, "amr",
                                    ifelse(EAS >= 0.8, "eas",
                                          ifelse(EUR >= 0.8, "eur",
                                                ifelse(SAS >= 0.8, "sas", "admixed")))))) %>%
filter(sample %in% all_of(samples.to.keep))


#Create vectors that have all of the lines by ancestry
afr.lines = rfmix.data %>% filter(max_ancestry %in% "afr") %>% pull(sample)
amr.lines = rfmix.data %>% filter(max_ancestry %in% "amr") %>% pull(sample)
eas.lines = rfmix.data %>% filter(max_ancestry %in% "eas") %>% pull(sample)
eur.lines = rfmix.data %>% filter(max_ancestry %in% "eur") %>% pull(sample)
sas.lines = rfmix.data %>% filter(max_ancestry %in% "sas") %>% pull(sample)
admixed.lines = rfmix.data %>% filter(max_ancestry %in% "admixed") %>% pull(sample)


#Calculate the frequency of every SNP in the population
afr.frequency = final.ancestry.snp.genotyping.matrix %>% 
select(snp, all_of(afr.lines)) %>% 
melt(id = "snp") %>% 
mutate(value = as.numeric(value)) %>% 
group_by(snp) %>%
summarise("num_alleles_afr" = sum(value))

eas.frequency = final.ancestry.snp.genotyping.matrix %>% 
select(snp, all_of(eas.lines)) %>% 
melt(id = "snp") %>% 
mutate(value = as.numeric(value)) %>% 
group_by(snp) %>%
summarise("num_alleles_eas" = sum(value))

eur.frequency = final.ancestry.snp.genotyping.matrix %>% 
select(snp, all_of(eur.lines)) %>% 
melt(id = "snp") %>% 
mutate(value = as.numeric(value)) %>% 
group_by(snp) %>%
summarise("num_alleles_eur" = sum(value))

admixed.frequency = final.ancestry.snp.genotyping.matrix %>% 
select(snp, all_of(admixed.lines)) %>% 
melt(id = "snp") %>% 
mutate(value = as.numeric(value)) %>% 
group_by(snp) %>%
summarise("num_alleles_admixed" = sum(value))

total.frequency = final.ancestry.snp.genotyping.matrix %>% 
melt(id = "snp") %>% 
mutate(value = as.numeric(value)) %>% 
group_by(snp) %>%
summarise("num_alleles_total" = sum(value))


#Group everything into a single dataset
merged.frequency.dataset = afr.frequency %>%
inner_join(eas.frequency, by = "snp") %>%
inner_join(eur.frequency, by = "snp") %>%
inner_join(admixed.frequency, by = "snp") %>%
inner_join(total.frequency, by = "snp") %>%
mutate("afr_allele_fraction" = num_alleles_afr/(2*length(afr.lines))) %>%
mutate("eas_allele_fraction" = num_alleles_eas/(2*length(eas.lines))) %>%
mutate("eur_allele_fraction" = num_alleles_eur/(2*length(eur.lines))) %>%
mutate("admixed_allele_fraction" = num_alleles_admixed/(2*length(admixed.lines))) %>%
mutate("total_allele_fraction" = num_alleles_total/(2*nrow(rfmix.data)))


#Export the dataset
setwd('../data')
write.table(merged.frequency.dataset, "merged_frequency_dataset.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)