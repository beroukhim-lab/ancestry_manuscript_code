#Recode the chromosome names in the reference panel from "8" to "chr8"
cd ../data
for i in {1..22};
do
bcftools annotate ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
--rename-chrs ../data/hg19_to_hg38_chr_bridge \
-Oz \
-o 1kg_chr${i}.recodeChr.vcf.gz

tabix -p vcf 1kg_chr${i}.recodeChr.vcf.gz
done

#Also do it for chrX
bcftools annotate ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
--rename-chrs ../data/hg19_to_hg38_chr_bridge \
-Oz \
-o 1kg_chrX.recodeChr.vcf.gz

tabix -p vcf 1kg_chrX.recodeChr.vcf.gz





#Remove unwanted samples from the reference panel.
#The unwanted samples are related, so they are probably going to skew our data.
cd ../data
awk '{print $1}' 1kg_sample_map > 1kg.samples.to.keep


#Now filter all of the bcf files so that they only contain these samples
cd ../data
for i in {1..22};
do
bcftools view -S 1kg.samples.to.keep --force-samples 1kg_chr${i}.recodeChr.vcf.gz -o 1kg.ogsamples.chr${i}.vcf.gz -Oz
done

#Also do chrX
bcftools view -S 1kg.samples.to.keep --force-samples 1kg_chrX.recodeChr.vcf.gz -o 1kg.ogsamples.chrX.vcf.gz -Oz



 


#Split the vcf file into each chromosome
cd ../data

for num in {1..22};
do
bcftools view ccle_all_called.vcf.gz --regions chr${num} -Oz -o chr${num}.split.all_chroms.maf.subset.reheader.vcf.gz
done

bcftools view ccle_all_called.vcf.gz --regions chrX -Oz -o chrX.split.all_chroms.maf.subset.reheader.vcf.gz





#Run RFMix
for i in {1..22};
do
rfmix -f chr${i}.split.all_chroms.maf.subset.reheader.vcf.gz \
-r ../data/1kg_chr${i}.recodeChr.vcf.gz \
-m ../data/1kg_sample_map \
-g ../data/chr${i}_genetic_map.txt \
-o chr${i}.rfmix.output \
--chromosome=chr${i}
done

#Also do chrX
rfmix -f chrX.split.all_chroms.maf.subset.reheader.vcf.gz \
-r ../data/1kg_chrX.recodeChr.vcf.gz \
-m ../data/1kg_sample_map \
-g ../data/chrX_genetic_map.txt \
-o chrX.rfmix.output \
--chromosome=chrX
