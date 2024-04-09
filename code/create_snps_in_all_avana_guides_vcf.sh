#Filter the ccle variant calls so that it only includes SNPs that map to at least one sgRNA targeting sequence
cd ../data
bcftools view -R ../data/Avana14_filtering.bed split.all_chroms.maf.subset.reheader.vcf.gz -o ../data/snps.in.all.avana.guides.vcf.gz -Oz
